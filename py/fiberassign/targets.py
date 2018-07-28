# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.targets
=====================

Functions for loading the target list

"""
from __future__ import absolute_import, division, print_function

import numpy as np

import fitsio

from desitarget.targetmask import desi_mask

from .utils import Logger, Timer

from ._internal import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                        TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                        Target, Targets, TargetTree, TargetsAvailable,
                        FibersAvailable)


def str_to_target_type(input):
    if input == "science":
        return TARGET_TYPE_SCIENCE
    elif input == "sky":
        return TARGET_TYPE_SKY
    elif input == "standard":
        return TARGET_TYPE_STANDARD
    elif input == "safe":
        return TARGET_TYPE_SAFE
    else:
        raise ValueError("unknown target type '{}'".format(input))
    return None


def load_target_file(tgs, tfile, typeforce=None, typecol="DESI_TARGET",
                     sciencemask=None, stdmask=None, skymask=None,
                     safemask=None, rowbuffer=100000):
    tm = Timer()
    tm.start()

    log = Logger()

    validtypes = [
        TARGET_TYPE_SCIENCE,
        TARGET_TYPE_SKY,
        TARGET_TYPE_STANDARD,
        TARGET_TYPE_SAFE
    ]
    if typeforce is not None:
        if typeforce not in validtypes:
            raise RuntimeError("Cannot force objects to be an invalid type")

    if sciencemask is None:
        sciencemask = 0
        sciencemask |= desi_mask["LRG"].mask
        sciencemask |= desi_mask["ELG"].mask
        sciencemask |= desi_mask["QSO"].mask
        # Note: BAD_SKY are treated as science targets with priority == 0
        sciencemask |= desi_mask["BAD_SKY"].mask
        sciencemask |= desi_mask["BGS_ANY"].mask
        sciencemask |= desi_mask["MWS_ANY"].mask
        sciencemask |= desi_mask["SECONDARY_ANY"].mask

    if stdmask is None:
        stdmask = 0
        stdmask |= desi_mask["STD_FAINT"].mask
        stdmask |= desi_mask["STD_WD"].mask
        stdmask |= desi_mask["STD_BRIGHT"].mask

    if skymask is None:
        skymask = 0
        skymask |= desi_mask["SKY"].mask

    if safemask is None:
        safemask = 0

    # Open file
    fits = fitsio.FITS(tfile, mode="r")

    # Total number of rows
    nrows = fits[1].get_nrows()

    log.info("Target file {} has {} rows.  Reading in chunks of {}"
             .format(tfile, nrows, rowbuffer))

    # Create buffers for column data
    d_obscond = np.empty(rowbuffer, dtype=np.int32)
    d_targetid = np.empty(rowbuffer, dtype=np.int64)
    d_ra = np.empty(rowbuffer, dtype=np.float64)
    d_dec = np.empty(rowbuffer, dtype=np.float64)
    d_type = np.empty(rowbuffer, dtype=np.uint8)
    d_nobs = np.empty(rowbuffer, dtype=np.int32)
    d_prior = np.empty(rowbuffer, dtype=np.int32)
    d_subprior = np.empty(rowbuffer, dtype=np.float64)

    offset = 0
    n = rowbuffer
    while offset < nrows:
        if offset + n > nrows:
            n = nrows - offset
        data = fits[1].read(rows=np.arange(offset, offset+n, dtype=np.int64))
        log.debug("Target file {} read rows {} - {}"
                  .format(tfile, offset, offset+n-1))
        d_targetid[0:n] = data["TARGETID"]
        d_ra[0:n] = data["RA"]
        d_dec[0:n] = data["DEC"]
        if typeforce is not None:
            d_type[0:n] = typeforce
        else:
            d_type[0:n] = 0
            d_type[np.where((data[typecol] & sciencemask) != 0)] \
                |= TARGET_TYPE_SCIENCE
            d_type[np.where((data[typecol] & stdmask) != 0)] \
                |= TARGET_TYPE_STANDARD
            d_type[np.where((data[typecol] & skymask) != 0)] \
                |= TARGET_TYPE_SKY
            d_type[np.where((data[typecol] & safemask) != 0)] \
                |= TARGET_TYPE_SAFE

        if "OBSCONDITIONS" in data.dtype.fields:
            d_obscond[0:n] = data["OBSCONDITIONS"]
        else:
            # Set obs conditions mask to be all bits
            d_obscond[0:n] = np.invert(np.zeros(n, dtype=np.int32))

        if "NUMOBS_MORE" in data.dtype.fields:
            d_nobs[0:n] = data["NUMOBS_MORE"]
        else:
            d_nobs[0:n] = np.zeros(n, dtype=np.int32)

        if "PRIORITY" in data.dtype.fields:
            d_prior[0:n] = data["PRIORITY"]
        else:
            d_prior[0:n] = np.zeros(n, dtype=np.int32)

        if "SUBPRIORITY" in data.dtype.fields:
            d_subprior[0:n] = data["SUBPRIORITY"]
        else:
            d_subprior[0:n] = np.zeros(n, dtype=np.float64)

        # Append the data to our targets list.  This will print a
        # warning if there are duplicate target IDs.
        tgs.append(d_targetid[0:n], d_ra[0:n], d_dec[0:n], d_nobs[0:n],
                   d_prior[0:n], d_subprior[0:n], d_obscond[0:n], d_type[0:n])

        offset += n

    tm.stop()
    tm.report("Read target file {}".format(tfile))

    return
