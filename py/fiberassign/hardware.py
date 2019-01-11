# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.hardware
=======================

Functions for loading information about the telescope hardware.

"""
from __future__ import absolute_import, division, print_function

import numpy as np

import fitsio
# from astropy.table import Table

import desimodel.io as dmio

from .utils import Logger

from ._internal import (Hardware, FIBER_STATE_OK, FIBER_STATE_STUCK,
                        FIBER_STATE_BROKEN)


def load_hardware(fiberpos_file=None, gfa_file=None, rundate=None,
                  status_file=None):
    """Create a hardware class representing properties of the telescope.

    Args:
        fiberpos_file (str):  Optional path to the fiber positioner FITS file.
            If not specified, desimodel is used to get the location of the
            default file.
        gfa_file (str):  Optional path to the GFA file.
        rundate ():  XXXX format time stamp
        status_file (str):  Path to fiber status file.  If not specified, all
            fibers are assumed good.

    Returns:
        (Hardware):  The hardware object.

    """
    log = Logger.get()

    # Read the fiber positioner data
    if fiberpos_file is None:
        fiberpos_file = dmio.findfile('focalplane/fiberpos-all.fits')
    log.info("Reading fiber positions from {}".format(fiberpos_file))

    fpdata = fitsio.read(fiberpos_file, ext=1)
    nfiber = len(fpdata)
    log.debug("  fiber position table has {} rows".format(nfiber))

    device_type = np.empty(nfiber, dtype="a8")
    device_type[:] = fpdata["DEVICE_TYPE"]

    # For non-science positioners, no fiber ID is assigned in the positioner
    # file.  This is a pain, since all our quantities are indexed by fiber ID.
    # Instead, we give every position a FIBER value which is unique and which
    # is sorted by LOCATION.  We start these fake FIBER values at 5000.
    fiber = np.copy(fpdata["FIBER"])
    location = np.copy(fpdata["LOCATION"])
    missing_fiber = [x for x, y in enumerate(fiber) if y < 0]
    if len(missing_fiber) > 0:
        missing_locsorted = np.sort(location[missing_fiber])
        locfiber = {y: x for x, y in enumerate(missing_locsorted)}
        fiber[missing_fiber] = [(5000 + locfiber[x]) for x in
                                location[missing_fiber]]

    # Read the status file...
    status = np.empty(nfiber, dtype=np.int32)
    status[:] = FIBER_STATE_OK

    hw = Hardware(fiber,
                  fpdata["PETAL"],
                  fpdata["SPECTRO"],
                  location,
                  fpdata["SLIT"],
                  fpdata["SLITBLOCK"],
                  fpdata["BLOCKFIBER"],
                  fpdata["DEVICE"],
                  device_type,
                  fpdata["X"],
                  fpdata["Y"],
                  fpdata["Z"],
                  fpdata["Q"],
                  fpdata["S"],
                  status)
    return hw
