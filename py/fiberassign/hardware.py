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
    log.debug("  fiber position table has {} rows"
              .format(len(fpdata["DEVICE_TYPE"])))
    fprows = np.where(fpdata["DEVICE_TYPE"] == b"POS")[0]
    log.debug("  fiber position table has {} rows with DEVICE_TYPE==POS"
              .format(len(fprows)))
    nfiber = len(fpdata[fprows])

    # Read the status file...
    status = np.empty(nfiber, dtype=np.int32)
    status[:] = FIBER_STATE_OK

    hw = Hardware(fpdata["FIBER"][fprows],
                  fpdata["PETAL"][fprows],
                  fpdata["SPECTRO"][fprows],
                  fpdata["LOCATION"][fprows],
                  fpdata["SLIT"][fprows],
                  fpdata["SLITBLOCK"][fprows],
                  fpdata["BLOCKFIBER"][fprows],
                  fpdata["DEVICE"][fprows],
                  fpdata["X"][fprows],
                  fpdata["Y"][fprows],
                  fpdata["Z"][fprows],
                  fpdata["Q"][fprows],
                  fpdata["S"][fprows],
                  status)
    return hw
