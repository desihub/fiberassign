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

    Currently this just parses the fiber status, but could also read files
    describing other properties of the fiber positioners, etc.

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

    hw = Hardware(fpdata["FIBER"][fprows], fpdata["PETAL"][fprows],
                  fpdata["SPECTRO"][fprows], fpdata["X"][fprows],
                  fpdata["Y"][fprows], fpdata["Z"][fprows], status)

    return hw
