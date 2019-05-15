# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.hardware
=======================

Functions for loading information about the telescope hardware.

"""
from __future__ import absolute_import, division, print_function

from datetime import datetime

import numpy as np

import fitsio

from astropy.table import Table

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
        rundate (str):  ISO 8601 format time stamp as a string in the
            format YYYY-MM-DDTHH:MM:SS.  If None, uses current time.
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
    pos_rows = np.where(fpdata["DEVICE_TYPE"].astype(str) == "POS")[0]
    etc_rows = np.where(fpdata["DEVICE_TYPE"].astype(str) == "ETC")[0]
    keep_rows = np.unique(np.concatenate((pos_rows, etc_rows)))

    nloc = len(keep_rows)
    log.debug("  fiber position table keeping {} rows".format(nloc))

    device_type = np.full(nloc, "OOPSBUG", dtype="a8")
    device_type[:] = fpdata["DEVICE_TYPE"][keep_rows]

    location = np.copy(fpdata["LOCATION"][keep_rows])
    fiber = np.copy(fpdata["FIBER"][keep_rows])

    # Read the status file...
    status = np.full(nloc, FIBER_STATE_OK, dtype=np.int32)

    if status_file is not None:
        runtime = None
        if rundate is None:
            runtime = datetime.utcnow()
        else:
            runtime = datetime.strptime(rundate, "%Y-%m-%dT%H:%M:%S")
        locindx = {y: x for x, y in enumerate(location)}
        statdata = Table.read(status_file, format="ascii.ecsv")
        for row in statdata:
            loc = row["LOCATION"]
            broken = row["BROKEN"]
            stuck = row["STUCK"]
            start = datetime.strptime(row["START_DATE"],
                                      "%Y-%m-%dT%H:%M:%S")
            stop = datetime.strptime(row["END_DATE"],
                                     "%Y-%m-%dT%H:%M:%S")
            # Check if this row applies to our current run
            if (runtime >= start) and (runtime < stop):
                # yep...
                if broken > 0:
                    status[locindx[loc]] |= FIBER_STATE_BROKEN
                if stuck > 0:
                    status[locindx[loc]] |= FIBER_STATE_STUCK

    hw = Hardware(location,
                  fpdata["PETAL"][keep_rows],
                  fpdata["DEVICE"][keep_rows],
                  fpdata["SLITBLOCK"][keep_rows],
                  fpdata["BLOCKFIBER"][keep_rows],
                  fpdata["SPECTRO"][keep_rows],
                  fiber,
                  fpdata["SLIT"][keep_rows],
                  device_type,
                  fpdata["X"][keep_rows],
                  fpdata["Y"][keep_rows],
                  fpdata["Z"][keep_rows],
                  fpdata["Q"][keep_rows],
                  fpdata["S"][keep_rows],
                  status)
    return hw
