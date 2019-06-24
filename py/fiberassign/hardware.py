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

import desimodel.io as dmio

from .utils import Logger

from ._internal import (Hardware, FIBER_STATE_OK, FIBER_STATE_STUCK,
                        FIBER_STATE_BROKEN, FIBER_STATE_UNASSIGNED,
                        Circle, Segments, Shape)


def load_hardware(focalplane=None, rundate=None):
    """Create a hardware class representing properties of the telescope.

    Args:
        focalplane (tuple):  Override the focalplane model.  If not None, this
            should be a tuple of the same data types returned by
            desimodel.io.load_focalplane()
        rundate (str):  ISO 8601 format time stamp as a string in the
            format YYYY-MM-DDTHH:MM:SS.  If None, uses current time.

    Returns:
        (Hardware):  The hardware object.

    """
    log = Logger.get()

    # The timestamp for this run.
    runtime = None
    if rundate is None:
        runtime = datetime.utcnow()
    else:
        runtime = datetime.strptime(rundate, "%Y-%m-%dT%H:%M:%S")

    # Get the focalplane information
    fp = None
    exclude = None
    state = None
    tmstr = None
    if focalplane is None:
        fp, exclude, state, tmstr = dmio.load_focalplane(runtime)
    else:
        fp, exclude, state = focalplane

    # We are only going to keep rows for LOCATIONs that are assigned to a
    # science or sky monitor positioner.

    log.info("Loaded focalplane for time stamp {}".format(runtime))

    pos_rows = np.where(fp["DEVICE_TYPE"].astype(str) == "POS")[0]
    etc_rows = np.where(fp["DEVICE_TYPE"].astype(str) == "ETC")[0]
    keep_rows = np.unique(np.concatenate((pos_rows, etc_rows)))

    nloc = len(keep_rows)
    log.debug(
        "  focalplane table keeping {} rows for POS and ETC devices"
        .format(nloc))

    device_type = np.full(nloc, "OOPSBUG", dtype="a8")
    device_type[:] = fp["DEVICE_TYPE"][keep_rows]

    locations = np.copy(fp["LOCATION"][keep_rows])

    # Map location to row in the table

    loc_to_fp = dict()
    for rw, loc in enumerate(fp["LOCATION"]):
        loc_to_fp[loc] = rw

    # FIXME:  Here we assume that the 32bit STATE column has the same bit
    # definitions as what is used by fiberassign (defined in hardware.h):
    # If this is not true, then re-map those values here inside the "state"
    # table loaded above.

    # Map location to row of the state table

    loc_to_state = dict()
    for rw, loc in enumerate(state["LOCATION"]):
        loc_to_state[loc] = rw

    # Convert the exclusion polygons into shapes.

    excl = dict()

    for nm, shp in exclude.items():
        excl[nm] = dict()
        for arm in ["theta", "phi"]:
            cr = list()
            for crc in shp[arm]["circles"]:
                cr.append(Circle(crc[0], crc[1]))
            sg = list()
            for sgm in shp[arm]["segments"]:
                sg.append(Segments(sgm))
            fshp = Shape((0.0, 0.0), cr, sg)
            excl[nm][arm] = fshp

    # For each positioner, select the exclusion polynomials.

    positioners = dict()

    for loc in locations:
        exclname = state["EXCLUSION"][loc_to_state[loc]]
        positioners[loc] = dict()
        positioners[loc]["theta"] = Shape(excl[exclname]["theta"])
        positioners[loc]["phi"] = Shape(excl[exclname]["phi"])

    hw = Hardware(locations,
                  fp["PETAL"][keep_rows],
                  fp["DEVICE"][keep_rows],
                  fp["SLITBLOCK"][keep_rows],
                  fp["BLOCKFIBER"][keep_rows],
                  fp["FIBER"][keep_rows],
                  device_type,
                  fp["OFFSET_X"][keep_rows],
                  fp["OFFSET_Y"][keep_rows],
                  np.array([state["STATE"][loc_to_state[x]]
                           for x in locations]),
                  np.array([fp["OFFSET_T"][loc_to_fp[x]] for x in locations]),
                  np.array([fp["MIN_T"][loc_to_fp[x]] for x in locations]),
                  np.array([fp["MAX_T"][loc_to_fp[x]] for x in locations]),
                  np.array([fp["LENGTH_R1"][loc_to_fp[x]] for x in locations]),
                  np.array([fp["OFFSET_P"][loc_to_fp[x]] for x in locations]),
                  np.array([fp["MIN_P"][loc_to_fp[x]] for x in locations]),
                  np.array([fp["MAX_P"][loc_to_fp[x]] for x in locations]),
                  np.array([fp["LENGTH_R2"][loc_to_fp[x]] for x in locations]),
                  [positioners[x]["theta"] for x in locations],
                  [positioners[x]["phi"] for x in locations])
    return hw
