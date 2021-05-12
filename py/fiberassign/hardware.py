# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.hardware
=======================

Functions for loading information about the telescope hardware.

"""
from __future__ import absolute_import, division, print_function

from datetime import datetime, timezone

import numpy as np

from scipy.interpolate import interp1d

import desimodel.io as dmio

from .utils import Logger

from ._internal import (
    Hardware,
    FIBER_STATE_OK,
    FIBER_STATE_UNASSIGNED,
    FIBER_STATE_STUCK,
    FIBER_STATE_BROKEN,
    FIBER_STATE_RESTRICT,
    Circle,
    Segments,
    Shape,
)


def load_hardware(focalplane=None, rundate=None):
    """Create a hardware class representing properties of the telescope.

    Args:
        focalplane (tuple):  Override the focalplane model.  If not None, this
            should be a tuple of the same data types returned by
            desimodel.io.load_focalplane()
        rundate (str):  ISO 8601 format time stamp as a string in the
            format YYYY-MM-DDTHH:MM:SS+-zz:zz.  If None, uses current time.

    Returns:
        (Hardware):  The hardware object.

    """
    log = Logger.get()

    # The timestamp for this run.
    runtime = None
    if rundate is None:
        runtime = datetime.now(tz=timezone.utc)
    else:
        try:
            runtime = datetime.strptime(rundate, "%Y-%m-%dT%H:%M:%S%z")
        except ValueError:
            runtime = datetime.strptime(rundate, "%Y-%m-%dT%H:%M:%S")
            msg = "Requested run date '{}' is not timezone-aware.  Assuming UTC.".format(runtime)
            log.warning(msg)
            runtime = runtime.replace(tzinfo=timezone.utc)
    runtimestr = None
    try:
        runtimestr = runtime.isoformat(timespec="seconds")
    except TypeError:
        runtimestr = runtime.isoformat()

    # Get the focalplane information
    fp = None
    exclude = None
    state = None
    create_time = "UNKNOWN"
    if focalplane is None:
        fp, exclude, state, create_time = dmio.load_focalplane(runtime)
    else:
        fp, exclude, state = focalplane

    # Get the plate scale
    platescale = dmio.load_platescale()

    # We are going to do a quadratic interpolation to the platescale on a fine grid,
    # and then use that for *linear* interpolation inside the compiled code.  The
    # default platescale data is on a one mm grid spacing.  We also do the same
    # interpolation of the arclength S(R).

    fine_radius = np.linspace(
        platescale["radius"][0], platescale["radius"][-1], num=10000, dtype=np.float64
    )
    fn = interp1d(platescale["radius"], platescale["theta"], kind="quadratic")
    fine_theta = fn(fine_radius).astype(np.float64)
    fn = interp1d(platescale["radius"], platescale["arclength"], kind="quadratic")
    fine_arc = fn(fine_radius).astype(np.float64)

    # We are only going to keep rows for LOCATIONs that are assigned to a
    # science or sky monitor positioner.

    log.info("Loaded focalplane for time stamp {}".format(runtime))

    pos_rows = np.where(fp["DEVICE_TYPE"].astype(str) == "POS")[0]
    etc_rows = np.where(fp["DEVICE_TYPE"].astype(str) == "ETC")[0]
    keep_rows = np.unique(np.concatenate((pos_rows, etc_rows)))

    nloc = len(keep_rows)
    log.debug("  focalplane table keeping {} rows for POS and ETC devices".format(nloc))

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
        for obj in shp.keys():
            cr = list()
            for crc in shp[obj]["circles"]:
                cr.append(Circle(crc[0], crc[1]))
            sg = list()
            for sgm in shp[obj]["segments"]:
                sg.append(Segments(sgm))
            fshp = Shape((0.0, 0.0), cr, sg)
            excl[nm][obj] = fshp

    # For each positioner, select the exclusion polynomials.

    positioners = dict()

    for loc in locations:
        exclname = state["EXCLUSION"][loc_to_state[loc]]
        positioners[loc] = dict()
        positioners[loc]["theta"] = Shape(excl[exclname]["theta"])
        positioners[loc]["phi"] = Shape(excl[exclname]["phi"])
        if "gfa" in excl[exclname]:
            positioners[loc]["gfa"] = Shape(excl[exclname]["gfa"])
        else:
            positioners[loc]["gfa"] = Shape()
        if "petal" in excl[exclname]:
            positioners[loc]["petal"] = Shape(excl[exclname]["petal"])
        else:
            positioners[loc]["petal"] = Shape()

    hw = None
    if "MIN_P" in state.colnames:
        # This is a new-format focalplane model (after desimodel PR #143)
        hw = Hardware(
            runtimestr,
            locations,
            fp["PETAL"][keep_rows],
            fp["DEVICE"][keep_rows],
            fp["SLITBLOCK"][keep_rows],
            fp["BLOCKFIBER"][keep_rows],
            fp["FIBER"][keep_rows],
            device_type,
            fp["OFFSET_X"][keep_rows],
            fp["OFFSET_Y"][keep_rows],
            np.array([state["STATE"][loc_to_state[x]] for x in locations]),
            np.array([fp["OFFSET_T"][loc_to_fp[x]] for x in locations]),
            np.array([state["MIN_T"][loc_to_state[x]] for x in locations]),
            np.array([state["MAX_T"][loc_to_state[x]] for x in locations]),
            np.array([state["POS_T"][loc_to_state[x]] for x in locations]),
            np.array([fp["LENGTH_R1"][loc_to_fp[x]] for x in locations]),
            np.array([fp["OFFSET_P"][loc_to_fp[x]] for x in locations]),
            np.array([state["MIN_P"][loc_to_state[x]] for x in locations]),
            np.array([state["MAX_P"][loc_to_state[x]] for x in locations]),
            np.array([state["POS_P"][loc_to_state[x]] for x in locations]),
            np.array([fp["LENGTH_R2"][loc_to_fp[x]] for x in locations]),
            fine_radius,
            fine_theta,
            fine_arc,
            [positioners[x]["theta"] for x in locations],
            [positioners[x]["phi"] for x in locations],
            [positioners[x]["gfa"] for x in locations],
            [positioners[x]["petal"] for x in locations],
        )
    else:
        # This is an old-format focalplane model (prior to desimodel PR #143).  For
        # stuck positioners, we want to specify a default POS_T / POS_P to use.
        # These old models did not include any information about that, so we use the
        # minimum Theta value and either the maximum Phi value or PI, whichever is
        # smaller
        fake_pos_p = np.zeros(len(locations), dtype=np.float64)
        fake_pos_t = np.zeros(len(locations), dtype=np.float64)
        for ilid, lid in enumerate(locations):
            pt = fp["MIN_T"][loc_to_fp[lid]] + fp["OFFSET_T"][loc_to_fp[lid]]
            pp = fp["MAX_P"][loc_to_fp[lid]] + fp["OFFSET_P"][loc_to_fp[lid]]
            if pp > 180.0:
                pp = 180.0
            fake_pos_p[ilid] = pp
            fake_pos_t[ilid] = pt
        hw = Hardware(
            runtimestr,
            locations,
            fp["PETAL"][keep_rows],
            fp["DEVICE"][keep_rows],
            fp["SLITBLOCK"][keep_rows],
            fp["BLOCKFIBER"][keep_rows],
            fp["FIBER"][keep_rows],
            device_type,
            fp["OFFSET_X"][keep_rows],
            fp["OFFSET_Y"][keep_rows],
            np.array([state["STATE"][loc_to_state[x]] for x in locations]),
            np.array([fp["OFFSET_T"][loc_to_fp[x]] for x in locations]),
            np.array([fp["MIN_T"][loc_to_fp[x]] for x in locations]),
            np.array([fp["MAX_T"][loc_to_fp[x]] for x in locations]),
            fake_pos_t,
            np.array([fp["LENGTH_R1"][loc_to_fp[x]] for x in locations]),
            np.array([fp["OFFSET_P"][loc_to_fp[x]] for x in locations]),
            np.array([fp["MIN_P"][loc_to_fp[x]] for x in locations]),
            np.array([fp["MAX_P"][loc_to_fp[x]] for x in locations]),
            fake_pos_p,
            np.array([fp["LENGTH_R2"][loc_to_fp[x]] for x in locations]),
            fine_radius,
            fine_theta,
            fine_arc,
            [positioners[x]["theta"] for x in locations],
            [positioners[x]["phi"] for x in locations],
            [positioners[x]["gfa"] for x in locations],
            [positioners[x]["petal"] for x in locations],
        )
    return hw

def radec2xy(hw, tile_ra, tile_dec, tile_obstime, tile_obstheta, tile_obsha,
             ra, dec, use_cs5, threads=0):
    #xy = hw.radec2xy_multi(
    #    tile_ra, tile_dec, tile_obstheta, ra, dec, use_cs5, threads=0
    #)
    #x = np.array([x for x,y in xy])
    #y = np.array([y for x,y in xy])
    from astropy.time import Time
    from desimeter.fiberassign import fiberassign_radec2xy_cs5, fiberassign_radec2xy_flat
    # Note that MJD is only used for precession, so no need for
    # high precision.
    t = Time(tile_obstime, format='isot')
    mjd = t.mjd

    # Don't pass adc[12]: Let desimeter use its pm-alike routines
    if use_cs5:
        x, y = fiberassign_radec2xy_cs5(ra, dec, tile_ra, tile_dec, mjd,
                                        tile_obsha, tile_obstheta)
    else:
        x, y = fiberassign_radec2xy_flat(ra, dec, tile_ra, tile_dec, mjd,
                                         tile_obsha, tile_obstheta)
    return x,y

def xy2radec(hw, tile_ra, tile_dec, tile_obstime, tile_obstheta, tile_obsha,
             x, y, use_cs5, threads=0):
    # radec = hw.xy2radec_multi(
    #     tile_ra, tile_dec, tile_obstheta, x, y, use_cs5, threads
    #     )
    # ra  = np.array([r for r,d in radec])
    # dec = np.array([d for r,d in radec])
    from desimeter.fiberassign import fiberassign_cs5_xy2radec, fiberassign_flat_xy2radec
    from astropy.time import Time
    t = Time(tile_obstime, format='isot')
    mjd = t.mjd
    if use_cs5:
        ra,dec = fiberassign_cs5_xy2radec(x, y, tile_ra, tile_dec, mjd,
                                          tile_obsha, tile_obstheta)
    else:
        ra,dec = fiberassign_flat_xy2radec(x, y, tile_ra, tile_dec, mjd,
                                           tile_obsha, tile_obstheta)
    return ra,dec

def xy2cs5(x, y):
    # There's a change in terminology between the focal-plane team and
    # the outside world here...
    from desimeter.transform.pos2ptl import ptl2flat
    return flat2ptl(x, y)
