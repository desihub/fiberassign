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

def expand_closed_curve(xx, yy, margin):
    '''
    For the --margin-{pos,petal,gfa} options, we can add a buffer zone
    around the positioner keep-out polygons.  This function implements
    the geometry to achieve this.

    Given a RIGHT-HANDED closed polygon xx,yy and margin, returns new
    x,y coordinates expanded by a margin of `margin`.

    (By right-handed, I mean that the points are listed
    counter-clockwise, and if you walk the boundary, the inside of the
    shape is to the left; the expanded points will be to the right.)

    If the order of the polygon is reversed, the "expanded" points
    will actually be on the *inside* of the polygon.  Setting the
    margin negative counteracts this.

    Note that we strictly require a closed curve.  Collinear polygon
    segments will cause problems!

    '''
    ex, ey = [],[]

    # These are closed curves (last point = first point)
    # (this isn't strictly required by the fundamental algorithm, but is assumed
    # in the way we select previous and next points in the loop below.)
    if (xx[0] != xx[-1]) or (yy[0] != yy[-1]):
        log = Logger.get()
        log.warning('Expected exclusion polygons to be closed curves; got x, y = %s, %s' % (str(xx), str(yy)))
        return xx, yy
    assert(xx[0] == xx[-1])
    assert(yy[0] == yy[-1])

    N = len(xx)

    for j in range(N):
        # We go through the points, and for each point we consider the previous
        # and next point.  The "expanded" point will be defined according to the
        # two edges (vector) coming from the point.
        i = j - 1
        if i == -1:
            # wrap around, skipping repeated point
            i = N-2
        k = j + 1
        if k == N:
            k = 1

        x1 = xx[i]
        y1 = yy[i]
        x2 = xx[j]
        y2 = yy[j]
        x3 = xx[k]
        y3 = yy[k]

        # Vectors to and from the central point.
        vx1 = x2 - x1
        vy1 = y2 - y1
        vx2 = x3 - x2
        vy2 = y3 - y2
        # We can't handle repeated points!
        assert(not(vx2 == 0. and vy2 == 0.))

        # Get the angle between the vectors -- our expanded point is going to
        # be halfway between these two vectors.
        cross1 = vx1 * vy2 - vx2 * vy1
        vv1 = np.hypot(vx1,vy1)
        vv2 = np.hypot(vx2,vy2)
        theta = np.arcsin(np.clip(cross1 / (vv1 * vv2), -1., 1.))
        # Detect sharp (>90 degree) turns
        dot = vx1*vx2 + vy1*vy2
        # the angle of the expanded point is relative to vector 1
        a = np.arctan2(vy1, vx1)
        if dot < 0:
            # sharp turn -- the theta=arcsin aliases the angle, which is outside
            # the range of [-pi/2, +pi/2]; adjust theta to the aliased angle.
            if theta > 0:
                theta =  np.pi - theta
            else:
                theta = -np.pi - theta
        da = np.pi/2. + theta/2.
        # This places the point further from the original keeps the vectors parallel to their originals
        stretch = 1./np.cos(theta/2.)
        dx = -margin * stretch * np.cos(a + da)
        dy = -margin * stretch * np.sin(a + da)

        ex.append(x2 + dx)
        ey.append(y2 + dy)

    return ex, ey


def load_hardware(focalplane=None, rundate=None,
                  add_margins={}):
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

    # Slightly reformat the 'add_margins' dict: 'pos' gets copied to
    # 'theta' and 'phi'.
    margins = {}
    if 'pos' in add_margins:
        margins['theta'] = add_margins['pos']
        margins['phi']   = add_margins['pos']
    if 'gfa' in add_margins:
        margins['gfa'] = add_margins['gfa']
    # PLUS -- and this is a HACK --
    # because the 'petal' exclusion polygons are listed in the db dump files
    # in "left-handed" order, unlike the other polygons, we must NEGATE the
    # margin!
    if 'petal' in add_margins:
        margins['petal'] = -add_margins['petal']

    # Convert the exclusion polygons into shapes (as required)
    excl = dict()
    # cache expanded polygons (because many of the polygons are actually duplicates)
    expanded = {}
    nhit = 0
    nexp = 0
    def get_exclusions(exclname):

        nonlocal nhit
        nonlocal nexp

        e = excl.get(exclname)
        if e is not None:
            return e
        shp = exclude[exclname]
        excl[exclname] = dict()
        for obj in shp.keys():
            cr = list()
            for crc in shp[obj]["circles"]:
                cr.append(Circle(crc[0], crc[1]))
            sg = list()
            for sgm in shp[obj]["segments"]:
                if obj in margins:
                    key = (obj, sgm)
                    if key in expanded:
                        nhit += 1
                        sgm = expanded[key]
                    else:
                        sx = [x for x,y in sgm]
                        sy = [y for x,y in sgm]
                        ex,ey = expand_closed_curve(sx, sy, margins[obj])
                        sgm = list(zip(ex, ey))
                        expanded[key] = sgm
                        nexp += 1
                sg.append(Segments(sgm))
            fshp = Shape((0.0, 0.0), cr, sg)
            excl[exclname][obj] = fshp
        return excl[exclname]

    # For each positioner, select the exclusion polynomials.
    positioners = dict()

    for loc in locations:
        exclname = state["EXCLUSION"][loc_to_state[loc]]
        positioners[loc] = dict()
        posexcl = get_exclusions(exclname)
        positioners[loc]["theta"] = Shape(posexcl["theta"])
        positioners[loc]["phi"] = Shape(posexcl["phi"])
        if "gfa" in posexcl:
            positioners[loc]["gfa"] = Shape(posexcl["gfa"])
        else:
            positioners[loc]["gfa"] = Shape()
        if "petal" in posexcl:
            positioners[loc]["petal"] = Shape(posexcl["petal"])
        else:
            positioners[loc]["petal"] = Shape()

    print('Polygons expanded:', nexp)
    print('N cached:', nhit)

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
            add_margins
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
            add_margins
        )
    return hw

def radec2xy(hw, tile_ra, tile_dec, tile_obstime, tile_obstheta, tile_obsha,
             ra, dec, use_cs5, threads=0):
    '''
    For the tile pointed at (tilera, tiledec), project the (ra, dec)
    value into X/Y mm.

    Args:
      hw: Hardware object
      tile_ra (float): Tile RA
      tile_dec (float): Tile Dec
      tile_obstime (string): Tile observation time, YYYY-MM-DDTHH:MM:SS / Astropy "isot" format.
      tile_obstheta (float): Tile "fieldrot" rotation angle.
      tile_obsha (float): Tile designed Hour Angle, in degrees.
      ra (numpy array): RA to project, in degrees
      dec (numpy array): Dec to project, in degrees
      use_CS5 (bool):  If True, return CS5 coordinates, else curved.
      threads=0 (int): currently unused; for backward compatibility.

    Returns:
      x, y: numpy arrays: the (X, Y) projected locations.
    '''
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
    '''
    For the tile pointed at (tilera, tiledec), compute the RA,Dec
    pointing of the specified X/Y location in millimeters.

    Args:
      hw: Hardware object
      tile_obstime (string): Tile observation time, YYYY-MM-DDTHH:MM:SS / Astropy "isot" format.
      tile_obstheta (float): Tile "fieldrot" rotation angle.
      tile_obsha (float): Tile designed Hour Angle, in degrees.
      x (numpy array): X position in mm.
      y (numpy array): Y position in mm.
      use_CS5 (bool):  If True, assume X,Y are CS5 coordinates, else curved.
      threads=0 (int): currently unused; for backward compatibility.

    Returns:
      ra, dec (numpy arrays): the (RA, Dec) values of the focalplane locations,
                              in degrees.
    '''
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
    '''
    Converts from curved focal-plane X,Y coordinates in mm into CS5
    coordinates in mm.

    Args:
    x (numpy array): X coord (mm)
    y (numpy array): Y coord (mm)

    Returns:
    cs5x (numpy array): CS5 X coord (mm)
    cs5y (numpy array): CS5 Y coord (mm)
    '''
    # There's a change in terminology between the focal-plane team and
    # the outside world here...
    from desimeter.transform.pos2ptl import flat2ptl
    return flat2ptl(x, y)
