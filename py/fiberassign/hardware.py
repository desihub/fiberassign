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
import astropy.table

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

def get_default_exclusion_margins():
    '''
    We add an additional margin around the exclusion zones in the hardware descriptions.
    This function returns the defaults used by the fiberassign scripts.

    Returns:
    margins (dict): with keys "pos", "gfa" and "petal", to floating-point margins in millimeters.
    '''
    return dict(pos=0.05,
                petal=0.4,
                gfa=0.4)

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

# Global cache of exclusions.  (exclusions are geometric regions describing the keep-out
# areas around positioners, GFAs, and petal edges)
_cache_shape_to_excl = dict()

def get_exclusions(exclude, exclname, margins, local_cache=None):
    '''
    Fetches an exclusion region from a data file, applying additional margins as required.

    A "local_cache", valid for a single "exclude" file, can be passed in.

    Args:
    * exclude: the exclusions file, eg from desimodel.io.load_focalplane.
    * exclname: the string name to look up
    * margins: dict of additional margins to apply
    * local_cache: dict
    '''
    global _cache_shape_to_excl

    if local_cache is not None:
        # check the local cache, where names are unique
        e = local_cache.get(exclname)
        if e is not None:
            return e

    shp = exclude[exclname]

    # NOTE that different exclusions files can assign
    # different values to the same key!  ie, code 0000fbb5303ebbc2 appears in
    #   desi-exclusion_2021-03-17T23:20:01.json and
    #   desi-exclusion_2021-06-25T22:38:59+00:00.json
    # with totally different values.
    # We therefore cannot use the "exclname" as the (global) cache key.
    #
    # I tried caching the "expand_closed_curve" computation using a hash of the
    # shape itself as the cache key, but that turned out not to help -- computing
    # the key was too expensive?
    #
    # For further speedups, I suspect that moving expand_closed_curve to C++ would be
    # very effective.

    rtn = dict()
    for obj in shp.keys():
        cr = list()
        for crc in shp[obj]["circles"]:
            cr.append(Circle(crc[0], crc[1]))
        sg = list()
        for sgm in shp[obj]["segments"]:
            if obj in margins:
                # Create a hashable version of the line-segment list
                key = (obj, tuple(tuple(xy) for xy in sgm))
                cached = _cache_shape_to_excl.get(key)
                if cached is not None:
                    sgm = cached
                else:
                    sx = [x for x,y in sgm]
                    sy = [y for x,y in sgm]
                    ex,ey = expand_closed_curve(sx, sy, margins[obj])
                    sgm = list(zip(ex, ey))
                    _cache_shape_to_excl[key] = sgm
            sg.append(Segments(sgm))
        fshp = Shape((0.0, 0.0), cr, sg)
        rtn[obj] = fshp
    if local_cache is not None:
        local_cache[exclname] = rtn
    return rtn

def load_hardware(focalplane=None, rundate=None, add_margins={},
                  get_time_range=False):
    """Create a hardware class representing properties of the telescope.

    Args:
        focalplane (tuple):  Override the focalplane model.  If not None, this
            should be a tuple of the same data types returned by
            desimodel.io.load_focalplane()
        rundate (str):  ISO 8601 format time stamp as a string in the
            format YYYY-MM-DDTHH:MM:SS+-zz:zz.  If None, uses current time.
        add_margins (dict): additional margins to add around positioners, GFAs,
            and petals.  Dict with keys "pos", "gfa", "petal" and values of
            millimeters.
        get_time_range (bool): if True, return (hw, time_lo, time_hi), where
            time_lo and time_hi are datetime objects corresponding to the first and
            last dates when the hardware was in this state.

    Returns:
        if get_time_range is True: (hardware, time_lo, time_hi)
        else: hardware

        (Hardware):  The hardware object.
    """
    args = load_hardware_args(focalplane=focalplane, rundate=rundate, add_margins=add_margins)
    args, time_lo, time_hi = args
    hw = Hardware(*args)
    if get_time_range:
        return hw, time_lo, time_hi
    return hw

def load_hardware_args(focalplane=None, rundate=None, add_margins={}):
    """Reads the arguments needed to create a Hardware object representing the
    properties of the instrument at a given date.  Also returns the range of
    dates when this hardware configuration is valid.

    Args:
        focalplane (tuple):  Override the focalplane model.  If not None, this
            should be a tuple of the same data types returned by
            desimodel.io.load_focalplane()
        rundate (str):  ISO 8601 format time stamp as a string in the
            format YYYY-MM-DDTHH:MM:SS+-zz:zz.  If None, uses current time.
        add_margins (dict): additional margins to add around positioners, GFAs,
            and petals.  Dict with keys "pos", "gfa", "petal" and values of
            millimeters.

    Returns:
        args (tuple): The arguments to be passed to the Hardware constructor
        time_lo (datetime): The earliest date when this hardware configuration is valid
        time_hi (datetime): The latest date when this hardware configuration is valid
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
    time_lo = time_hi = None
    if focalplane is None:
        fp, exclude, state, create_time, time_lo, time_hi = dmio.load_focalplane(runtime, get_time_range=True)
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


    # For each positioner, select the exclusion polynomials.
    positioners = dict()
    excl_cache = dict()
    for loc in locations:
        exclname = state["EXCLUSION"][loc_to_state[loc]]
        positioners[loc] = dict()
        posexcl = get_exclusions(exclude, exclname, margins, local_cache=excl_cache)
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
    del excl_cache

    if "MIN_P" in state.colnames:
        # This is a new-format focalplane model (after desimodel PR #143)
        Istate = np.array([loc_to_state[x] for x in locations])
        Ifp    = np.array([loc_to_fp   [x] for x in locations])
        args = (
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
            state['STATE'][Istate],
            fp['OFFSET_T'][Ifp],
            state['MIN_T'][Istate],
            state['MAX_T'][Istate],
            state['POS_T'][Istate],
            fp['LENGTH_R1'][Ifp],
            fp['OFFSET_P'][Ifp],
            state['MIN_P'][Istate],
            state['MAX_P'][Istate],
            state['POS_P'][Istate],
            fp['LENGTH_R2'][Ifp],
            fine_radius,
            fine_theta,
            fine_arc,
            [positioners[x]["theta"] for x in locations],
            [positioners[x]["phi"]   for x in locations],
            [positioners[x]["gfa"]   for x in locations],
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

        args = (
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

    aa = []
    for a in args:
        if type(a) == astropy.table.Column:
            a = a.data
        aa.append(a)

    return aa, time_lo, time_hi

def radec2xy(hw, tile_ra, tile_dec, tile_obstime, tile_obstheta, tile_obsha,
             ra, dec, use_cs5, threads=0, use_hardcoded_polmis_rotmat=True):
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
    from inspect import getfullargspec

    # Check for the precomputed polar misalignment option in desispec.
    precomputed_matrix = 'use_hardcoded_polmis_rotmat' in getfullargspec(fiberassign_radec2xy_cs5).args

    # Note that MJD is only used for precession, so no need for
    # high precision.
    t = Time(tile_obstime, format='isot')
    mjd = t.mjd

    # Don't pass adc[12]: Let desimeter use its pm-alike routines
    if precomputed_matrix:
        if use_cs5:
            x, y = fiberassign_radec2xy_cs5(ra, dec, tile_ra, tile_dec, mjd,
                                            tile_obsha, tile_obstheta, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)
        else:
            x, y = fiberassign_radec2xy_flat(ra, dec, tile_ra, tile_dec, mjd,
                                             tile_obsha, tile_obstheta, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)
    else:
        log = Logger.get()
        log.warning('fiberassign_radec2xy_(cs5|flat) without use_hardcoded_polmis_rotmat option is deprecated')

        if use_cs5:
            x, y = fiberassign_radec2xy_cs5(ra, dec, tile_ra, tile_dec, mjd,
                                            tile_obsha, tile_obstheta)
        else:
            x, y = fiberassign_radec2xy_flat(ra, dec, tile_ra, tile_dec, mjd,
                                             tile_obsha, tile_obstheta)
    return x,y

def xy2radec(hw, tile_ra, tile_dec, tile_obstime, tile_obstheta, tile_obsha,
             x, y, use_cs5, threads=0, use_hardcoded_polmis_rotmat=True):
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
    from inspect import getfullargspec
    from astropy.time import Time

    # Check if we can use the precomputed matrix option.
    precomputed_matrix = 'use_hardcoded_polmis_rotmat' in getfullargspec(fiberassign_cs5_xy2radec).args

    t = Time(tile_obstime, format='isot')
    mjd = t.mjd
    if precomputed_matrix:
        if use_cs5:
            ra,dec = fiberassign_cs5_xy2radec(x, y, tile_ra, tile_dec, mjd,
                                              tile_obsha, tile_obstheta, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)
        else:
            ra,dec = fiberassign_flat_xy2radec(x, y, tile_ra, tile_dec, mjd,
                                               tile_obsha, tile_obstheta, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)
    else:
        log = Logger.get()
        log.warning('fiberassign_(cs5|flat)_xy2radec without use_hardcoded_polmis_rotmat option is deprecated')

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

def cs52xy(x, y):
    '''
    Converts from CS5 coordinates (mm) into curved focal-plane X,Y coordinates in mm.

    Args:
    cs5x (numpy array): CS5 X coord (mm)
    cs5y (numpy array): CS5 Y coord (mm)

    Returns:
    x (numpy array): X coord (mm)
    y (numpy array): Y coord (mm)
    '''
    # There's a change in terminology between the focal-plane team and
    # the outside world here...
    from desimeter.transform.pos2ptl import ptl2flat
    return ptl2flat(x, y)
