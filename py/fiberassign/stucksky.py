import sys
import os
from datetime import datetime
import numpy as np
import fitsio

from desitarget.skybricks import Skybricks
from desitarget.skyhealpixs import Skyhealpixs
from astropy.time import Time

def stuck_on_sky(hw, tiles, lookup_sky_source, rundate=None):
    '''
    Will STUCK positioners land on good SKY locations for the given set of tiles?

    Args:
        hw: hardware properties loaded by fiberassign.hardware.load_hardware()
        tiles: tiles loaded by fiberassign.tiles.load_tiles()
        lookup_sky_source: "ls" or "gaia"

    Returns a nested dict:
        stuck_sky[tileid][loc] = bool_good_sky
    '''
    from fiberassign.utils import Logger, get_date_cutoff
    from fiberassign.hardware import (FIBER_STATE_STUCK, FIBER_STATE_BROKEN,
                                      xy2radec)

    log = Logger.get()

    # AR allowed lookup_sky_source
    if lookup_sky_source not in ["ls", "gaia"]:
        log.error("lookup_sky_source = {} not in ['ls', 'gaia']; exiting".format(lookup_sky_source))
        sys.exit(1)

    # Only load the Skybricks object the first time we need it (there might not be any stuck positioners, eg in tests)
    # AR keep similar approach for Skyhealpixs...
    if lookup_sky_source == "ls":
        skybricks = None
        log.info("lookup_sky_source == 'ls', will look for $SKYBRICKS_DIR")
    if lookup_sky_source == "gaia":
        skyhealpixs = None
        log.info("lookup_sky_source == 'gaia', will look for $SKYHEALPIXS_DIR")

    stuck_sky = dict()
    if rundate is None:
        etc_fibers_are_stuck = False
    else:
        etc_cutoff = get_date_cutoff('rundate', 'etc_stuck')
        etc_cutoff = Time(datetime.strptime(etc_cutoff, "%Y-%m-%dT%H:%M:%S%z")).mjd
        rundate = Time(datetime.strptime(rundate,  "%Y-%m-%dT%H:%M:%S%z")).mjd
        etc_fibers_are_stuck = rundate > etc_cutoff

    for tile_id, tile_ra, tile_dec, tile_obstime, tile_theta, tile_obsha in zip(
            tiles.id, tiles.ra, tiles.dec, tiles.obstime, tiles.obstheta, tiles.obshourang):
        stuck_sky[tile_id] = dict()

        # Stuck locations and their angles
        # (grab the hw dictionaries once -- these are python wrappers over C++ so not simple accessors)
        state = hw.state
        devtype = hw.loc_device_type
        stuck_loc = []
        for loc in hw.locations:
            if (devtype[loc] == 'ETC') and etc_fibers_are_stuck:
                stuck_loc.append(loc)
            elif devtype[loc] != 'POS':
                continue
            if (state[loc] & (FIBER_STATE_STUCK | FIBER_STATE_BROKEN)) == FIBER_STATE_STUCK:
                stuck_loc.append(loc)
        if len(stuck_loc) == 0:
            log.debug('Tile %i: %i positioners are stuck/broken' % (tile_id, len(stuck_loc)))
            continue
        theta_pos = hw.loc_theta_pos
        theta_off = hw.loc_theta_offset
        phi_pos = hw.loc_phi_pos
        phi_off = hw.loc_phi_offset
        stuck_theta = [theta_pos[loc] + theta_off[loc] for loc in stuck_loc]
        stuck_phi   = [phi_pos  [loc] + phi_off  [loc] for loc in stuck_loc]

        # Convert positioner angle orientations to curved focal surface X / Y (not CS5)
        # Note:  we could add some methods to the python bindings to vectorize this or make it less clunky...
        stuck_x = np.zeros(len(stuck_loc))
        stuck_y = np.zeros(len(stuck_loc))
        curved_mm = hw.loc_pos_curved_mm
        theta_arm = hw.loc_theta_arm
        phi_arm   = hw.loc_phi_arm
        theta_min = hw.loc_theta_min
        theta_max = hw.loc_theta_max
        phi_min   = hw.loc_phi_min
        phi_max   = hw.loc_phi_max
        for iloc, (loc, theta, phi) in enumerate(zip(stuck_loc, stuck_theta, stuck_phi)):
            loc_x, loc_y = hw.thetaphi_to_xy(
                curved_mm[loc],
                theta,
                phi,
                theta_arm[loc],
                phi_arm[loc],
                theta_off[loc],
                phi_off[loc],
                theta_min[loc],
                phi_min[loc],
                theta_max[loc],
                phi_max[loc],
                True
            )
            stuck_x[iloc] = loc_x
            stuck_y[iloc] = loc_y

        loc_ra,loc_dec = xy2radec(
            hw, tile_ra, tile_dec, tile_obstime, tile_theta, tile_obsha,
            stuck_x, stuck_y, False, 0
        )

        # AR using the legacysurveys good locations
        if lookup_sky_source == "ls":
            if skybricks is None:
                try:
                    skybricks = Skybricks()
                except:
                    log.warning('Environment variable SKYBRICKS_DIR is not set; not looking '
                                'up whether stuck positioners land on good sky')
                    return
            good_sky = skybricks.lookup_tile(tile_ra, tile_dec, hw.focalplane_radius_deg,
                                             loc_ra, loc_dec)

        # AR using the gaia good locations
        if lookup_sky_source == "gaia":
            if skyhealpixs is None:
                try:
                    skyhealpixs = Skyhealpixs()
                except:
                    log.warning('Environment variable SKYHEALPIXS_DIR is not set; not looking '
                                'up whether stuck positioners land on good sky')
                    return
            good_sky = skyhealpixs.lookup_position(loc_ra, loc_dec)

        stuck_isetc = np.array([devtype[loc] == 'ETC' for loc in stuck_loc])
        log.info('%i of %i stuck positioners land on good sky locations' %
                 (np.sum(good_sky & ~stuck_isetc), np.sum(~stuck_isetc)))
        if etc_fibers_are_stuck:
            log.info('%i of %i ETC positioners land on good sky locations' %
                     (np.sum(good_sky & stuck_isetc), np.sum(stuck_isetc)))
        for loc,good in zip(stuck_loc, good_sky):
            stuck_sky[tile_id][loc] = good

    return stuck_sky


def stuck_on_sky_from_fafns(fafns):
    '''
    Retrieve the information if STUCK positioners land on good SKY locations from a list of fiberassign-TILEID.fits.gz files.

    Args:
        fafns: comma-separated list of full paths to fiberassign-TILEID.fits.gz files.

    Returns a nested dict:
        stuck_sky[tileid][loc] = bool_good_sky
    '''

    from fiberassign.utils import Logger, get_date_cutoff
    from fiberassign.hardware import FIBER_STATE_STUCK, FIBER_STATE_BROKEN
    from fiberassign.targets import TARGET_TYPE_SKY

    stuck_sky = dict()

    for fafn in fafns.split(","):

        hdr = fitsio.read_header(fafn, 0)
        tile_id = hdr["TILEID"]
        rundate = hdr["RUNDATE"]

        if rundate is None:
            etc_fibers_are_stuck = False
        else:
            etc_cutoff = get_date_cutoff('rundate', 'etc_stuck')
            etc_cutoff = Time(datetime.strptime(etc_cutoff, "%Y-%m-%dT%H:%M:%S%z")).mjd
            rundate = Time(datetime.strptime(rundate,  "%Y-%m-%dT%H:%M:%S%z")).mjd
            etc_fibers_are_stuck = rundate > etc_cutoff

        stuck_sky[tile_id] = dict()
        # FIBERASSIGN
        d = fitsio.read(fafn, "FIBERASSIGN", columns=["LOCATION", "FIBER", "FIBERSTATUS", "OBJTYPE"])
        sel = (d["FIBERSTATUS"] & (FIBER_STATE_STUCK | FIBER_STATE_BROKEN)) == FIBER_STATE_STUCK
        for loc, good in zip(d["LOCATION"][sel], d["OBJTYPE"][sel] == "SKY"):
            stuck_sky[tile_id][loc] = good
        # ETC
        if etc_fibers_are_stuck:
            d = fitsio.read(fafn, "SKY_MONITOR", columns=["LOCATION", "FA_TYPE"])
            for loc, good in zip(d["LOCATION"], (d["FA_TYPE"] & TARGET_TYPE_SKY) > 0):
                stuck_sky[tile_id][loc] = good
        # re-order by increasing locations, to reproduce stuck_on_sky()
        stuck_sky[tile_id] = dict(sorted(stuck_sky[tile_id].items()))

    return stuck_sky
