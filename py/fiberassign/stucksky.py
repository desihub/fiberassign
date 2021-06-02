import os
import numpy as np

from desitarget.skybricks import Skybricks

def stuck_on_sky(hw, tiles):
    '''
    Will STUCK positioners land on good SKY locations for the given set of tiles?

    Returns a nested dict:
        stuck_sky[tileid][loc] = bool_good_sky
    '''
    from fiberassign.utils import Logger
    from fiberassign.hardware import (FIBER_STATE_STUCK, FIBER_STATE_BROKEN,
                                      xy2radec)

    log = Logger.get()

    # Only load the Skybricks object the first time we need it (there might not be any stuck positioners, eg in tests)
    skybricks = None

    stuck_sky = dict()
    for tile_id, tile_ra, tile_dec, tile_obstime, tile_theta, tile_obsha in zip(
            tiles.id, tiles.ra, tiles.dec, tiles.obstime, tiles.obstheta, tiles.obshourang):
        stuck_sky[tile_id] = dict()

        # Stuck locations and their angles
        # (grab the hw dictionaries once -- these are python wrappers over C++ so not simple accessors)
        state = hw.state
        devtype = hw.loc_device_type
        stuck_loc = [loc for loc in hw.locations
                     if (((state[loc] & FIBER_STATE_STUCK) != 0) and
                         ((state[loc] & FIBER_STATE_BROKEN) == 0) and
                         (devtype[loc] == 'POS'))]
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

        if skybricks is None:
            try:
                skybricks = Skybricks()
            except:
                log.warning('Environment variable SKYBRICKS_DIR is not set; not looking '
                            'up whether stuck positioners land on good sky')
                return

        good_sky = skybricks.lookup_tile(tile_ra, tile_dec, hw.focalplane_radius_deg,
                                         loc_ra, loc_dec)
        log.info('%i of %i stuck positioners land on good sky locations' % (np.sum(good_sky), len(good_sky)))
        for loc,good in zip(stuck_loc, good_sky):
            stuck_sky[tile_id][loc] = good

    return stuck_sky

