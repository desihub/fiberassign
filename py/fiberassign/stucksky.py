    
def stuck_on_sky(hw, tiles):
    '''
    Will STUCK positioners land on good SKY locations for the given set of tiles?

    Returns a nested dict:
        stuck_sky[tileid][loc] = bool_good_sky
    '''
    import os
    import numpy as np
    from scipy.spatial import KDTree
    import fitsio
    from astropy.wcs import WCS
    from fiberassign.hardware import FIBER_STATE_STUCK
    from fiberassign.utils import Logger

    skybricks_dir = os.environ.get('SKYBRICKS_DIR', None)
    if skybricks_dir is None:
        log = Logger.get()
        log.warning('Environment variable SKYBRICKS_DIR is not set; not looking up whether '
                    'stuck positioners land on good sky')
        return
    skybricks_fn = os.path.join(skybricks_dir, 'skybricks-exist.fits')
    skybricks = fitsio.read(skybricks_fn, upper=True)
    skykd = _radec2kd(skybricks['RA'], skybricks['DEC'])

    tilekd = _radec2kd(tiles.ra, tiles.dec)

    # skybricks are 1 x 1 deg.
    brickrad = (1. * np.sqrt(2.) / 2.)
    searchrad = 1.01 * (hw.focalplane_radius_deg + brickrad)

    # here, convert search radius to radians -- an overestimate vs
    # unit-sphere distance, but that's the safe direction.
    sky_indices = tilekd.query_ball_tree(skykd, np.deg2rad(searchrad))

    stuck_sky = dict()
    for tile_id, tile_ra, tile_dec, tile_theta, sky_inds in zip(
            tiles.id, tiles.ra, tiles.dec, tiles.obstheta, sky_indices):
        stuck_sky[tile_id] = dict()
    
        # Stuck locations and their angles
        stuck_loc = [x for x in hw.locations if hw.state[x] & FIBER_STATE_STUCK]
        stuck_theta = [hw.loc_theta_pos[x] for x in stuck_loc]
        stuck_phi   = [hw.loc_phi_pos  [x] for x in stuck_loc]

        # Convert positioner angle orientations to curved focal surface X / Y (not CS5)
        # Note:  we could add some methods to the python bindings to vectorize this or make it less clunky...
        stuck_x = np.zeros(len(stuck_loc))
        stuck_y = np.zeros(len(stuck_loc))
        for iloc, (loc, theta, phi) in enumerate(zip(stuck_loc, stuck_theta, stuck_phi)):
            loc_x, loc_y = hw.thetaphi_to_xy(
                hw.loc_pos_curved_mm[loc],
                theta,
                phi,
                hw.loc_theta_arm[loc],
                hw.loc_phi_arm[loc],
                hw.loc_theta_offset[loc],
                hw.loc_phi_offset[loc],
                hw.loc_theta_min[loc],
                hw.loc_phi_min[loc],
                hw.loc_theta_max[loc],
                hw.loc_phi_max[loc]
            )
            stuck_x[iloc] = loc_x
            stuck_y[iloc] = loc_y
        
        loc_radec = hw.xy2radec_multi(
            tile_ra,
            tile_dec,
            tile_theta,
            stuck_x,
            stuck_y,
            False,
            0
        )
        loc_ra  = np.array([r for r,d in loc_radec])
        loc_dec = np.array([d for r,d in loc_radec])

        print('Tile', tile_ra, tile_dec)
        
        good_sky = np.zeros(len(loc_ra), bool)
        # Check possibly-overlapping skybricks.
        for i in sky_inds:
            loc_in = np.flatnonzero(
                (loc_ra  >= skybricks['RA1'][i]) *
                (loc_ra  <  skybricks['RA2'][i]) *
                (loc_dec >= skybricks['DEC1'][i]) *
                (loc_dec <  skybricks['DEC2'][i]))
            print(len(loc_in), 'fibers overlap brick', skybricks['BRICKNAME'][i])
            if len(loc_in) == 0:
                continue

            fn = os.path.join(skybricks_dir,
                              'sky-%s.fits.gz' % skybricks['BRICKNAME'][i])
            # FIXME -- debugging
            if not os.path.exists(fn):
                print('Missing:', fn)
                continue

            skymap,hdr = fitsio.read(fn, header=True)
            H,W = skymap.shape

            # load WCS
            w = WCS(naxis=2)
            w.wcs.ctype = [hdr['CTYPE1'], hdr['CTYPE2']]
            w.wcs.crpix = [hdr['CRPIX1'], hdr['CRPIX2']]
            w.wcs.crval = [hdr['CRVAL1'], hdr['CRVAL2']]
            w.wcs.cd = [[hdr['CD1_1'], hdr['CD1_2']],
                        [hdr['CD2_1'], hdr['CD2_2']]]
            x,y = w.wcs_world2pix(loc_ra[loc_in], loc_dec[loc_in], 0)
            x = np.round(x).astype(int)
            y = np.round(y).astype(int)
            # we have margins that should ensure this...
            assert(np.all(x >= 0))
            assert(np.all(x <  W))
            assert(np.all(y >= 0))
            assert(np.all(y <  H))

            # FIXME -- look at surrounding pixels too??
            good_sky[loc_in] = (skymap[y, x] == 0)

        print(np.sum(good_sky), 'stuck positioners land on good sky')

        for loc,good in zip(stuck_loc, good_sky):
            stuck_sky[tile_id][loc] = good

    return stuck_sky

def _radec2kd(ra, dec):
    from scipy.spatial import KDTree
    import numpy as np

    rr = np.deg2rad(ra)
    dd = np.deg2rad(dec)
    xyz = np.zeros((len(rr),3))
    xyz[:,0] = np.cos(rr) * np.cos(dd)
    xyz[:,1] = np.sin(rr) * np.cos(dd)
    xyz[:,2] = np.sin(dd)
    return KDTree(xyz)
                
            
if __name__ == '__main__':
    from fiberassign.tiles import load_tiles
    from fiberassign.hardware import load_hardware

    rundate = '2021-04-22T00:00:00'
    tiles = load_tiles('tiles.fits')
    hw = load_hardware(rundate=rundate)

    stuck_on_sky(hw, tiles)
