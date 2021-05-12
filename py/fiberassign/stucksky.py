import os
import numpy as np

class Skybricks(object):
    '''
    This class handles dynamic lookup of whether a given (RA,Dec) should make a good location for a sky fiber.
    '''    
    def __init__(self, skybricks_dir=None):
        '''
        Create a Skybricks object, reading metadata.

        *skybricks_dir*: directory to find skybricks data files; if None, will read from SKYBRICKS_DIR environment variable.
        '''
        import fitsio
        if skybricks_dir is None:
            skybricks_dir = os.environ.get('SKYBRICKS_DIR', None)
        if skybricks_dir is None:
            raise RuntimeError('Environment variable SKYBRICKS_DIR is not set; needed to look up dynamic sky fiber positions')
        self.skybricks_dir = skybricks_dir
        skybricks_fn = os.path.join(self.skybricks_dir, 'skybricks-exist.fits')
        self.skybricks = fitsio.read(skybricks_fn, upper=True)
        self.skykd = _radec2kd(self.skybricks['RA'], self.skybricks['DEC'])

    def lookup_tile(self, tilera, tiledec, tileradius,
                    ras, decs):
        '''
        Given scalar (tilera, tiledec) RA,Dec (in degrees) and tile radius (in degrees),
        for a set of RA,Dec points within that tile, does the lookup for all those positions,
        returning a *good_sky* boolean array.

        *tilera*: scalar RA in deg of tile center
        *tiledec*: scalar Dec in deg of tile center
        *tileradius*: scalar radiun in deg of the tile
        *ras*: numpy array of RAs in deg
        *decs*: numpy array of Decs in deg

        Returns:
        *good_sky*: boolean numpy array of length(ras), would these RA,Decs make good sky fibers?
        '''
        import fitsio
        from astropy.wcs import WCS
        from fiberassign.utils import Logger

        log = Logger.get()
        # skybricks are 1 x 1 deg.
        brickrad = (1. * np.sqrt(2.) / 2.)
        searchrad = 1.01 * (tileradius + brickrad)
        # here, convert search radius to radians -- an overestimate vs
        # unit-sphere distance, but that's the safe direction.
        searchrad = np.deg2rad(searchrad)
        tilexyz = _radec2xyz([tilera], [tiledec])
        sky_inds = self.skykd.query_ball_point(tilexyz[0,:], searchrad)
        good_sky = np.zeros(len(ras), bool)
        # Check possibly-overlapping skybricks.
        for i in sky_inds:
            # Do any of the query points overlap in the brick's RA,DEC
            # unique-area bounding-box?
            I = np.flatnonzero(
                (ras  >= self.skybricks['RA1'][i]) *
                (ras  <  self.skybricks['RA2'][i]) *
                (decs >= self.skybricks['DEC1'][i]) *
                (decs <  self.skybricks['DEC2'][i]))
            log.debug('Skybricks: %i locations overlap skybrick %s' %
                      (len(I), self.skybricks['BRICKNAME'][i]))
            if len(I) == 0:
                continue

            # Read skybrick file
            fn = os.path.join(self.skybricks_dir,
                              'sky-%s.fits.fz' % self.skybricks['BRICKNAME'][i])
            if not os.path.exists(fn):
                log.warning('Missing "skybrick" file: %s' % fn)
                continue
            skymap,hdr = fitsio.read(fn, header=True)
            H,W = skymap.shape
            # create WCS object
            w = WCS(naxis=2)
            w.wcs.ctype = [hdr['CTYPE1'], hdr['CTYPE2']]
            w.wcs.crpix = [hdr['CRPIX1'], hdr['CRPIX2']]
            w.wcs.crval = [hdr['CRVAL1'], hdr['CRVAL2']]
            w.wcs.cd = [[hdr['CD1_1'], hdr['CD1_2']],
                        [hdr['CD2_1'], hdr['CD2_2']]]
            x,y = w.wcs_world2pix(ras[I], decs[I], 0)
            x = np.round(x).astype(int)
            y = np.round(y).astype(int)
            # we have margins that should ensure this...
            if not (np.all(x >= 0) and np.all(x <  W) and
                    np.all(y >= 0) and np.all(y <  H)):
                raise RuntimeError('Skybrick %s: locations project outside the brick bounds' %
                                   (self.skybricks['BRICKNAME'][i]))

            # FIXME -- look at surrounding pixels too??
            good_sky[I] = (skymap[y, x] == 0)
        return good_sky


def stuck_on_sky(hw, tiles):
    '''
    Will STUCK positioners land on good SKY locations for the given set of tiles?

    Returns a nested dict:
        stuck_sky[tileid][loc] = bool_good_sky
    '''
    from fiberassign.utils import Logger
    from fiberassign.hardware import FIBER_STATE_STUCK, FIBER_STATE_BROKEN

    from fiberassign.targets import xy2radec

    log = Logger.get()

    try:
        skybricks = Skybricks()
    except:
        log.warning('Environment variable SKYBRICKS_DIR is not set; not looking '
                    'up whether stuck positioners land on good sky')
        return

    stuck_sky = dict()
    for tile_id, tile_ra, tile_dec, tile_obstime, tile_theta, tile_obsha in zip(
            tiles.id, tiles.ra, tiles.dec, tiles.obstime, tiles.obstheta, tiles.obshourang):
        stuck_sky[tile_id] = dict()

        # Stuck locations and their angles
        stuck_loc = [loc for loc in hw.locations
                     if (((hw.state[loc] & FIBER_STATE_STUCK) != 0) and
                         ((hw.state[loc] & FIBER_STATE_BROKEN) == 0) and
                         (hw.loc_device_type[loc] == 'POS'))]
        stuck_theta = [hw.loc_theta_pos[loc] + hw.loc_theta_offset[loc] for loc in stuck_loc]
        stuck_phi   = [hw.loc_phi_pos  [loc] + hw.loc_phi_offset  [loc] for loc in stuck_loc]

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
                hw.loc_phi_max[loc],
                True
            )
            stuck_x[iloc] = loc_x
            stuck_y[iloc] = loc_y

        loc_ra,loc_dec = xy2radec(
            hw, tile_ra, tile_dec, tile_obstime, tile_theta, tile_obsha,
            stuck_x, stuck_y, False, 0
        )

        good_sky = skybricks.lookup_tile(tile_ra, tile_dec, hw.focalplane_radius_deg,
                                         loc_ra, loc_dec)
        log.info('%i of %i stuck positioners land on good sky locations' % (np.sum(good_sky), len(good_sky)))
        for loc,good in zip(stuck_loc, good_sky):
            stuck_sky[tile_id][loc] = good

    return stuck_sky

def _radec2kd(ra, dec):
    from scipy.spatial import KDTree
    xyz = _radec2xyz(ra, dec)
    return KDTree(xyz)

def _radec2xyz(ra, dec):
    rr = np.deg2rad(ra)
    dd = np.deg2rad(dec)
    xyz = np.zeros((len(rr),3))
    xyz[:,0] = np.cos(rr) * np.cos(dd)
    xyz[:,1] = np.sin(rr) * np.cos(dd)
    xyz[:,2] = np.sin(dd)
    return xyz
