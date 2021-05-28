# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.gfa
===============

Functions for computing coverage of GFAs
"""

import numpy as np
import fitsio
from astropy.table import Table
import desimodel.focalplane.gfa

from ._internal import Tiles
from .utils import Logger, Timer


def isolated(ra, dec, mindist=5.0):
    '''
    Returns bool array for whether targets are isolated
    Args:
        ra: array, RA in degrees
        dec: array, dec in degrees
        mindist: minimum separation distance in arcsec
    '''
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    cx = SkyCoord(ra*u.deg, dec*u.deg)
    index, separation, dist3d = cx.match_to_catalog_sky(cx, nthneighbor=2)
    return separation.to(u.arcsec) > mindist*u.arcsec

def is_gaia_variable(tab):
    '''
    Decide whether each star is flux variable based on Gaia variability proxy
    Args:
        tab: table with columns GAIA_PHOT_G_MEAN_MAG, GAIA_PHOT_G_N_OBS,
             GAIA_PHOT_G_MEAN_FLUX_OVER_ERROR
    Returns:
        boolean array listing whether each star is variable (True) or
        not (False)

    The Gaia flux variability proxy can only be computed when the column
    GAIA_PHOT_G_N_OBS is available. This column first became available
    in the desitarget 'gfas' files when these were transitioned to
    use Gaia eDR3 (even though the column PHOT_G_N_OBS was not a new eDR3
    addition to the native Gaia catalogs themselves).

    If any of the columns needed to compute the variability flag are not
    present, then all stars will be labeled as not variable.

    One reference for the variability metric and cut adopted here is:
        https://arxiv.org/pdf/2009.07746.pdf
    Specifically Equation 2.
    '''

    colnames = tab.colnames

    if ('GAIA_PHOT_G_MEAN_MAG' not in colnames) or ('GAIA_PHOT_G_N_OBS' not in colnames) or ('GAIA_PHOT_G_MEAN_FLUX_OVER_ERROR' not in colnames):
        return np.zeros(len(tab), dtype=bool)

    # avoid division by zero in cases with REF_CAT = T2, where
    # GAIA_PHOT_G_N_OBS and GAIA_PHOT_G_MEAN_FLUX_OVER_ERROR will both be
    # 0-valued in the desitarget 'gfas' files
    # in cases with REF_CAT = T2, this setup will give amplproxyg = 0
    amplproxyg = np.sqrt(tab['GAIA_PHOT_G_N_OBS'])/(tab['GAIA_PHOT_G_MEAN_FLUX_OVER_ERROR'] + (tab['GAIA_PHOT_G_MEAN_FLUX_OVER_ERROR'] == 0).astype(float))

    is_variable = np.logical_and(amplproxyg > 0.06, tab['GAIA_PHOT_G_MEAN_MAG'] < 18.3)

    return is_variable

def gaia_synth_r_flux(tab):
    '''
    Generate synthetic r band flux based on Gaia photometry
    Args:
        tab: table with columns GAIA_PHOT_G_MEAN_MAG, GAIA_PHOT_BP_MEAN_MAG,
             GAIA_PHOT_RP_MEAN_MAG, REF_CAT

    Returns:
        array of synthetic Gaia-based r-band fluxes

    This code has largely been borrowed from legacypipe/reference.py
    Uses Rongpu Zhou's transformations from (G, BP-RP) to DECam r
    '''

    gaiadr = "dr2" if "G2" in tab["REF_CAT"] else "edr3"

    color = tab["GAIA_PHOT_BP_MEAN_MAG"] - tab["GAIA_PHOT_RP_MEAN_MAG"]

    # the 1.0.0 targeting gfas files have a mixture of zero and NaN
    # placeholder values for BP and RP
    badcolor = np.logical_not(np.isfinite(color)) | (tab["GAIA_PHOT_BP_MEAN_MAG"] == 0) | (tab["GAIA_PHOT_RP_MEAN_MAG"] == 0)

    # color clipping parameters to match transform training color range
    # from Rongpu:
    color_lo = {"dr2": -0.6, "edr3": -0.5}
    color_hi = {"dr2": 4.1, "edr3": 4.5}

    # clip to reasonable range for the polynomial fit
    color = np.clip(color, color_lo[gaiadr], color_hi[gaiadr])

    coeffs = {"dr2": [0.1139078673, -0.2868955307, 0.0013196434, 0.1029151074,
                      0.1196710702, -0.3729031390, 0.1859874242, 0.1370162451,
                     -0.1808580848, 0.0803219195, -0.0180218196, 0.0020584707,
                     -0.0000953486],
              "edr3": [0.1431278873, -0.2999797766, -0.0553379742, 0.1544273115,
                       0.3068634689, -0.9499143903, 0.9769739362, -0.4926704528,
                       0.1272539574, -0.0133178183, -0.0008153813, 0.0003094116,
                      -0.0000198891]}

    coeffs = coeffs[gaiadr]

    synth_r_mag = np.array(tab["GAIA_PHOT_G_MEAN_MAG"])

    for order,c in enumerate(coeffs):
            synth_r_mag += c * color**order

    synth_r_flux = np.power(10, (22.5-synth_r_mag)/2.5)

    # -99 is the FLUX_R placeholder that GFA_TARGETS uses
    synth_r_flux[badcolor] = -99.0

    return synth_r_flux

def get_gfa_targets(tiles, gfafile, faintlim=99):
     
    """Returns a list of tables of GFA targets on each tile

    Args:
        tiles: table with columns TILEID, RA, DEC; or Tiles object
        targets: table of targets with columsn RA, DEC

    Returns:
        list of tables (one row per input tile) with the subset of targets
        that are covered by GFAs on each tile.  Each table has additional
        `GFA_LOC` column indicating 0-9 which GFA was covered.

    Note that a given target could be covered by GFAs on more than one tile.

    Output is a list of astropy Tables; inputs can be numpy structured arrays
    or astropy Tables
    """
    log = Logger.get()
    tm = Timer()
    tm.start()

    # Convert tiles to vanilla numpy array if needed
    if isinstance(tiles, Tiles):
        tx = np.zeros(len(tiles.ra),
                      dtype=[("RA", "f8"), ("DEC", "f8"), ("TILEID", "i4")])
        tx["RA"] = tiles.ra
        tx["DEC"] = tiles.dec
        tx["TILEID"] = tiles.id
        tiles = tx

    # Load potential GFA targets and GFA locations
    targets = fitsio.read(gfafile)
    gfa = desimodel.focalplane.gfa.GFALocations(scale=2)

    # Pre-filter what GFA targets cover what tiles with some buffer.
    # find_points_in_tiles returns a list of lists;
    # convert to dictionary of lists keyed by tileid
    log.info("Finding overlap of {} GFA targets on {} tiles".format(
        len(targets), len(tiles)))
    gfa_tile_indices = dict()
    ii = desimodel.footprint.find_points_in_tiles(
        tiles, targets["RA"], targets["DEC"], radius=1.8)
    for i, tileid in enumerate(tiles["TILEID"]):
        gfa_tile_indices[tileid] = ii[i]

    gfa_targets = list()

    log.info("Generating GFA targets tables")
    for telra, teldec, tileid in zip(tiles["RA"], tiles["DEC"],
                                     tiles["TILEID"]):
        tmp = gfa.targets_on_gfa(telra, teldec,
                                 targets[gfa_tile_indices[tileid]])
        t = Table(tmp)
        
        # Rename some columns for downstream clarity and consistency
        for oldname, newname in [
                ("TYPE", "MORPHTYPE"),
                ("RA", "TARGET_RA"),
                ("DEC", "TARGET_DEC"),
                ("RA_IVAR", "TARGET_RA_IVAR"),
                ("DEC_IVAR", "TARGET_DEC_IVAR")]:
            if oldname in t.colnames:
                t.rename_column(oldname, newname)

        # Select which targets are good for ETC / GUIDE / FOCUS
        # 0 == good

        flag = np.zeros(len(t), dtype="i2")

        #- Not PSF-like
        isPSF = (t["MORPHTYPE"] == "PSF ") | (t["MORPHTYPE"] == "GPSF") | (t["MORPHTYPE"] == "PSF")
        flag[~isPSF] |= 2**0

        #- Not Isolated
        if len(tmp) > 1:
            notIsolated = ~isolated(tmp['RA'], tmp['DEC'])
            flag[notIsolated] |= 2**1

        #- Questionable astrometry / proper motion
        tych = (0 < t['REF_ID']) 
        tych &= ( t['REF_ID'] < 1e10)
        flag[tych] |= 2**2

        #- Too faint
        faint = t['GAIA_PHOT_G_MEAN_MAG'] > faintlim
        flag[faint] |= 2**3

        # AR not passing the Gaia AEN criterion (PM correction done for AEN targets only)
        g = t["GAIA_PHOT_G_MEAN_MAG"]
        aen = t["GAIA_ASTROMETRIC_EXCESS_NOISE"]
        isaen = np.logical_or(
            (g <= 19.0) * (aen < 10.0 ** 0.5),
            (g >= 19.0) * (aen < 10.0 ** (0.5 + 0.2 * (g - 19.0))),
        )
        flag[~isaen] |= 2**4
        
        
        if len(flag)-np.count_nonzero(flag) == 0:
            log.error("ERROR: no good GFA targets for "
                      "ETC/GUIDE/FOCUS on tile {}".format(tileid))

        t["GUIDE_FLAG"] = flag
        t["FOCUS_FLAG"] = flag

        # variable star flagging specifically for ETC
        is_variable = is_gaia_variable(t)
        flag[is_variable] |= 2**5

        t["ETC_FLAG"] = flag

        # patch in Gaia-based synthetic r flux for use by ETC
        t["FLUX_R"] = gaia_synth_r_flux(t)

        gfa_targets.append(t)

    tm.stop()
    tm.report("  Identifying GFA targets")

    return gfa_targets
