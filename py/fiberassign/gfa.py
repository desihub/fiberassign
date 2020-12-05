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


def get_gfa_targets(tiles, gfafile,faintlim=99):
     
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

        
        t["ETC_FLAG"] = flag
        t["GUIDE_FLAG"] = flag
        t["FOCUS_FLAG"] = flag

        gfa_targets.append(t)

    tm.stop()
    tm.report("  Identifying GFA targets")

    return gfa_targets
