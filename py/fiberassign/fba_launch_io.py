"""
fiberassign.fba_launch_io
=============
Utility functions for fba_launch
"""
from __future__ import absolute_import, division

import os
import sys
import subprocess
import numpy as np
from glob import glob
from astropy.io import fits
from astropy.table import Table
import fitsio
import desitarget
from desitarget.io import read_targets_in_tiles, write_targets, write_skies
from desitarget.mtl import inflate_ledger
from desitarget.targetmask import desi_mask, obsconditions
from desitarget.targets import set_obsconditions
import desimodel
from desimodel.footprint import is_point_in_desi
import fiberassign
from fiberassign.scripts.assign import parse_assign, run_assign_full
from fiberassign.assign import merge_results, minimal_target_columns
from time import time
from datetime import datetime, timezone
from astropy.time import Time
import tempfile
import shutil
from fiberassign.utils import Logger


def assert_isoformat_utc(time_str):
    """
    Asserts if a date formats as "YYYY-MM-DDThh:mm:ss+00:00".
    
    Args:
        time_str: string with a date
    Returns:
        boolean asserting if time_str formats as "YYYY-MM-DDThh:mm:ss+00:00"
    """
    try:
        test_time = datetime.strptime(time_str, "%Y-%m-%dT%H:%M:%S%z")
    except ValueError:
        return False
    # AR/SB it parses as an ISO string, now just check UTC timezone +00:00 and not +0000
    return time_str.endswith("+00:00")


def custom_read_targets_in_tiles(
    targdirs,
    tiles,
    quick=True,
    mtl=False,
    unique=True,
    isodate=None,
    log=None,
    step="",
    start=None,
):
    """
    Wrapper to desitarget.io.read_targets_in_tiles, allowing multiple folders
    and sanity check TARGETID if more than one folder.

    Args:
        targdirs: list of folders
        tiles: tiles object (as required by desitarget.io.read_targets_in_tiles)
        quick, mtl, unique, isodate (optional): same as desitarget.io.read_targets_in_tiles arguments
        log (optional): Logger object
        step (optional): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start (optional): start time for log (in seconds; output of time.time()
    Returns:
        array of targets in the passed tiles.
    """
    if log is None:
        log = Logger.get()
    if start is None:
        start = time()
    # AR reading
    log.info(
        "{:.1f}s\t{}\treading input targets from {}".format(
            time() - start, step, targdirs
        )
    )
    if len(targdirs) == 1:
        d = read_targets_in_tiles(
            targdirs[0],
            tiles=tiles,
            quick=quick,
            mtl=mtl,
            unique=unique,
            isodate=isodate,
        )
    else:
        ds = [
            read_targets_in_tiles(
                targdir,
                tiles=tiles,
                quick=quick,
                mtl=mtl,
                unique=unique,
                isodate=isodate,
            )
            for targdir in targdirs
        ]
        # AR merging
        d = np.concatenate(ds)
        # AR remove duplicates based on TARGETID (so duplicates not identified if in mixed surveys)
        ii_m1 = np.where(d["TARGETID"] == -1)[0]
        ii_nm1 = np.where(d["TARGETID"] != -1)[0]
        _, ii = np.unique(d["TARGETID"][ii_nm1], return_index=True)
        ii_nm1 = ii_nm1[ii]
        if len(ii_m1) + len(ii_nm1) != len(d):
            log.info(
                "{:.1f}s\t{}\tremoving {}/{} duplicates".format(
                    time() - start, step, len(d) - len(ii_m1) - len(ii_nm1), len(d)
                )
            )
            d = d[ii_m1.tolist() + ii_nm1.tolist()]
    return d


def mv_write_targets_out(infn, targdir, outfn, log=None, step="", start=None):
    """
    Moves the file created by desitarget.io.write_targets
    and removes folder created by desitarget.io.write_targets    
    
    Args:
        infn: filename output by desitarget.io.write_targets
        targdir: folder provided as desitarget.io.write_targets input
        outfn: desired renaming of infn
        log (optional): Logger object
        step (optional): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start (optional): start time for log (in seconds; output of time.time()
    """
    if log is None:
        log = Logger.get()
    if start is None:
        start = time()
    # AR renaming
    _ = shutil.move(infn, outfn)
    log.info("{:.1f}s\t{}\trenaming {} to {}".format(time() - start, step, infn, outfn))
    # AR removing folders
    if targdir[-1] != "/":
        targdir = "{}/".format(targdir)
    tmpdirs = infn.replace(targdir, "").split("/")[:-1]
    for i in range(len(tmpdirs))[::-1]:
        os.rmdir(os.path.join(*[targdir] + tmpdirs[: i + 1]))


def get_isaen(g, aen):
    """
    Gaia AEN criterion to identify stars
    
    Args:
        g: Gaia GAIA_PHOT_G_MEAN_MAG
        aen: GAIA GAIA_ASTROMETRIC_EXCESS_NOISE
    Returns:
        boolean array, True for identified stars
    """
    return np.logical_or(
        (g <= 19.0) * (aen < 10.0 ** 0.5),
        (g >= 19.0) * (aen < 10.0 ** (0.5 + 0.2 * (g - 19.0))),
    )


def get_nowradec(ra, dec, pmra, pmdec, parallax, ref_year, pmtime_utc_str, scnd=False):
    """
    Apply proper motion correction
    
    Args:
        ra: numpy array of RAs (deg)
        dec: numpy array of DECs (deg)
        pmra: numpy array of projected proper-motion in RA (mas/year)
        pmdec: numpy array of projected proper-motion in DEC (mas/year)
        parallax: numpy array of parallax (mas)
        ref_year: reference epoch (e.g. 2015.5 for Gaia/DR2)
        pmtime_utc_str: date to update position to (format: YYYY-MM-DDThh:mm:ss+00:00)
        scnd (optional, defaults to False): secondary target? (boolean; if True, sets parallax=0)
    Returns:
        ra: numpy array of RAs updated to pmtime_utc_str (deg)
        dec: numpy array of DECs updated to pmtime_utc_str (deg)
    Notes:
        Courtesy of DL; adapted from legacypipe.survey
        Originally named radec_at_mjd()
    """
    # AR pmtime_utc : UTC time of the new ref_epoch; "%Y-%m-%dT%H:%M:%S%z", e.g. "2021-04-21T00:00:00+00:00"
    # AR scnd=True -> parallax is set to 0, i.e. not used
    """
    Units:
    - matches Gaia DR1/DR2
    - pmra,pmdec are in mas/yr.
      pmra is in angular speed (ie, has a cos(dec) factor)
    - parallax is in mas.
    Returns: RA,Dec
    """
    equinox = 53084.28  # mjd of the spring equinox in 2004
    equinox_jyear = Time(equinox, format="mjd").jyear
    axistilt = 23.44  # degrees
    arcsecperrad = 3600.0 * 180.0 / np.pi
    # AR pmtime
    pmtime_utc = datetime.strptime(pmtime_utc_str, "%Y-%m-%dT%H:%M:%S%z")
    pmtime_utc_jyear = Time(pmtime_utc).jyear
    pmtime_utc_mjd = Time(pmtime_utc).mjd

    def xyztoradec(xyz):
        assert len(xyz.shape) == 2
        ra = np.arctan2(xyz[:, 1], xyz[:, 0])  # AR added "np." in front of arctan2...
        ra += 2 * np.pi * (ra < 0)
        norm = np.sqrt(np.sum(xyz ** 2, axis=1))
        dec = np.arcsin(xyz[:, 2] / norm)
        return np.rad2deg(ra), np.rad2deg(dec)

    def radectoxyz(ra_deg, dec_deg):  # AR changed inputs from ra,dec to ra_deg,dec_deg
        ra = np.deg2rad(ra_deg)
        dec = np.deg2rad(dec_deg)
        cosd = np.cos(dec)
        return np.vstack((cosd * np.cos(ra), cosd * np.sin(ra), np.sin(dec))).T

    dt = pmtime_utc_jyear - ref_year
    cosdec = np.cos(np.deg2rad(dec))
    dec = dec + dt * pmdec / (3600.0 * 1000.0)
    ra = ra + (dt * pmra / (3600.0 * 1000.0)) / cosdec
    parallax = np.atleast_1d(parallax)
    # AR discards parallax for scnd=True
    if scnd == True:
        parallax *= 0.0
    I = np.flatnonzero(parallax)
    if len(I):
        suntheta = 2.0 * np.pi * np.fmod(pmtime_utc_jyear - equinox_jyear, 1.0)
        # Finite differences on the unit sphere -- xyztoradec handles
        # points that are not exactly on the surface of the sphere.
        axis = np.deg2rad(axistilt)
        scale = parallax[I] / 1000.0 / arcsecperrad
        xyz = radectoxyz(ra[I], dec[I])
        xyz[:, 0] += scale * np.cos(suntheta)
        xyz[:, 1] += scale * np.sin(suntheta) * np.cos(axis)
        xyz[:, 2] += scale * np.sin(suntheta) * np.sin(axis)
        r, d = xyztoradec(xyz)
        ra[I] = r
        dec[I] = d
    return ra, dec


def force_finite_pm(
    d, pmra_key="PMRA", pmdec_key="PMDEC", log=None, step="", start=None
):
    """
    Replaces NaN PMRA, PMDEC by 0    
    
    Args:
        d: array with at least proper-motion columns
        pmra_key (optional, defaults to PMRA): column name for PMRA
        pmdec_key (optional, defaults to PMDEC): column name for PMDEC
        log (optional): Logger object
        step (optional): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start (optional): start time for log (in seconds; output of time.time()

    Returns:
        d: same as input d, but NaN proper motions replaced by 0
        
    """
    if log is None:
        log = Logger.get()
    if start is None:
        start = time()

    for key in [pmra_key, pmdec_key]:
        keep = ~np.isfinite(d[key])
        if keep.sum() > 0:
            d[key][keep] = 0.0
            log.info(
                "{:.1f}s\t{}\t replacing NaN by 0 for {} targets".format(
                    time() - start, step, keep.sum()
                )
            )
    return d


def force_nonzero_refepoch(
    d,
    force_ref_epoch,
    ref_epoch_key="REF_EPOCH",
    pmra_key="PMRA",
    pmdec_key="PMDEC",
    log=None,
    step="",
    start=None,
):
    """
    Replaces 0 by force_ref_epoch in ref_epoch
    
    Args:
        d: array with at least proper-motion columns
        force_ref_epoch: float, ref_epoch to replace 0 by
        ref_epoch_key (optional, defaults to REF_EPOCH): column name for the ref_epoch
        pmra_key (optional, defaults to PMRA): column name for PMRA
        pmdec_key (optional, defaults to PMDEC): column name for PMDEC
        log (optional): Logger object
        step (optional): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start (optional): start time for log (in seconds; output of time.time()        
    Returns:
        d: same as input d, but 0 ref_epochs replaced by force_ref_epoch

    Notes:
        Will exit with error if ref_epoch=0, but pmra or pmdec != 0
        
    """
    if log is None:
        log = Logger.get()
    if start is None:
        start = time()
    keep = d[ref_epoch_key] == 0
    n = (d[pmra_key][keep] != 0) | (d[pmra_key][keep] != 0)
    if n > 0:
        log.error(
            "{:.1f}s\t{}\t{} targets have {}=0 but {} or {} != 0; exiting".format(
                time() - start, step, n, ref_epoch_key, pmra_key, pmdec_key,
            )
        )
        sys.exit(1)
    d[ref_epoch_key][keep] = ref_epochs[force_ref_epoch]
    log.info(
        "{:.1f}s\t{}\tsetting {}={} for {} objects with {}=0".format(
            time() - start,
            step,
            ref_epoch_key,
            force_ref_epoch,
            keep.sum(),
            ref_epoch_key,
        )
    )
    return d


def update_nowradec(
    d,
    pmtime_utc_str,
    ra_key="RA",
    dec_key="DEC",
    pmra_key="PMRA",
    pmdec_key="PMDEC",
    parallax_key="PARALLAX",
    ref_epoch_key="REF_EPOCH",
    gaiag_key="GAIA_PHOT_G_MEAN_MAG",
    gaiaaen_key="GAIA_ASTROMETRIC_EXCESS_NOISE",
    scnd=False,
    log=None,
    step="",
    start=None,
):
    """
    Update (RA, DEC, REF_EPOCH) using proper motion
    
    Args:
        d: array with at least proper-motion columns
        pmtime_utc_str: date to update position to (format: YYYY-MM-DDThh:mm:ss+00:00)
        ra_key (optional, defaults to RA): column name for RA
        dec_key (optional, defaults to DEC): column name for DEC
        pmra_key (optional, defaults to PMRA): column name for PMRA
        pmdec_key (optional, defaults to PMDEC): column name for PMDEC
        parallax_key (optional, defaults to PARALLAX): column name for PARALLAX
        ref_epoch_key (optional, defaults to REF_EPOCH): column name for the REF_EPOCH
        gaia_key (optional, defaults to GAIA_PHOT_G_MEAN_MAG): column name for Gaia g-mag
        gaiaaen_key (optional, defaults to GAIA_ASTROMETRIC_EXCESS_NOISE): column name for Gaia GAIA_ASTROMETRIC_EXCESS_NOISE
        scnd (optional, defaults to False): secondary target? (boolean);
              if False, update for REF_EPOCH>0 + AEN only
              if True, update for REF_EPOCH>0 + finite(PMRA,PMDEC) ; forces PARALLAX=0
        log (optional): Logger object
        step (optional): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start (optional): start time for log (in seconds; output of time.time()        
        
    Returns:
        d: same as input, but with RA, DEC updated to pmtime_utc_str

    Notes:
        REF_EPOCH is updated for *all* objects
    """
    if log is None:
        log = Logger.get()
    if start is None:
        start = time()
    # AR
    pmtime_utc = datetime.strptime(pmtime_utc_str, "%Y-%m-%dT%H:%M:%S%z")
    pmtime_utc_jyear = Time(pmtime_utc).jyear
    # AR computing positions at pmtime_utc_str using Gaia PMRA, PMDEC
    nowra, nowdec = get_nowradec(
        d[ra_key],
        d[dec_key],
        d[pmra_key],
        d[pmdec_key],
        d[parallax_key],
        d[ref_epoch_key],
        pmtime_utc_str,
        scnd=scnd,
    )
    if scnd == True:
        # AR secondary: REF_EPOCH>0
        keep = d["REF_EPOCH"] > 0
    else:
        # AR targets with REF_EPOCH>0 and passing the AEN criterion
        keep = (d["REF_EPOCH"] > 0) & (get_isaen(d[gaiag_key], d[gaiaaen_key]))
    # AR storing changes to report extrema in the log
    dra = nowra - d[ra_key]
    ddec = nowdec - d[dec_key]
    # AR updating positions to pmtime_utc_str for targets passing the AEN criterion
    d[ra_key][keep] = nowra[keep]
    d[dec_key][keep] = nowdec[keep]
    log.info(
        "{:.1f}s\t{}\tupdating RA,DEC at {} with PM for {:.0f}/{:.0f} targets passing AEN; maximum changes: RA={:.1f},{:.1f} arcsec, DEC={:.1f},{:.1f} arcsec".format(
            time() - start,
            step,
            pmtime_utc_jyear,
            keep.sum(),
            len(keep),
            3600.0 * dra.min(),
            3600.0 * dra.max(),
            3600 * ddec.min(),
            3600.0 * ddec.max(),
        )
    )
    # AR updating REF_EPOCH for *all* objects (for PlateMaker)
    d[ref_epoch_key] = pmtime_utc_jyear
    log.info(
        "{:.1f}s\tupdating REF_EPOCH to {} for all {} targets".format(
            time() - start, pmtime_utc_jyear, len(keep)
        )
    )
    return d
