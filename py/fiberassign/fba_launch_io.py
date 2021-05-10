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
from desitarget.gaiamatch import gaia_psflike
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


# AR default REF_EPOCH for PMRA=PMDEC=REF_EPOCH=0 objects
gaia_ref_epochs = {"dr2": 2015.5}


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
    n = ((d[pmra_key][keep] != 0) | (d[pmra_key][keep] != 0)).sum()
    if n > 0:
        log.error(
            "{:.1f}s\t{}\t{} targets have {}=0 but {} or {} != 0; exiting".format(
                time() - start, step, n, ref_epoch_key, pmra_key, pmdec_key,
            )
        )
        sys.exit(1)
    d[ref_epoch_key][keep] = force_ref_epoch
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
    gaiadr,
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
        gaiadr: Gaia dr ("dr2" or "edr3")
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
        keep = d["REF_EPOCH"] > 0
        # AR gaia_psflike arguments changed at desitarget-0.58.0
        if desitarget.__version__ < "0.58.0":
            keep &= gaia_psflike(d[gaiag_key], d[gaiaaen_key])
        else:
            keep &= gaia_psflike(d[gaiag_key], d[gaiaaen_key], dr=gaiadr)
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


def assert_env_vars(
    required_env_vars=[
        "DESI_ROOT",
        "DESI_TARGET",
        "DESIMODEL",
        "DESI_SURVEYOPS",
        "SKYBRICKS_DIR",
    ],
    log=None,
    step="settings",
    start=None,
):
    """
    Assert the environment variables required by fba_launch 
    
    Args:
        required_env_vars (optional, defaults to ["DESI_ROOT",
        "DESI_TARGET",
        "DESIMODEL",
        "DESI_SURVEYOPS",
        "SKYBRICKS_DIR",]): list of environment variables required by fba_launch
        log (optional): Logger object
        step (optional): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start (optional): start time for log (in seconds; output of time.time()        
        
    Notes:
        will exit with error if some assertions are not verified 
    """
    if log is None:
        log = Logger.get()
    if start is None:
        start = time()

    # AR safe: DESI environment variables
    for required_env_var in required_env_vars:
        if os.getenv(required_env_var) is None:
            log.error(
                "{:.1f}s\t{}\tenvironment variable {} not defined; exiting".format(
                    time() - start, step, required_env_var
                )
            )
            sys.exit(1)


def assert_arg_dates(
    args,
    dates=["pmtime_utc_str", "rundate", "mtltime"],
    log=None,
    step="settings",
    start=None,
):
    """
    Assert the fba_launch date arguments are correctly formatted ("YYYY-MM-DDThh:mm:ss+00:00")
    
    Args:
        args: fba_launch parser.parse_args() output
        dates (optional, defaults to ["pmtime_utc_str", "rundate", "mtltime"]): list of date fba_launch argument names to check
        log (optional): Logger object
        step (optional): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start (optional): start time for log (in seconds; output of time.time()        
        
    Notes:
        will exit with error if some assertions are not verified 
    """
    if log is None:
        log = Logger.get()
    if start is None:
        start = time()

    # AR dates properly formatted?
    for kwargs in args._get_kwargs():
        if kwargs[0] in dates:
            if not assert_isoformat_utc(kwargs[1]):
                log.error(
                    "{:.1f}s\t{}\t{}={} is not yyyy-mm-ddThh:mm:ss+00:00; exiting".format(
                        time() - start, step, kwargs[0], kwargs[1],
                    )
                )
                sys.exit(1)


def assert_svn_tileid(
    tileid, forcetileid="n", log=None, step="settings", start=None,
):
    """
    Asserts if TILEID already exists in the SVN tile folder
    
    Args:
        tileid: TILEID to check (int)
        forcetileid (optional, defaults to "n"): "y" or "n";
                if "n", will trigger a warning + an error
                if "y", e will trigger a warning only
        log (optional): Logger object
        step (optional): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start (optional): start time for log (in seconds; output of time.time()        
       
    """
    if log is None:
        log = Logger.get()
    if start is None:
        start = time()

    svn_trunk = os.path.join(os.getenv("DESI_TARGET"), "fiberassign/tiles/trunk")
    # AR needs a wildcard to verify .fits and fits.gz files
    # AR as the gzipping was not done before ~SV1
    svn_fn = os.path.join(
        svn_trunk, "{:06d}".format(tileid)[:3], "fiberassign-{:06d}.fits".format(tileid)
    )
    if os.path.isfile(svn_fn) | os.path.isfile("{}.gz".format(svn_fn)):
        log.warning(
            "{:.1f}s\t{}\tTILEID={} already exists in SVN folder {}".format(
                time() - start, step, tileid, svn_trunk
            )
        )
        if forcetileid == "y":
            log.warning(
                "{:.1f}s\t{}\tproceeding as forcetileid == y".format(
                    time() - start, step
                )
            )
        else:
            log.error(
                "{:.1f}s\tsettings\texiting as forcetileid == n".format(time() - start)
            )
            sys.exit(1)
    else:
        log.info(
            "{:.1f}s\t{}\tTILEID={} does not exist in SVN folder {}; proceeding".format(
                time() - start, step, tileid, svn_trunk
            )
        )


def print_config_infos(
    required_env_vars=[
        "DESI_ROOT",
        "DESI_TARGET",
        "DESIMODEL",
        "DESI_SURVEYOPS",
        "SKYBRICKS_DIR",
    ],
    log=None,
    step="settings",
    start=None,
):
    """
    Print various configuration informations (machine, modules version/path, DESI environment variables).
    
    Args:
        required_env_vars (optional, defaults to ["DESI_ROOT",
        "DESI_TARGET",
        "DESIMODEL",
        "DESI_SURVEYOPS",
        "SKYBRICKS_DIR",]): list of environment variables required by fba_launch
        log (optional): Logger object
        step (optional): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start (optional): start time for log (in seconds; output of time.time()        
    """
    if log is None:
        log = Logger.get()
    if start is None:
        start = time()

    # AR machine
    log.info(
        "{:.1f}s\t{}\tHOSTNAME={}".format(time() - start, step, os.getenv("HOSTNAME"))
    )

    # AR fiberassign, desitarget, desimodel code version/path
    for module, name in zip(
        [fiberassign, desitarget, desimodel], ["fiberassign", "desitarget", "desimodel"]
    ):
        log.info(
            "{:.1f}s\t{}\trunning with {} code version: {}".format(
                time() - start, step, name, module.__version__
            )
        )
        log.info(
            "{:.1f}s\t{}\trunning with {} code path: {}".format(
                time() - start, step, name, module.__path__
            )
        )

    # AR DESI environment variables
    for required_env_var in required_env_vars:
        log.info(
            "{:.1f}s\t{}\t{}={}".format(
                time() - start, step, required_env_var, os.getenv(required_env_var)
            )
        )


def get_desitarget_paths(
    dtver,
    survey,
    program,
    dr="dr9",
    gaiadr="gaiadr2",
    log=None,
    step="settings",
    start=None,
):
    """
    Obtain the folder/file full paths for desitarget products
    
    Args:
        dtver: desitarget catalog version (string; e.g., "0.57.0")
        survey: survey (string; e.g. "sv1", "sv2", "sv3", "main")
        program: "dark", "bright", or "backup" (string)
        dr (optional, defaults to "dr9"): legacypipe dr (string)
        gaiadr (optional, defaults to "gaiadr2"): gaia dr (string)
        log (optional): Logger object
        step (optional): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start (optional): start time for log (in seconds; output of time.time()        
        
    Returns:
        Dictionary with the following keys:
        - sky: sky folder
        - skysupp: skysupp folder
        - gfa: GFA folder
        - targ: targets folder (static catalogs, with all columns)
        - mtl: MTL folder
        - scnd: secondary fits catalog (static)
        - scndmtl: MTL folder for secondary targets
        - too: ToO ecsv catalog

    Notes:
        if survey not in ["sv1", "sv2", "sv3", "main"]
        or program not in ["dark", "bright", or "backup"], will return a warning only
        same warning only if the built paths/files do not exist.
    """
    if log is None:
        log = Logger.get()
    if start is None:
        start = time()

    # AR expected survey, program?
    exp_surveys = ["sv1", "sv2", "sv3", "main"]
    exp_programs = ["dark", "bright", "backup"]
    if survey.lower() not in exp_surveys:
        log.warning(
            "{:.1f}s\t{}\tunexpected survey={} ({}; proceeding anyway)".format(
                time() - start, step, survey.lower(), exp_surveys
            )
        )
    if program.lower() not in exp_programs:
        log.warning(
            "{:.1f}s\t{}\tunexpected program={} ({}; proceeding anyway)".format(
                time() - start, step, program.lower(), exp_programs
            )
        )

    # AR folder architecture is now the same at NERSC/KPNO (https://github.com/desihub/fiberassign/issues/302)
    mydirs = {}
    mydirs["sky"] = os.path.join(
        os.getenv("DESI_TARGET"), "catalogs", dr, dtver, "skies"
    )
    mydirs["skysupp"] = os.path.join(
        os.getenv("DESI_TARGET"), "catalogs", gaiadr, dtver, "skies-supp"
    )
    mydirs["gfa"] = os.path.join(
        os.getenv("DESI_TARGET"), "catalogs", dr, dtver, "gfas"
    )
    if program.lower() == "backup":
        dtcat = gaiadr
    else:
        dtcat = dr
    mydirs["targ"] = os.path.join(
        os.getenv("DESI_TARGET"),
        "catalogs",
        dtcat,
        dtver,
        "targets",
        survey.lower(),
        "resolve",
        program.lower(),
    )
    mydirs["mtl"] = os.path.join(
        os.getenv("DESI_SURVEYOPS"), "mtl", survey.lower(), program.lower(),
    )
    # AR secondary (dark, bright; no secondary for backup)
    if program.lower() in ["dark", "bright"]:
        mydirs["scnd"] = os.path.join(
            os.getenv("DESI_TARGET"),
            "catalogs",
            dr,
            dtver,
            "targets",
            survey.lower(),
            "secondary",
            program.lower(),
            "{}targets-{}-secondary.fits".format(survey.lower(), program.lower()),
        )
        mydirs["scndmtl"] = os.path.join(
            os.getenv("DESI_SURVEYOPS"),
            "mtl",
            survey.lower(),
            "secondary",
            program.lower(),
        )

    # AR ToO (same for dark, bright)
    mydirs["too"] = os.path.join(
        os.getenv("DESI_SURVEYOPS"), "mtl", survey.lower(), "ToO", "ToO.ecsv",
    )

    # AR log
    for key in list(mydirs.keys()):
        log.info(
            "{:.1f}s\t{}\tdirectory for {}: {}".format(
                time() - start, step, key, mydirs[key]
            )
        )
        if not os.path.exists(mydirs[key]):
            log.warning(
                "{:.1f}s\t{}\tdirectory for {}: {} does not exist".format(
                    time() - start, step, key, mydirs[key]
                )
            )

    return mydirs


def create_tile(
    tileid,
    tilera,
    tiledec,
    outfn,
    survey,
    obscon="DARK|GRAY|BRIGHT|BACKUP",
    log=None,
    step="",
    start=None,
):
    """
    Create a tiles fits file.
    
    Args:
        tileid: TILEID (int)
        tilera: tile center R.A. (float)
        tiledec: tile center Dec. (float)
        outfn: fits file name to be written
        survey: survey (string; e.g. "sv1", "sv2", "sv3", "main")
        obscon (optional, defaults to "DARK|GRAY|BRIGHT|BACKUP"): tile allowed observing conditions (string)
        log (optional): Logger object
        step (optional): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start (optional): start time for log (in seconds; output of time.time()
    """
    hdr = fitsio.FITSHDR()
    log.info("")
    log.info("")
    log.info("{:.1f}s\t{}\tTIMESTAMP={}".format(time() - start, step, Time.now().isot))
    log.info("{:.1f}s\t{}\tstart generating {}".format(time() - start, step, outfn))
    log.info(
        "{:.1f}s\t{}\ttileid={}, tilera={}, tiledec={}, survey={}, obscon={}".format(
            time() - start, step, tileid, tilera, tiledec, survey, obscon
        )
    )
    d = np.zeros(
        1,
        dtype=[
            ("TILEID", "i4"),
            ("RA", "f8"),
            ("DEC", "f8"),
            ("OBSCONDITIONS", "i4"),
            ("IN_DESI", "i2"),
            ("PROGRAM", "S6"),
        ],
    )
    d["TILEID"] = tileid
    d["RA"] = tilera
    d["DEC"] = tiledec
    d["IN_DESI"] = 1  # AR forcing 1;
    # AR otherwise the default onlydesi=True option in
    # AR desimodel.io.load_tiles() discards tiles outside the desi footprint,
    # AR so return no tiles for the dithered tiles outside desi
    d["PROGRAM"] = survey.upper()  # AR custom... SV2, SV3, MAIN
    log.info("{:.1f}s\t{}\ttile obscon={}".format(time() - start, step, obscon))
    d["OBSCONDITIONS"] = obsconditions.mask(obscon)
    fitsio.write(
        outfn, d, extname="TILES", header=hdr, clobber=True,
    )
    log.info("{:.1f}s\t{}\t{} written".format(time() - start, step, outfn,))


def create_sky(
    tilesfn,
    skydir,
    outfn,
    suppskydir=None,
    tmpoutdir=tempfile.mkdtemp(),
    log=None,
    step="",
    start=None,
):
    """
    Create a sky fits file.
    
    Args:
        tilesfn: path to a tiles fits file (string)
        skydir: desitarget sky folder (string)
        outfn: fits file name to be written (string)
        suppskydir (optional, defaults to None): desitarget suppsky folder (string)
        tmpoutdir (optional, defaults to a temporary directory): temporary directory where
                write_skies will write (creating some sub-directories)
        log (optional): Logger object
        step (optional): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start (optional): start time for log (in seconds; output of time.time()
    """
    log.info("")
    log.info("")
    log.info("{:.1f}s\t{}\tTIMESTAMP={}".format(time() - start, step, Time.now().isot))
    log.info("{:.1f}s\t{}\tstart generating {}".format(time() - start, step, outfn))
    # AR sky: read targets
    tiles = fits.open(tilesfn)[1].data
    skydirs = [skydir]
    if suppskydir is not None:
        skydirs.append(suppskydir)
    d = custom_read_targets_in_tiles(
        skydirs, tiles, quick=True, mtl=False, log=log, step=step, start=start
    )
    n, tmpfn = write_skies(tmpoutdir, d, indir=skydir, indir2=suppskydir,)
    _ = mv_write_targets_out(tmpfn, tmpoutdir, outfn, log=log, step=step, start=start)
    log.info("{:.1f}s\t{}\t{} written".format(time() - start, step, outfn))


def create_gfa(
    tilesfn,
    gfadir,
    survey,
    gaiadr,
    pmcorr,
    outfn,
    tmpoutdir=tempfile.mkdtemp(),
    pmtime_utc_str=None,
    log=None,
    step="",
    start=None,
):
    """
    Create a GFA fits file.
    
    Args:
        tilesfn: path to a tiles fits file (string)
        gfadir: desitarget GFA folder (string)
        survey: survey (string; e.g. "sv1", "sv2", "sv3", "main")
        gaiadr: Gaia dr ("dr2" or "edr3")
        pmcorr: apply proper-motion correction? ("y" or "n")
        outfn: fits file name to be written (string)
        tmpoutdir (optional, defaults to a temporary directory): temporary directory where
                write_skies will write (creating some sub-directories)
        pmtime_utc_str (optional, defaults to None): UTC time use to compute
                new coordinates after applying proper motion since REF_EPOCH
                (string formatted as "yyyy-mm-ddThh:mm:ss+00:00")
        log (optional): Logger object
        step (optional): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start (optional): start time for log (in seconds; output of time.time()

    Notes:
        if pmcorr="y", then pmtime_utc_str needs to be set; will trigger an error otherwise.
    """
    log.info("")
    log.info("")
    log.info("{:.1f}s\t{}\tTIMESTAMP={}".format(time() - start, step, Time.now().isot))
    log.info("{:.1f}s\t{}\tstart generating {}".format(time() - start, step, outfn))
    # AR gfa: read targets
    tiles = fits.open(tilesfn)[1].data
    d = custom_read_targets_in_tiles(
        [gfadir], tiles, quick=True, mtl=False, log=log, step=step, start=start
    )
    log.info(
        "{:.1f}s\t{}\tkeeping {} targets to {}".format(
            time() - start, step, len(d), outfn
        )
    )
    # AR gfa: PMRA, PMDEC: convert NaN to zeros
    d = force_finite_pm(d, log=log, step=step, start=start)
    # AR gfa: update RA, DEC, REF_EPOCH using proper motion?
    if pmcorr == "y":
        if pmtime_utc_str is None:
            log.error(
                "{:.1f}s\t{}\tneed to provide pmtime_utc_str, as proper-correction is requested; exiting".format(
                    time() - start, step,
                )
            )
            sys.exti(1)
        d = update_nowradec(d, gaiadr, pmtime_utc_str, log=log, step=step, start=start)
    else:
        log.info(
            "{:.1f}s\t{}\t*not* applying proper-motion correction".format(
                time() - start, step
            )
        )
        # AR Replaces 0 by force_ref_epoch in ref_epoch
        d = force_nonzero_refepoch(
            d, gaia_ref_epochs[gaiadr], log=log, step=step, start=start
        )

    # AR gfa: write fits
    n, tmpfn = write_targets(tmpoutdir, d, indir=gfadir, survey=survey)
    _ = mv_write_targets_out(tmpfn, tmpoutdir, outfn, log=log, step=step, start=start)
    # AR gfa: update header
    fd = fitsio.FITS(outfn, "rw")
    fd["TARGETS"].write_key("COMMENT", "RA,DEC updated with PM for AEN objects")
    fd["TARGETS"].write_key("COMMENT", "REF_EPOCH updated for all objects")
    fd.close()
    log.info("{:.1f}s\t{}\t{} written".format(time() - start, step, outfn))
