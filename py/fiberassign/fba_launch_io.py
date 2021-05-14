"""
fiberassign.fba_launch_io
=============
Utility functions for fba_launch
"""
from __future__ import absolute_import, division

# system
import os
import subprocess
import sys
import tempfile
import shutil
import re

# time
from time import time
from datetime import datetime, timedelta

#
import numpy as np
import fitsio

# astropy
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy import units
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time

# desitarget
import desitarget
from desitarget.gaiamatch import gaia_psflike
from desitarget.io import read_targets_in_tiles, write_targets, write_skies
from desitarget.mtl import inflate_ledger
from desitarget.targetmask import desi_mask, obsconditions
from desitarget.targets import set_obsconditions
from desitarget.geomask import match

# desimodel
import desimodel
from desimodel.footprint import is_point_in_desi

# fiberassign
import fiberassign
from fiberassign.scripts.assign import parse_assign, run_assign_full
from fiberassign.assign import merge_results, minimal_target_columns
from fiberassign.utils import Logger

# matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib
import matplotlib.image as mpimg

# AR default REF_EPOCH for PMRA=PMDEC=REF_EPOCH=0 objects
gaia_ref_epochs = {"dr2": 2015.5}

# AR tile radius in degrees
tile_radius_deg = 1.628
# AR approx. tile area in degrees
tile_area = np.pi * tile_radius_deg ** 2


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


def get_svn_version(svn_dir):
    """
    Gets the SVN revision number of an SVN folder.

    Args:
        svn_dir: SVN folder path (string)

    Returns:
        svnver: SVN revision number of svn_dir, or "unknown" if not an svn checkout

    Notes:
        `svn_dir` can contain environment variables to expand, e.g. "$DESIMODEL/data"
    """
    cmd = ["svn", "info", "--show-item", "revision", os.path.expandvars(svn_dir)]
    try:
        svn_ver = (
            subprocess.check_output(cmd, stderr=subprocess.DEVNULL).strip().decode()
        )
    except subprocess.CalledProcessError:
        svn_ver = "unknown"

    return svn_ver


def get_program_latest_timestamp(
    program, log=Logger.get(), step="", start=time(),
):
    """
    Get the latest timestamp for a given program from the MTL per-tile file.

    Args:
        program: ideally "dark", "bright", or "backup" (string)
                though if different will return None
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()

    Returns:
        if some entries for input program: UTC YYYY-MM-DDThh:mm:ss+00:00 formatted timestamp (string)
        else: None

    Notes:
        if the per-tile MTL files does not exist or has zero entries,
        TBD: currently add +1min because of a mismatch between the ledgers and the per-tile file.
        TBD: still see if a +1min or +1s is desirable
    """
    # AR check DESI_SURVEYOPS is defined
    assert_env_vars(
        required_env_vars=["DESI_SURVEYOPS"], log=log, step=step, start=start,
    )

    # AR defaults to None (returned if no file or no selected rows)
    timestamp = None

    # AR check if the per-tile file is here
    # AR no need to check the scnd-mtl-done-tiles.ecsv file,
    # AR     as we restrict to a given program ([desi-survey 2434])
    fn = os.path.join(os.getenv("DESI_SURVEYOPS"), "mtl", "mtl-done-tiles.ecsv")
    if os.path.isfile(fn):
        d = Table.read(fn)
        keep = d["PROGRAM"] == program.upper()
        if keep.sum() > 0:
            d = d[keep]
            # AR taking the latest timestamp
            tm = np.unique(d["TIMESTAMP"])[-1]
            # AR does not end with +NN:MM timezone?
            if re.search('\+\d{2}:\d{2}$', tm) is None:
                tm = "{}+00:00".format(tm)
            tm = datetime.strptime(tm, "%Y-%m-%dT%H:%M:%S%z")
            # AR TBD: we currently add one minute; can be removed once
            # AR TBD  update is done on the desitarget side
            tm += timedelta(minutes=1)
            timestamp = tm.isoformat(timespec="seconds")
    return timestamp


def custom_read_targets_in_tiles(
    targdirs,
    tiles,
    quick=True,
    mtl=False,
    unique=True,
    isodate=None,
    log=Logger.get(),
    step="",
    start=time(),
):
    """
    Wrapper to desitarget.io.read_targets_in_tiles, allowing multiple folders
    and sanity check TARGETID if more than one folder.

    Args:
        targdirs: list of folders
        tiles: tiles object (as required by desitarget.io.read_targets_in_tiles)
        quick, mtl, unique, isodate (optional): same as desitarget.io.read_targets_in_tiles arguments
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()
    Returns:
        array of targets in the passed tiles.
    """
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


def mv_write_targets_out(infn, targdir, outfn, log=Logger.get(), step="", start=time()):
    """
    Moves the file created by desitarget.io.write_targets
    and removes folder created by desitarget.io.write_targets    
    
    Args:
        infn: filename output by desitarget.io.write_targets
        targdir: folder provided as desitarget.io.write_targets input
        outfn: desired renaming of infn
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()
    """
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
    d, pmra_key="PMRA", pmdec_key="PMDEC", log=Logger.get(), step="", start=time()
):
    """
    Replaces NaN PMRA, PMDEC by 0    
    
    Args:
        d: array with at least proper-motion columns
        pmra_key (optional, defaults to PMRA): column name for PMRA
        pmdec_key (optional, defaults to PMDEC): column name for PMDEC
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()

    Returns:
        d: same as input d, but NaN proper motions replaced by 0
        
    """
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
    log=Logger.get(),
    step="",
    start=time(),
):
    """
    Replaces 0 by force_ref_epoch in ref_epoch
    
    Args:
        d: array with at least proper-motion columns
        force_ref_epoch: float, ref_epoch to replace 0 by
        ref_epoch_key (optional, defaults to REF_EPOCH): column name for the ref_epoch
        pmra_key (optional, defaults to PMRA): column name for PMRA
        pmdec_key (optional, defaults to PMDEC): column name for PMDEC
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()        
    Returns:
        d: same as input d, but 0 ref_epochs replaced by force_ref_epoch

    Notes:
        Will exit with error if ref_epoch=0, but pmra or pmdec != 0
        
    """
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
    log=Logger.get(),
    step="",
    start=time(),
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
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()        
        
    Returns:
        d: same as input, but with RA, DEC updated to pmtime_utc_str

    Notes:
        REF_EPOCH is updated for *all* objects
    """
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
    log=Logger.get(),
    step="settings",
    start=time(),
):
    """
    Assert the environment variables required by fba_launch 
    
    Args:
        required_env_vars (optional, defaults to ["DESI_ROOT",
        "DESI_TARGET",
        "DESIMODEL",
        "DESI_SURVEYOPS",
        "SKYBRICKS_DIR",]): list of environment variables required by fba_launch
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()        
        
    Notes:
        will exit with error if some assertions are not verified 
    """
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
    log=Logger.get(),
    step="settings",
    start=time(),
):
    """
    Assert the fba_launch date arguments are correctly formatted ("YYYY-MM-DDThh:mm:ss+00:00")
    
    Args:
        args: fba_launch parser.parse_args() output
        dates (optional, defaults to ["pmtime_utc_str", "rundate", "mtltime"]): list of date fba_launch argument names to check
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()        
        
    Notes:
        will exit with error if some assertions are not verified 
    """
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
    tileid, forcetileid="n", log=Logger.get(), step="settings", start=time(),
):
    """
    Asserts if TILEID already exists in the SVN tile folder
    
    Args:
        tileid: TILEID to check (int)
        forcetileid (optional, defaults to "n"): "y" or "n";
                if "n", will trigger a warning + an error
                if "y", e will trigger a warning only
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()        
       
    """
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
    log=Logger.get(),
    step="settings",
    start=time(),
):
    """
    Print various configuration informations (machine, modules version/path, DESI environment variables).
    
    Args:
        required_env_vars (optional, defaults to ["DESI_ROOT",
        "DESI_TARGET",
        "DESIMODEL",
        "DESI_SURVEYOPS",
        "SKYBRICKS_DIR",]): list of environment variables required by fba_launch
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()        
    """
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
    log=Logger.get(),
    step="settings",
    start=time(),
):
    """
    Obtain the folder/file full paths for desitarget products
    
    Args:
        dtver: desitarget catalog version (string; e.g., "0.57.0")
        survey: survey (string; e.g. "sv1", "sv2", "sv3", "main")
        program: "dark", "bright", or "backup" (string)
        dr (optional, defaults to "dr9"): legacypipe dr (string)
        gaiadr (optional, defaults to "gaiadr2"): gaia dr (string)
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()        
        
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
        if survey.lower() == "main":
            basename = "targets-{}-secondary.fits".format(program.lower())
        else:
            basename = "{}targets-{}-secondary.fits".format(survey.lower(), program.lower())
        mydirs["scnd"] = os.path.join(
            os.getenv("DESI_TARGET"),
            "catalogs",
            dr,
            dtver,
            "targets",
            survey.lower(),
            "secondary",
            program.lower(),
            basename,
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
    log=Logger.get(),
    step="",
    start=time(),
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
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()
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
    log=Logger.get(),
    step="",
    start=time(),
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
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()
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


def create_targ_nomtl(
    tilesfn,
    targdir,
    survey,
    gaiadr,
    pmcorr,
    outfn,
    tmpoutdir=tempfile.mkdtemp(),
    pmtime_utc_str=None,
    quick=True,
    log=Logger.get(),
    step="",
    start=time(),
):
    """
    Create a target fits file, with solely using desitarget catalogs, no MTL.
        e.g. for the GFA, but could be used for other purposes.
    
    Args:
        tilesfn: path to a tiles fits file (string)
        targdir: desitarget target folder (string)
        survey: survey (string; e.g. "sv1", "sv2", "sv3", "main")
        gaiadr: Gaia dr ("dr2" or "edr3")
        pmcorr: apply proper-motion correction? ("y" or "n")
        outfn: fits file name to be written (string)
        tmpoutdir (optional, defaults to a temporary directory): temporary directory where
                write_skies will write (creating some sub-directories)
        pmtime_utc_str (optional, defaults to None): UTC time use to compute
                new coordinates after applying proper motion since REF_EPOCH
                (string formatted as "yyyy-mm-ddThh:mm:ss+00:00")
        quick (optional, defaults to True): boolean, arguments of desitarget.io.read_targets_in_tiles()
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()

    Notes:
        if pmcorr="y", then pmtime_utc_str needs to be set; will trigger an error otherwise.
    """
    log.info("")
    log.info("")
    log.info("{:.1f}s\t{}\tTIMESTAMP={}".format(time() - start, step, Time.now().isot))
    log.info("{:.1f}s\t{}\tstart generating {}".format(time() - start, step, outfn))
    # AR targ_nomtl: read targets
    tiles = fits.open(tilesfn)[1].data
    d = custom_read_targets_in_tiles(
        [targdir], tiles, quick=quick, mtl=False, log=log, step=step, start=start
    )
    log.info(
        "{:.1f}s\t{}\tkeeping {} targets to {}".format(
            time() - start, step, len(d), outfn
        )
    )
    # AR targ_nomtl: PMRA, PMDEC: convert NaN to zeros
    d = force_finite_pm(d, log=log, step=step, start=start)
    # AR targ_nomtl: update RA, DEC, REF_EPOCH using proper motion?
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
        # AR targ_nomtl: Replaces 0 by force_ref_epoch in ref_epoch
        d = force_nonzero_refepoch(
            d, gaia_ref_epochs[gaiadr], log=log, step=step, start=start
        )

    # AR targ_nomtl: write fits
    n, tmpfn = write_targets(tmpoutdir, d, indir=targdir, survey=survey)
    _ = mv_write_targets_out(tmpfn, tmpoutdir, outfn, log=log, step=step, start=start)
    # AR targ_nomtl: update header if pmcorr = "y"
    if pmcorr == "y":
        fd = fitsio.FITS(outfn, "rw")
        fd["TARGETS"].write_key("COMMENT", "RA,DEC updated with PM for AEN objects")
        fd["TARGETS"].write_key("COMMENT", "REF_EPOCH updated for all objects")
        fd.close()
    log.info("{:.1f}s\t{}\t{} written".format(time() - start, step, outfn))


def create_mtl(
    tilesfn,
    mtldir,
    mtltime,
    targdir,
    survey,
    gaiadr,
    pmcorr,
    outfn,
    tmpoutdir=tempfile.mkdtemp(),
    pmtime_utc_str=None,
    log=Logger.get(),
    step="",
    start=time(),
):
    """
    Create a (primary or secondary) target fits file, based on MTL ledgers (and complementary columns from desitarget targets files).
    
    Args:
        tilesfn: path to a tiles fits file (string)
        mtldir: desisurveyops MTL folder (string)
        mtltime: MTL isodate (string formatted as yyyy-mm-ddThh:mm:ss+00:00)
        targdir: desitarget targets folder (or file name if secondary) (string)
        survey: survey (string; e.g. "sv1", "sv2", "sv3", "main")
        gaiadr: Gaia dr ("dr2" or "edr3")
        pmcorr: apply proper-motion correction? ("y" or "n")
        outfn: fits file name to be written (string)
        tmpoutdir (optional, defaults to a temporary directory): temporary directory where
                write_skies will write (creating some sub-directories)
        pmtime_utc_str (optional, defaults to None): UTC time use to compute
                new coordinates after applying proper motion since REF_EPOCH
                (string formatted as "yyyy-mm-ddThh:mm:ss+00:00")
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()

    Notes:
        if pmcorr="y", then pmtime_utc_str needs to be set; will trigger an error otherwise.
        for sv3-backup, we remove BACKUP_BRIGHT targets.
        TBD : if secondary targets, we currently disable the inflate_ledger(), as it
                seems to not currently work.
                hence if secondary and pmcorr="y", the code will crash, as the 
                GAIA_ASTROMETRIC_EXCESS_NOISE column will be missing; though we do not
                expect this configuration to happen, so it should be fine for now.
    """
    log.info("")
    log.info("")
    log.info("{:.1f}s\t{}\tTIMESTAMP={}".format(time() - start, step, Time.now().isot))
    log.info("{:.1f}s\t{}\tstart generating {}".format(time() - start, step, outfn))
    tiles = fits.open(tilesfn)[1].data
    # AR mtl: storing the timestamp at which we queried MTL
    log.info("{:.1f}s\t{}\tmtltime={}".format(time() - start, step, mtltime))

    # AR mtl: read mtl
    d = custom_read_targets_in_tiles(
        [mtldir],
        tiles,
        quick=False,
        mtl=True,
        unique=True,
        isodate=mtltime,
        log=log,
        step=step,
        start=start,
    )

    # AR mtl: removing by hand BACKUP_BRIGHT for sv3/BACKUP
    # AR mtl: using an indirect way to find if program=backup,
    # AR mtl:   to avoid the need of an extra program argument
    # AR mtl:   for sv3, there is no secondary-backup, so no ambiguity
    if (survey == "sv3") & ("backup" in mtldir):
        from desitarget.sv3.sv3_targetmask import mws_mask

        keep = (d["SV3_MWS_TARGET"] & mws_mask["BACKUP_BRIGHT"]) == 0
        log.info(
            "{:.1f}s\t{}\tremoving {}/{} BACKUP_BRIGHT targets".format(
                time() - start, step, len(d) - keep.sum(), len(d)
            )
        )
        d = d[keep]

    # AR mtl: add columns not present in ledgers
    # AR mtl: need to provide exact list (if columns=None, inflate_ledger()
    # AR mtl:    overwrites existing columns)
    # AR mtl: TBD : we currently disable it for secondary targets
    # AR mtl:       using an indirect way to find if secondary,
    # AR mtl:       to avoid the need of an extra program argument
    if "secondary" not in mtldir:
        columns = [key for key in minimal_target_columns if key not in d.dtype.names]
        # AR mtl: also add GAIA_ASTROMETRIC_EXCESS_NOISE, in case args.pmcorr == "y"
        if pmcorr == "y":
            columns += ["GAIA_ASTROMETRIC_EXCESS_NOISE"]
        log.info(
            "{:.1f}s\t{}\tadding {} from {}".format(
                time() - start, step, ",".join(columns), targdir
            )
        )
        d = inflate_ledger(
            d, targdir, columns=columns, header=False, strictcols=False, quick=True
        )

    # AR mtl: PMRA, PMDEC: convert NaN to zeros
    d = force_finite_pm(d, log=log, step=step, start=start)
    # AR mtl: update RA, DEC, REF_EPOCH using proper motion?
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
    # AR mtl: write fits
    n, tmpfn = write_targets(tmpoutdir, d, indir=mtldir, indir2=targdir, survey=survey)
    _ = mv_write_targets_out(tmpfn, tmpoutdir, outfn, log=log, step=step, start=start,)

    # AR mtl: update header if pmcorr = "y"
    if pmcorr == "y":
        fd = fitsio.FITS(outfn, "rw")
        fd["TARGETS"].write_key("COMMENT", "RA,DEC updated with PM for AEN objects")
        fd["TARGETS"].write_key("COMMENT", "REF_EPOCH updated for all objects")
        fd.close()
    log.info("{:.1f}s\t{}\t{} written".format(time() - start, step, outfn))


def create_too(
    tilesfn,
    toofn,
    mjd_min,
    mjd_max,
    survey,
    gaiadr,
    pmcorr,
    outfn,
    tmpoutdir=tempfile.mkdtemp(),
    pmtime_utc_str=None,
    too_tile=False,
    log=Logger.get(),
    step="",
    start=time(),
):
    """
    Create a ToO target fits file, with selecting targets in a MJD time window.
    
    Args:
        tilesfn: path to a tiles fits file (string)
        toofn: ToO file name (string)
        mjd_min, mjd_max (floats): we keep targets with MJD_BEGIN < mjd_max and MJD_END > mjd_min 
        survey: survey (string; e.g. "sv1", "sv2", "sv3", "main")
        gaiadr: Gaia dr ("dr2" or "edr3")
        pmcorr: apply proper-motion correction? ("y" or "n")
        outfn: fits file name to be written (string)
        tmpoutdir (optional, defaults to a temporary directory): temporary directory where
                write_skies will write (creating some sub-directories)
        pmtime_utc_str (optional, defaults to None): UTC time use to compute
                new coordinates after applying proper motion since REF_EPOCH
                (string formatted as "yyyy-mm-ddThh:mm:ss+00:00")
        too_tile (optional, defaults to False): if False, we only keep TOO_TYPE!="TILE",
                if True, we do not cut on TOO_TYPE, hence keeping both TOO_TYPE="FIBER" *and*
                TOO_TYPE="TILE" for ToO dedicated tiles (boolean)
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()

    Notes:
        if pmcorr="y", then pmtime_utc_str needs to be set; will trigger an error otherwise.
        TBD : the MJD window to accept targets; currently in fba_launch, we set a month
                from the tile design date;
                it surely needs to be updated/refined once operations are more clear.
        some steps in common with create_mtl().
    """
    log.info("")
    log.info("")
    log.info("{:.1f}s\t{}\tTIMESTAMP={}".format(time() - start, step, Time.now().isot))
    log.info("{:.1f}s\t{}\tstart generating {}".format(time() - start, step, outfn))
    tiles = fits.open(tilesfn)[1].data

    # AR too: read too file
    # AR cut on:
    # AR - tiles
    # AR - mjd (! TBD !)
    d = Table.read(toofn)
    keep = is_point_in_desi(tiles, d["RA"], d["DEC"])
    if not too_tile:
        keep &= d["TOO_TYPE"] != "TILE"
    keep &= (d["MJD_BEGIN"] < mjd_max) & (d["MJD_END"] > mjd_min)
    log.info(
        "{:.1f}s\t{}\tkeeping {}/{} targets in tiles, with TOO_TYPE={}, and in the MJD time window: {}, {}".format(
            time() - start, step, keep.sum(), len(keep), "TILE,FIBER" if too_tile else "FIBER", mjd_min, mjd_max
        )
    )

    if keep.sum() > 0:
        d = d[keep]
        # AR too: PMRA, PMDEC: convert NaN to zeros
        d = force_finite_pm(d, log=log, step=step, start=start)

        # AR too: update RA, DEC, REF_EPOCH using proper motion
        if pmcorr == "y":
            if pmtime_utc_str is None:
                log.error(
                    "{:.1f}s\t{}\tneed to provide pmtime_utc_str, as proper-correction is requested; exiting".format(
                        time() - start, step,
                    )
                )
                sys.exti(1)
            d = update_nowradec(
                d, gaiadr, pmtime_utc_str, log=log, step=step, start=start
            )
        else:
            log.info(
                "{:.1f}s\t{}\t*not* applying proper-motion correction".format(
                    time() - start, step
                )
            )
            # AR single REF_EPOCH needed
            # AR TBD currently all targets have PMRA=PMDEC=0,
            # AR TBD so it s fine to just change all REF_EPOCH
            d["REF_EPOCH"] = np.zeros(len(d))
            # AR Replaces 0 by force_ref_epoch in ref_epoch
            d = force_nonzero_refepoch(
                d, gaia_ref_epochs[gaiadr], log=log, step=step, start=start
            )

        # AR mtl: write fits
        n, tmpfn = write_targets(tmpoutdir, d.as_array(), indir=toofn, survey=survey)
        _ = mv_write_targets_out(
            tmpfn, tmpoutdir, outfn, log=log, step=step, start=start,
        )
        # AR mtl: update header if pmcorr = "y"
        if pmcorr == "y":
            fd = fitsio.FITS(outfn, "rw")
            fd["TARGETS"].write_key("COMMENT", "RA,DEC updated with PM for AEN objects")
            fd["TARGETS"].write_key("COMMENT", "REF_EPOCH updated for all objects")
            fd.close()
        log.info("{:.1f}s\t{}\t{} written".format(time() - start, step, outfn))

    else:
        log.info(
            "{:.1f}s\t{}\tno too kept too targets, no {} written".format(
                time() - start, step, outfn
            )
        )


def launch_onetile_fa(
    args,
    tilesfn,
    targfns,
    fbafn,
    fiberassignfn,
    skyfn=None,
    gfafn=None,
    log=Logger.get(),
    step="",
    start=time(),
):
    """
    Runs the fiber assignment (run_assign_full) and merges the results (merge_results) for a single tile.
    
    Args:
        args: fba_launch-like parser.parse_args() output
            should contain at least:
                - rundate
                - sky_per_petal
                - standards_per_petal
                - sky_per_slitblock
        tilesfn: path to the input tiles fits file (string)
        targfns: paths to the input targets fits files, e.g. targ, scnd, too (either a string if only one file, or a list of strings)
        fbafn: path to the output fba-TILEID.fits file (string)
        fiberassignfn: path to the output fiberassign-TILEID.fits file (string)
        skyfn (optional, defaults to None): path to a sky fits file (string)
        gfafn (optional, defaults to None): path to a gfa fits file (string)
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()

    Notes:
        no sanity checks done on inputs; assumed to be done elsewhere
        assumes the output directory is the same for fbafn and fiberassignfn
        we keep a generic "args" input, so that any later added argument in fba_launch does not
            requires a change in the launch_fa() call format.
        fba_launch-like adding information in the header is done in another function, update_fiberassign_header
        TBD: be careful if working in the SVN-directory; maybe add additional safety lines? 
    """
    log.info("")
    log.info("")
    log.info("{:.1f}s\t{}\tTIMESTAMP={}".format(time() - start, step, Time.now().isot))
    log.info("{:.1f}s\t{}\tstart running fiber assignment".format(time() - start, step))

    # AR convert targfns to list if string (i.e. only one input file)
    if isinstance(targfns, str):
        targfns = [targfns]

    # AR tileid
    tileid = fits.open(tilesfn)[1].data["TILEID"][0]

    # AR output directory (picking the one of fbafn)
    outdir = os.path.dirname(fbafn)

    # AR safe: delete possibly existing fba-{tileid}.fits and fiberassign-{tileid}.fits
    # AR TBD: add additional safety check if running in SVN folder?
    if os.path.isfile(fbafn):
        os.remove(fbafn)
    if os.path.isfile(fiberassignfn):
        os.remove(fiberassignfn)

    # AR preparing fba_run inputs
    opts = [
        "--targets",
    ]
    for targfn in targfns:
        opts += [
            targfn,
        ]
    opts += [
        "--overwrite",
        "--write_all_targets",
        "--dir",
        outdir,
        "--footprint",
        tilesfn,
        "--rundate",
        args.rundate,
        "--sky_per_petal",
        args.sky_per_petal,
        "--standards_per_petal",
        args.standards_per_petal,
        "--sky_per_slitblock",
        str(args.sky_per_slitblock),
        "--ha",
        str(args.ha),
    ]
    if skyfn is not None:
        opts += [
            "--sky",
            skyfn,
        ]
    if gfafn is not None:
        opts += [
            "--gfafile",
            gfafn,
        ]
    log.info(
        "{:.1f}s\t{}\ttileid={:06d}: running raw fiber assignment (run_assign_full) with opts={}".format(
            time() - start, step, tileid, " ; ".join(opts)
        )
    )
    ag = parse_assign(opts)
    run_assign_full(ag)

    # AR merging
    # AR not using run_merge(), because it looks for all fba-TILEID.fits file
    # AR in the out directory...
    ag = {}
    ag["tiles"] = [tileid]
    ag["columns"] = None
    ag["targets"] = [gfafn] + targfns
    if skyfn is not None:
        ag["sky"] = [skyfn]
    else:
        ag["sky"] = []
    ag["result_dir"] = outdir
    ag["copy_fba"] = False
    tmparr = []
    for key in list(ag.keys()):
        tmparr += ["{} = {}".format(key, ag[key])]
    log.info(
        "{:.1f}s\t{}\ttileid={:06d}: merging input target data (merge_results) with argument={}".format(
            time() - start, step, tileid, " ; ".join(tmparr)
        )
    )
    merge_results(
        ag["targets"],
        ag["sky"],
        ag["tiles"],
        result_dir=ag["result_dir"],
        columns=ag["columns"],
        copy_fba=ag["copy_fba"],
    )


def update_fiberassign_header(
    fiberassignfn,
    args,
    mydirs,
    hdr_survey,
    hdr_faprgrm,
    faflavor,
    ebv,
    obscon,
    fascript,
    log=Logger.get(),
    step="",
    start=time(),
):
    """
    Adds various information in the fiberassign-TILEID.fits PRIMARY header.
    
    Args:
        fiberassignfn: path to fiberassign-TILEID.fits file (string)
        args: fba_launch-like parser.parse_args() output
            should contain at least (see fba_launch arguments):
                - outdir, survey, program, rundate, pmcorr, pmtime_utc_str
                - faprgrm, mtltime, goaltime, goaltype, sbprof, mintfrac
            will also be used to store in FAARGS the list of input arguments of the fba_launch call.
        mydirs: dictionary with the desitarget paths; ideally:
                - sky: sky folder
                - skysupp: skysupp folder
                - gfa: GFA folder
                - targ: targets folder (static catalogs, with all columns)
                - mtl: MTL folder
                - scnd: secondary fits catalog (static)
                - scndmtl: MTL folder for secondary targets
                - too: ToO ecsv catalog
        hdr_survey: value for the SURVEY keyword (string)
        hdr_faprgrm: value for the FAPRGRM keyword (string)
        faflavor: usually {survey}{program} in lower cases (string)
        ebv: median EBV over the tile targets (float)
        obscon: tile allowed observing conditions (string; e.g. "DARK|GRAY|BRIGHT|BACKUP") 
        fascript: fba_launch-like script used to designed the tile; in case of different scripts for dedicated tiles
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()

    Notes:
        no check is done on mydirs.
        hdr_survey, hdr_faprgrm: for the "regular" surveys (e.g., sv3, main), those will be args.survey, args.program
            but for dedicated survey, they will (have to) be different. 
        faflavor has to be {hdr_survey}{hdr_faprgrm}; will exit with an error if not;
            keeping this to be sure it is not forgotten to be done for dedicated programs.
    """
    # AR sanity check on faflavor
    if faflavor != "{}{}".format(hdr_survey, hdr_faprgrm):
        log.error(
            "{:.1f}s\t{}\tfaflavor={} inconsitent with hdr_survey={} and hdr_faprgrm={}; exiting".format(
                time() - start, step, faflavor, hdr_survey, hdr_faprgrm,
            )
        )
        sys.exit(1)

    # AR propagating some settings into the PRIMARY header
    fd = fitsio.FITS(fiberassignfn, "rw")

    # AR faflavor
    fd["PRIMARY"].write_key("FAFLAVOR", faflavor)

    # AR folders, with replacing $DESI_ROOT by DESIROOT
    desiroot = os.getenv("DESI_ROOT")
    fd["PRIMARY"].write_key("DESIROOT", desiroot)
    for key in np.sort(list(mydirs.keys())):
        if (key == "mtl") & (isinstance(mydirs["mtl"], list)):
            # AR header keywords: MTL, MTL2, MTL3, etc
            # AR probably to be deprecate for sv2
            suffixs = [""] + np.arange(2, len(mydirs["mtl"]) + 1).astype(str).tolist()
            for mtldir, suffix in zip(mydirs["mtl"], suffixs):
                fd["PRIMARY"].write_key(
                    "mtl{}".format(suffix), mtldir.replace(desiroot, "DESIROOT"),
                )
        else:
            fd["PRIMARY"].write_key(key, mydirs[key].replace(desiroot, "DESIROOT"))
    # AR storing some specific arguments
    # AR plus a (long) FAARGS keyword with storing arguments to re-run the fiber assignment
    # AR     we exclude from FAARGS outdir, forcetiled, and any None argument
    tmparr = []
    for kwargs in args._get_kwargs():
        if (kwargs[0].lower() not in ["outdir", "forcetileid"]) & (
            kwargs[1] is not None
        ):
            tmparr += ["--{} {}".format(kwargs[0], kwargs[1])]
    fd["PRIMARY"].write_key(
        "faargs", " ".join(tmparr),
    )
    # AR some keywords
    fd["PRIMARY"].write_key("outdir", args.outdir)
    fd["PRIMARY"].write_key("survey", hdr_survey)  # AR not args.survey!
    fd["PRIMARY"].write_key("rundate", args.rundate)
    fd["PRIMARY"].write_key("pmcorr", args.pmcorr)
    fd["PRIMARY"].write_key("pmtime", args.pmtime_utc_str)
    fd["PRIMARY"].write_key("faprgrm", hdr_faprgrm)  # AR not args.program!
    fd["PRIMARY"].write_key("mtltime", args.mtltime)
    fd["PRIMARY"].write_key("obscon", obscon)
    # AR informations for NTS
    # AR SBPROF from https://desi.lbl.gov/trac/wiki/SurveyOps/SurveySpeed#NominalFiberfracValues
    # AR version 35
    fd["PRIMARY"].write_key("goaltime", args.goaltime)
    fd["PRIMARY"].write_key("goaltype", args.program)
    fd["PRIMARY"].write_key("ebvfac", 10.0 ** (2.165 * np.median(ebv) / 2.5))
    fd["PRIMARY"].write_key("sbprof", args.sbprof)
    fd["PRIMARY"].write_key("mintfrac", args.mintfrac)
    # AR fba_launch-like script name used to designed the tile
    fd["PRIMARY"].write_key("fascript", fascript)
    # AR SVN revision number
    fd["PRIMARY"].write_key(
        "svndm", get_svn_version(os.path.join(os.getenv("DESIMODEL"), "data"))
    )
    fd["PRIMARY"].write_key(
        "svnmtl", get_svn_version(os.path.join(os.getenv("DESI_SURVEYOPS"), "mtl"))
    )
    fd.close()


def secure_gzip(
    fiberassignfn, log=Logger.get(), step="", start=time(),
):
    """
    Secure gzipping of the fiberassign-TILEID.fits file.
    
    Args:
        fiberassignfn: path to fiberassign-TILEID.fits file (string)
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()
    """
    log.info("")
    log.info("")
    log.info("{:.1f}s\t{}\tTIMESTAMP={}".format(time() - start, step, Time.now().isot))
    if os.path.isfile("{}.gz".format(fiberassignfn)):
        os.remove("{}.gz".format(fiberassignfn))
        log.info(
            "{:.1f}s\t{}\tdeleting existing {}.gz".format(
                time() - start, step, fiberassignfn
            )
        )
    os.system("gzip {}".format(fiberassignfn))
    log.info("{:.1f}s\t{}\tgzipping {}".format(time() - start, step, fiberassignfn))


def get_dt_masks(
    survey, log=Logger.get(), step="", start=time(),
):
    """
    Get the desitarget masks for a survey.
    
    Args:
        survey: survey name: "sv1", "sv2", "sv3" or "main") (string)
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()

    Returns:
        yaml_masks: dictionary with storing in the
            "DESI_TARGET", "BGS_TARGET", "MWS_TARGET", "SCND_TARGET" keys
            the corresponding desitarget YAML masks for survey
        wd_mskkeys: list of keys identifying the WDs
        wd_msks: list of masks identifying the WDs
        std_mskkeys: list of keys identifying the STDs
        std_msks: list of masks identifying the STDs

    Notes:
        close to desitarget.targets.main_cmx_or_sv,
        but using a dictionary, more adapted to this code
    """
    if survey == "sv1":
        from desitarget.sv1 import sv1_targetmask as targetmask
    elif survey == "sv2":
        from desitarget.sv2 import sv2_targetmask as targetmask
    elif survey == "sv3":
        from desitarget.sv3 import sv3_targetmask as targetmask
    elif survey == "main":
        from desitarget import targetmask
        from fiberassign.targets import default_main_stdmask
    else:
        log.error(
            "{:.1f}s\t{}\tsurvey={} is not in sv1, sv2, sv3 or main; exiting".format(
                time() - start, step, survey
            )
        )
        sys.exit(1)

    # AR YAML masks
    yaml_masks = {
        "DESI_TARGET": targetmask.desi_mask,
        "BGS_TARGET": targetmask.bgs_mask,
        "MWS_TARGET": targetmask.mws_mask,
        "SCND_TARGET": targetmask.scnd_mask,
    }

    # AR WD masks
    wd_mskkeys, wd_msks = [], []
    for mskkey in ["DESI_TARGET", "MWS_TARGET"]:
        wd_mskkeys += [mskkey for key in yaml_masks[mskkey].names() if "_WD" in key]
        wd_msks += [key for key in yaml_masks[mskkey].names() if "_WD" in key]

    # AR STD masks
    std_mskkeys, std_msks = [], []
    for mskkey in ["DESI_TARGET", "MWS_TARGET"]:
        std_mskkeys += [mskkey for key in yaml_masks[mskkey].names() if "STD" in key]
        std_msks += [key for key in yaml_masks[mskkey].names() if "STD" in key]

    return yaml_masks, wd_mskkeys, wd_msks, std_mskkeys, std_msks


def get_qa_tracers(
    survey, program, log=Logger.get(), step="", start=time(),
):
    """
    Returns the tracers for which we provide QA plots of fiber assignment.
    
    Args:
        survey: survey name: "sv1", "sv2", "sv3" or "main") (string)
        program: "DARK", "BRIGHT", or "BACKUP" (string)
        log (optional, defaults to Logger.get()): Logger object                                                                                                                                            
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()
                                                                                                                                                                                 
    Returns:
        trmskkeys: list of keys to select the mask on (list of strings)
        trmsks: list of mask names (list of strings)
    """
    if program == "DARK":
        trmskkeys = ["DESI_TARGET", "DESI_TARGET", "DESI_TARGET"]
        trmsks = ["LRG", "ELG", "QSO"]
    elif program == "BRIGHT":
        trmskkeys = ["BGS_TARGET", "BGS_TARGET"]
        trmsks = ["BGS_BRIGHT", "BGS_FAINT"]
        trmskkeys += ["MWS_TARGET", "MWS_TARGET"]
        if survey == "sv1":
            trmsks += ["MWS_MAIN_BROAD", "MWS_NEARBY"]
        else:
            trmsks += ["MWS_BROAD", "MWS_NEARBY"]
    elif program == "BACKUP":
        trmskkeys = ["MWS_TARGET", "MWS_TARGET", "MWS_TARGET"]
        trmsks = ["BACKUP_BRIGHT", "BACKUP_FAINT", "BACKUP_VERY_FAINT"]
    else:
        log.error(
            "{:.1f}s\t{}\tprogram={} not in DARK, BRIGHT, or BACKUP; exiting".format(
                time() - start, step, program
            )
        )
        sys.exit(1)

    return trmskkeys, trmsks


def get_parent_assign_quants(
    survey,
    targfns,
    fiberassignfn,
    tilera,
    tiledec,
    log=Logger.get(),
    step="",
    start=time(),
):
    """
    Stores the parent and assigned targets properties (desitarget columns).
    
    Args:
        survey: survey (string; e.g. "sv1", "sv2", "sv3", "main")
        targfns: paths to the input targets fits files, e.g. targ, scnd, too (either a string if only one file, or a list of strings)
        fiberassignfn: path to the output fiberassign-TILEID.fits file (string)
        tilera: tile center R.A. (float)
        tiledec: tile center Dec. (float)
        log (optional, defaults to Logger.get()): Logger object                                                                                                                                            
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()

    Returns:
        parent: dictionary of the parent target sample, with each key being some desitarget column
        assign: same as parent, with similar row-ordering (filling with zeros or NaNs if not assigned)
        dras: dictionary with projected distance (degrees) along R.A. to the center of the tile (np.array of floats),
            for each of the following subsample: "parent", "assign", "sky", "bad", "wd", "std" (all assigned subsamples,
            except parent)
        ddecs: same as dras, for projected distances along Dec.
        petals: dictionary with PETAL_LOC (np.array of floats) for each of the assigned "sky", "bad", "wd", "std" subsamples
        nassign: dictionary with the number of assigned fibers for each of the assigned "SKY", "BAD", "TGT", "WD", "STD"  subsamples
    """
    # AR convert targfns to list if string (i.e. only one input file)
    if isinstance(targfns, str):
        targfns = [targfns]

    # AR initializing dictionaires
    parent, assign, dras, ddecs, petals, nassign = {}, {}, {}, {}, {}, {}

    # AR YAML and WD and STD masks
    yaml_masks, wd_mskkeys, wd_msks, std_mskkeys, std_msks = get_dt_masks(
        survey, log=log, step=step, start=start,
    )

    # AR keys we use (plus few for assign)
    keys = [
        "TARGETID",
        "FLUX_G",
        "FLUX_R",
        "FIBERTOTFLUX_R",
        "FLUX_Z",
        "FLUX_W1",
        "FLUX_W2",
        "EBV",
        "GAIA_PHOT_G_MEAN_MAG",
        "RA",
        "DEC",
        "DESI_TARGET",
        "BGS_TARGET",
        "MWS_TARGET",
        "SCND_TARGET",
    ]

    # AR parent
    for key in keys:
        parent[key] = []
    for targfn in targfns:
        d = fits.open(targfn)[1].data
        for key in keys:
            if key in ["DESI_TARGET", "BGS_TARGET", "MWS_TARGET", "SCND_TARGET",]:
                if survey.lower()[:2] == "sv":
                    key_orig = "{}_{}".format(survey.upper(), key)
                else:
                    key_orig = key
                if key_orig in d.dtype.names:
                    parent[key] += d[key_orig].tolist()
                else:
                    parent[key] += [0 for x in d["RA"]]
            # AR flux, ebv for secondary
            elif key not in d.dtype.names:
                parent[key] += [0.0 for x in d["RA"]]
            else:
                parent[key] += d[key].tolist()
    for key in keys:
        parent[key] = np.array(parent[key])
    dras["parent"], ddecs["parent"] = get_tpos(
        tilera, tiledec, parent["RA"], parent["DEC"]
    )

    # AR fiberassign
    d = fits.open(fiberassignfn)[1].data
    # AR
    for key in ["SKY", "BAD", "TGT"]:
        nassign[key] = (d["OBJTYPE"] == key).sum()
    # AR SKY
    keep = d["OBJTYPE"] == "SKY"
    dras["sky"], ddecs["sky"] = get_tpos(
        tilera, tiledec, d["TARGET_RA"][keep], d["TARGET_DEC"][keep]
    )
    petals["sky"] = d["PETAL_LOC"][keep]
    # AR BAD
    keep = d["OBJTYPE"] == "BAD"
    dras["bad"], ddecs["bad"] = get_tpos(
        tilera, tiledec, d["TARGET_RA"][keep], d["TARGET_DEC"][keep]
    )
    petals["bad"] = d["PETAL_LOC"][keep]
    # AR TGT
    # AR arrays twinning the parent ordering, with nans/zeros
    # AR e.g. SV2_DESI_TARGET -> DESI_TARGET
    d = d[d["OBJTYPE"] == "TGT"]
    # AR TARGETIDs are unique in both arrays, so we can use geomask
    iip, ii = match(parent["TARGETID"], d["TARGETID"])
    keys = [key for key in keys if key != "RA" and key != "DEC"]
    keys += [
        "TARGET_RA",
        "TARGET_DEC",
        "PETAL_LOC",
    ]
    for key in keys:
        if key in [
            "TARGETID",
            "DESI_TARGET",
            "BGS_TARGET",
            "MWS_TARGET",
            "SCND_TARGET",
        ]:
            assign[key] = np.zeros(len(parent["TARGETID"]), dtype=int)
            if (key != "TARGETID") & (survey.lower()[:2] == "sv"):
                assign[key][iip] = d["{}_{}".format(survey.upper(), key)][ii]
            else:
                assign[key][iip] = d[key][ii]
        else:
            assign[key] = np.nan + np.zeros(len(parent["TARGETID"]))
            assign[key][iip] = d[key][ii]
    dras["assign"], ddecs["assign"] = get_tpos(
        tilera, tiledec, assign["TARGET_RA"], assign["TARGET_DEC"]
    )

    # AR WD
    keep = np.zeros(len(assign["TARGET_RA"]), dtype=bool)
    for mskkey, msk in zip(wd_mskkeys, wd_msks):
        keep |= (assign[mskkey] & yaml_masks[mskkey][msk]) > 0
    dras["wd"], ddecs["wd"] = get_tpos(
        tilera, tiledec, assign["TARGET_RA"][keep], assign["TARGET_DEC"][keep]
    )
    petals["wd"] = assign["PETAL_LOC"][keep]
    nassign["WD"] = keep.sum()
    # AR STD
    keep = np.zeros(len(assign["TARGET_RA"]), dtype=bool)
    for mskkey, msk in zip(std_mskkeys, std_msks):
        keep |= (assign[mskkey] & yaml_masks[mskkey][msk]) > 0
    dras["std"], ddecs["std"] = get_tpos(
        tilera, tiledec, assign["TARGET_RA"][keep], assign["TARGET_DEC"][keep]
    )
    petals["std"] = assign["PETAL_LOC"][keep]
    nassign["STD"] = keep.sum()

    return parent, assign, dras, ddecs, petals, nassign


def print_assgn_parent_stats(
    survey, parent, assign, log=Logger.get(), step="", start=time(),
):
    """
    Prints for each mask the number of parent and assigned targets, and also the fraction of assigned targets.
    
    Args:
        survey: survey (string; e.g. "sv1", "sv2", "sv3", "main")
        parent: dictionary for the parent target sample (output by get_parent_assign_quants())
        assign: dictionary for the assigned target sample (output by get_parent_assign_quants()) 
        log (optional, defaults to Logger.get()): Logger object                                                                                                                                            
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()
    """
    # AR YAML and WD and STD masks
    yaml_masks, _, _, _, _ = get_dt_masks(survey, log=log, step=step, start=start,)

    # AR stats : assigned / parent
    log.info("======= ASSIGNMENT STATISTICS : START =======")
    log.info("# MASKKEY\tMASK\tPARENT\tASSIGN\tFRACTION")
    for mskkey in list(yaml_masks.keys()):
        if survey.lower()[:2] == "sv":
            mskkey_orig = "{}_{}".format(survey.upper(), mskkey)
        else:
            mskkey_orig = mskkey
        for msk in yaml_masks[mskkey].names():
            nparent = ((parent[mskkey] & yaml_masks[mskkey][msk]) > 0).sum()
            nassign = ((assign[mskkey] & yaml_masks[mskkey][msk]) > 0).sum()
            if nparent == 0:
                frac = 0.0
            else:
                frac = nassign / nparent
            log.info(
                "{}\t{}\t{}\t{}\t{:.2f}".format(
                    mskkey_orig, msk, nparent, nassign, frac
                )
            )
    log.info("======= ASSIGNMENT STATISTICS : END =======")


def get_ext_coeffs(band):
    """
    Returns the extinction coefficient for a given band.

    Args:
        band: band name: "G", "R", "Z", "W1", or "W2" (string)

    Returns:
        ext: extinction coefficient (float)

    Note:
        https://www.legacysurvey.org/dr9/catalogs/#galactic-extinction-coefficients
    """
    exts = {"G": 3.214, "R": 2.165, "Z": 1.211, "W1": 0.184, "W2": 0.113}
    return exts[band]


def flux2mag(flux, band=None, ebv=None):
    """
    Converts a flux to a (optionally extinction-corrected) magnitude

    Args:
        flux: flux in Nmgy (np.array of floats)
        band (optional, defaults to None): band name: "G", "R", "Z", "W1", or "W2" (string)
        ebv (optional, defaults to None): EBV values (np.array of floats)

    Returns:
        mag: AB magnitudes (np.array of floats); extinction-corrected if band and ebv not None

    Notes:
        flux < 0 values are converted to NaN in magnitudes
    """
    keep = flux > 0
    mag = np.nan + np.zeros(len(flux))
    mag[keep] = 22.5 - 2.5 * np.log10(flux[keep])
    if ebv is not None:
        mag -= get_ext_coeffs(band) * ebv
    return mag


def qa_print_infos(
    ax,
    survey,
    program,
    faflavor,
    tileid,
    tilera,
    tiledec,
    obscon,
    rundate,
    parent,
    assign,
    log=Logger.get(),
    step="",
    start=time(),
):
    """
    Print general fiber assignment infos on the QA plot.
    
    Args:
        ax: pyplot object
        survey: "sv1", "sv2", "sv3" or "main" (string)
        program: "DARK", "BRIGHT", or "BACKUP" (string)
        faflavor: usually {survey}{program} in lower cases (string)
        tileid: tile TILEID (int)
        tilera: tile center R.A. in degrees (float)
        tiledec: tile center Dec. in degrees (float)
        obscon: tile allowed observing conditions (string; e.g. "DARK|GRAY|BRIGHT|BACKUP") 
        rundate: used rundate (string)
        parent: dictionary for the parent target sample (output by get_parent_assign_quants())
        assign: dictionary for the assigned target sample (output by get_parent_assign_quants()) 
        log (optional, defaults to Logger.get()): Logger object                                                                                                                                            
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()
    """
    # AR hard-setting the plotted tracers
    # AR TBD: handle secondaries
    trmskkeys, trmsks = get_qa_tracers(
        survey, program, log=Logger.get(), step="", start=time(),
    )

    # AR masks
    yaml_masks, wd_mskkeys, wd_msks, std_mskkeys, std_msks = get_dt_masks(
        survey, log=log, step=step, start=start,
    )

    # AR infos : general
    x, y, dy, fs = 0.05, 0.95, -0.1, 10
    for t in [
        "flavor={}".format(faflavor),
        "TILEID={:06d}".format(tileid),
        "RA,DEC={:.3f},{:.3f}".format(tilera, tiledec),
        "obscon={}".format(obscon),
        "rundate={}".format(rundate),
        "",
    ]:
        ax.text(x, y, t.expandtabs(), fontsize=fs, transform=ax.transAxes)
        y += dy
    # AR infos: wd/std + tracers
    xs = [0.05, 0.65, 0.95, 1.20]
    has = ["left", "right", "right", "right"]
    tracers = []
    for mskkey, msk in zip(wd_mskkeys + std_mskkeys, wd_msks + std_msks):
        n = ((assign[mskkey] & yaml_masks[mskkey][msk]) > 0).sum()
        tracers += [[msk, "{}".format(n), "", ""]]
    tracers += [["", "", "", ""], ["MASK", "ASSGN", "PARENT", "FAFRAC"]]
    for msk, mskkey in zip(trmsks, trmskkeys):
        nparent = ((parent[mskkey] & yaml_masks[mskkey][msk]) > 0).sum()
        n = ((assign[mskkey] & yaml_masks[mskkey][msk]) > 0).sum()
        tracers += [
            [msk, "{}".format(n), "{}".format(nparent), "{:.2f}".format(n / nparent),]
        ]
    tracers += [["", "", "", ""]]
    for tracer in tracers:
        for i in range(4):
            ax.text(
                xs[i],
                y,
                tracer[i].expandtabs(),
                fontsize=fs,
                ha=has[i],
                transform=ax.transAxes,
            )
        y += dy

    # AR infos: brightest target and assigned object
    # AR infos: taking a default 16, in case new programs are added
    magthresh = 16.0
    if program == "DARK":
        magthres = 16.0
    if program == "BRIGHT":
        magthresh = 15.0
    if program == "BACKUP":
        magthresh = 15.0
    for sample, d in zip(["parent", "assgn"], [parent, assign]):
        ax.text(
            0.05, y, "Min. {} mag ".format(sample), fontsize=fs, transform=ax.transAxes
        )
        y += dy
        for mag, lab in zip(
            [d["GAIA_PHOT_G_MEAN_MAG"], flux2mag(d["FIBERTOTFLUX_R"])],
            ["GAIA_PHOT_G_MEAN_MAG)", "min(LS-R-FIBTOTMAG)"],
        ):
            magmin, color = "-", "k"
            keep = (np.isfinite(mag)) & (mag > 0)
            if keep.sum() > 0:
                magmin = mag[keep].min()
                if magmin < magthresh:
                    magmin, color = "{:.1f}".format(magmin), "r"
                else:
                    magmin, color = "{:.1f}".format(magmin), "k"
            ax.text(
                0.05,
                y,
                "{} = {}".format(lab, magmin),
                fontsize=fs,
                color=color,
                transform=ax.transAxes,
            )
            y += dy
        y += dy


def qa_print_petal_infos(
    ax, petals, assign,
):
    """
    Print general the assigned SKY, BAD, WD, STD, TGT per petal on the QA plot.
    
    Args:
        ax: pyplot object
        petals: dictionary with PETAL_LOC (np.array of floats) for each of the assigned "sky", "bad", "wd", "std" subsamples
        assign: dictionary for the assigned target sample (output by get_parent_assign_quants()) 
    """
    # AR stats per petal
    xs = [0.05, 0.25, 0.45, 0.65, 0.85, 1.05]
    ts = ["PETAL", "NSKY", "NBAD", "NWD", "NSTD", "NTGT"]
    y, dy = 0.95, -0.1
    fs = 10
    for i in range(6):
        ax.text(xs[i], y, ts[i], fontsize=fs, ha="center", transform=ax.transAxes)
    y += dy
    for p in range(10):
        if (petals["std"] == p).sum() == 0:
            color = "r"
        else:
            color = "k"
        ts = [
            "{:.0f}".format(p),
            "{:.0f}".format((petals["sky"] == p).sum()),
            "{:.0f}".format((petals["bad"] == p).sum()),
            "{:.0f}".format((petals["wd"] == p).sum()),
            "{:.0f}".format((petals["std"] == p).sum()),
            "{:.0f}".format((assign["PETAL_LOC"] == p).sum()),
        ]
        for i in range(6):
            ax.text(
                xs[i],
                y,
                ts[i],
                color=color,
                fontsize=fs,
                ha="center",
                transform=ax.transAxes,
            )
        y += dy
    # AR stats for all petals
    ts = [
        "ALL",
        "{:.0f}".format(len(petals["sky"])),
        "{:.0f}".format(len(petals["bad"])),
        "{:.0f}".format(len(petals["wd"])),
        "{:.0f}".format(len(petals["std"])),
        "{:.0f}".format(np.isfinite(assign["PETAL_LOC"]).sum()),
    ]
    for i in range(6):
        ax.text(
            xs[i],
            y,
            ts[i],
            color=color,
            fontsize=fs,
            ha="center",
            transform=ax.transAxes,
        )


def get_viewer_cutout(
    tileid,
    tilera,
    tiledec,
    tmpoutdir=tempfile.mkdtemp(),
    width_deg=4,
    pixscale=10,
    dr="dr9",
):
    """
    Downloads a cutout of the tile region from legacysurvey.org/viewer.
    
    Args:
        tileid: tile TILEID (int)
        tilera: tile center R.A. (float)
        tiledec: tile center Dec. (float)
        tmpoutdir (optional, defaults to a temporary directory): temporary directory where
        width_deg (optional, defaults to 4): width of the cutout in degrees (float)
        pixscale (optional, defaults to 10): pixel scale of the cutout
        dr (optional, default do "dr9"): imaging data release

    Returns:
        img: output of mpimg.imread() reading of the cutout (np.array of floats)
    """
    # AR cutout
    size = int(width_deg * 3600.0 / pixscale)
    layer = "ls-{}".format(dr)
    tmpstr = 'wget -q -O {}tmp-{}.jpeg "http://legacysurvey.org/viewer-dev/jpeg-cutout/?layer={}&ra={:.5f}&dec={:.5f}&pixscale={:.0f}&size={:.0f}"'.format(
        tmpoutdir, tileid, layer, tilera, tiledec, pixscale, size
    )
    # print(tmpstr)
    os.system(tmpstr)
    try:
        img = mpimg.imread("{}tmp-{}.jpeg".format(tmpoutdir, tileid))
        os.remove("{}tmp-{}.jpeg".format(tmpoutdir, tileid))
    except:
        img = np.zeros((size, size, 3))
    return img


def mycmap(name, n, cmin=0, cmax=1):
    """
    Defines a quantised color scheme.
    
    Args:
        name: matplotlib colormap name (used through: matplotlib.cm.get_cmap(name)) (string)
        n: number of different colors to be in the color scheme (int)
        cmin (optional, defaults to 0): flooring "color-value" (float)
        cmax (optional, defaults to 1): ceiling "color-value" (float)

    Returns:
        The quantised color map.

    Notes:
        https://matplotlib.org/examples/api/colorbar_only.html
    """
    cmaporig = matplotlib.cm.get_cmap(name)
    mycol = cmaporig(np.linspace(cmin, cmax, n))
    cmap = matplotlib.colors.ListedColormap(mycol)
    cmap.set_under(mycol[0])
    cmap.set_over(mycol[-1])
    return cmap


def get_tpos(tilera, tiledec, ras, decs):
    """
    Computes the projected distance of a set of coordinates to a tile center.
    
    Args:
        tilera: tile center R.A. in degrees (float)
        tiledec: tile center Dec. in degrees (float)
        ras: R.A. in degrees (np.array of floats)
        decs: Dec. in degrees (np.array of floats)
                                                                                                                                                                                 
    Returns:
        dras: projected distance (degrees) to the tile center along R.A. (np.array of floats)
        ddecs: projected distance (degrees) to the tile center along Dec. (np.array of floats)
    """
    tsky = SkyCoord(ra=tilera * units.deg, dec=tiledec * units.deg, frame="icrs")
    sky = SkyCoord(ra=ras * units.deg, dec=decs * units.deg, frame="icrs")
    spho = tsky.spherical_offsets_to(sky)
    return spho[0].value, spho[1].value


def deg2pix(dras, ddecs, width_deg, width_pix):
    """
    Converts (dras,ddecs) to (xs,ys) in cutout img pixels.

    Args:
        dras: projected distance (degrees) along R.A. to the center of the cutout (np.array of floats)
        ddecs: projected distance (degrees) along Dec. to the center of the cutout (np.array of floats)
        width_deg: width of the cutout in degrees (np.array of floats)
        width_pix: width of the cutout in pixels (np.array of floats)

    Returns:
        dxs: distance (pixels) along x to the center of the cutout (np.array of floats)
        dys: distance (pixels) along y to the center of the cutout (np.array of floats)

    Notes:
        not sure at the <1 pixel level...
    """
    dxs = width_pix * (0.5 - dras / width_deg)
    dys = width_pix * (0.5 + ddecs / width_deg)
    return dxs, dys


def plot_cutout(
    ax,
    img,
    width_deg,
    dras,
    ddecs,
    dopetal=False,
    c="w",
    alpha=None,
    txts=None,
    xtxts=None,
    ytxts=None,
    vmin=None,
    vmax=None,
    cmap=mycmap("coolwarm", 10, 0, 1),
):
    """
    Plots a ls-dr9 cutout, with overlaying targets coordinates.

    Args:
        ax: pyplot object
        img: mpimg.imread(ls-dr9-cutout)
        width_deg: width of the cutout in degrees (np.array of floats)
        dras: targets projected distance (degrees) along R.A. to the center of the cutout (np.array of floats)
        ddecs: targets projected distance (degrees) along Dec. to the center of the cutout (np.array of floats)
        dopetal (optional, defaults to False): overplot petals? (boolean)
        c (optional, defaults to "w"): color used to display targets (string)
        alpha (optional, defaults to None): pyplot alpha
        txts (optional, defaults to None): list of text to display (list of strings)
        xtxts (optional, defaults to None): list normalized x-positions of text to display (list of strings)
        ytxts (optional, defaults to None): list normalized y-positions of text to display (list of strings)
        vmin (optional, defaults to None): minimum value for the colorbar
        vmax (optional, defaults to None): maximum value for the colorbar
        cmap (optional, defaults to mycmap("coolwarm", 10, 0, 1)): colormap scheme
    """
    # AR txts, xtxts, ytxts : lists
    # AR setting transparency as a function of density /deg2
    if (dras is not None) & (alpha is None):
        tmpdens = np.array([0, 100, 500, 1000, 5000, 7500, 1e10],)
        tmpalph = np.array([1, 0.8, 0.5, 0.2, 0.1, 0.05, 0.025])
        alpha = tmpalph[
            np.where(tmpdens > len(dras) / (np.pi * tile_radius_deg ** 2))[0][0]
        ]
    width_pix = img.shape[0]
    ax.imshow(
        img,
        origin="upper",
        zorder=0,
        extent=[0, width_pix, 0, width_pix],
        aspect="equal",
    )
    ax.set_aspect("equal")
    ax.set_xlim(-0.5, width_pix + 0.5)
    ax.set_ylim(-0.5, width_pix + 0.5)

    # AR data points
    if dras is not None:
        # AR rescaling degrees to img pixels ; not sure at <1 pixel...
        dxs, dys = deg2pix(dras, ddecs, width_deg, width_pix)
        if isinstance(c, str):
            ax.scatter(dxs, dys, c=c, s=1, alpha=alpha)
        else:
            ax.scatter(dxs, dys, c=c, s=1, alpha=alpha, vmin=vmin, vmax=vmax, cmap=cm)
    # AR per petal infos
    if dopetal:
        for ang, p in zip(
            np.linspace(2 * np.pi, 0, 11), [7, 8, 9, 0, 1, 2, 3, 4, 5, 6]
        ):
            dxs, dys = deg2pix(
                np.array([0, tile_radius_deg * np.cos(ang)]),
                np.array([0, tile_radius_deg * np.sin(ang)]),
                width_deg,
                width_pix,
            )
            ax.plot(
                dxs, dys, c="r", lw=0.25, alpha=1.0, zorder=1,
            )
            anglab = ang + 0.1 * np.pi
            dxs, dys = deg2pix(
                1.1 * tile_radius_deg * np.cos(anglab),
                1.1 * tile_radius_deg * np.sin(anglab),
                width_deg,
                width_pix,
            )
            ax.text(
                dxs, dys, "{:.0f}".format(p), color="r", va="center", ha="center",
            )

    ax.axis("off")
    if txts is not None:
        for txt, xtxt, ytxt in zip(txts, xtxts, ytxts):
            ax.text(
                xtxt,
                ytxt,
                txt,
                color="y",
                fontweight="bold",
                fontsize=10,
                ha="center",
                va="top",
                transform=ax.transAxes,
            )


def plot_hist(ax, mags, magps, msk):
    """
    Plots a normalized histogram for the assigned magnitudes (xs) and the parent magnitudes (xps).

    Args:
        ax: pyplot object
        mags: assigned magnitudes (np.array of floats)
        magps : parent magnitudes (np.array of floats)
        msk: mask name of the plotted sample
    """
    #
    selp = np.isfinite(magps)
    sel = np.isfinite(mags)
    bins = np.linspace(magps[selp].min(), magps[selp].max(), 26)
    #
    cps, _, _ = ax.hist(
        magps[selp],
        bins=bins,
        histtype="step",
        alpha=0.3,
        lw=3,
        color="k",
        density=False,
        label="{} parent ({})".format(msk, len(magps)),
    )
    cs, _, _, = ax.hist(
        mags[sel],
        bins=bins,
        histtype="step",
        alpha=1.0,
        lw=1.0,
        color="k",
        density=False,
        label="{} assigned ({})".format(msk, len(mags)),
    )
    ax.set_ylabel("counts")
    ax.grid(True)
    # ax.legend(loc=2)
    axr = ax.twinx()
    axr.plot(
        0.5 * (bins[1:] + bins[:-1]),
        np.array(cs) / np.array(cps).astype(float),
        color="r",
        lw=0.5,
    )
    axr.yaxis.label.set_color("r")
    axr.tick_params(axis="y", colors="r")
    axr.set_ylabel("ratio", labelpad=0)
    axr.set_ylim(0, 1)
    txts = [msk, "assigned/parent = {}/{}".format(len(mags), len(magps))]
    xtxts = [0.5, 0.5]
    ytxts = [0.98, 0.90]
    for txt, xtxt, ytxt in zip(txts, xtxts, ytxts):
        ax.text(
            xtxt,
            ytxt,
            txt,
            color="k",
            fontweight="bold",
            fontsize=10,
            ha="center",
            va="top",
            transform=ax.transAxes,
        )


def get_qa_farange(fafrac, dfa=0.2):
    """
    Picks the plotted fiber assignment rate range for the QA plot.

    Args:
        fafrac: fiber assignment rate for the plotted sample (float)
        dfa (optional, defaults to 0.2): plotted range (float)

    Returns:
        famin: lower boundary of the plotted fiber assignment rate (float)
        famax: upper boundary of the plotted fiber assignment rate (float)
    """
    famin = np.max([0, np.round(fafrac - dfa / 2, 1)])
    famax = np.min([1, np.round(fafrac + dfa / 2, 1)])
    return famin, famax


def plot_hist_tracer(ax, survey, parent, assign, msk, mskkey):
    """
    Plots a normalized histogram for the assigned magnitudes (xs) and the parent magnitudes (xps).

    Args:
        ax: pyplot object
        survey: survey (string; e.g. "sv1", "sv2", "sv3", "main")
        parent: dictionary for the parent target sample (output by get_parent_assign_quants())                                                                                   
        assign: dictionary for the assigned target sample (output by get_parent_assign_quants()) 
        msk: mask name of the plotted sample (string)
        mskkey: key to select the mask on (string)
    """
    # AR YAML mask dictionary
    yaml_masks, _, _, _, _ = get_dt_masks(survey)

    # AR selecting the relevant tracer
    if mskkey in list(parent.keys()):
        mskpsel = (parent[mskkey] & yaml_masks[mskkey][msk]) > 0
    else:
        mskpsel = np.zeros(len(parent["TARGETID"]), dtype=bool)
    # AR if no parent target, just skip
    if mskpsel.sum() > 0:
        msksel = (assign[mskkey] & yaml_masks[mskkey][msk]) > 0
        famin, famax = get_qa_farange(msksel.sum() / float(mskpsel.sum()))
        # AR mag hist
        band = "R"
        if "MWS" in msk:
            band = "R"
        if "BGS" in msk:
            band = "R"
        if "LRG" in msk:
            band = "Z"
        if "ELG" in msk:
            band = "G"
        if "QSO" in msk:
            band = "R"
        #
        dohist = 0
        # AR if ls-dr9 flux is here, we plot that
        if ((parent["FLUX_{}".format(band)] > 0) & (mskpsel)).sum() > 0:
            dohist = 1
            # AR parent
            magp = flux2mag(
                parent["FLUX_{}".format(band)][mskpsel],
                band=band,
                ebv=parent["EBV"][mskpsel],
            )
            # AR assign
            mag = flux2mag(
                assign["FLUX_{}".format(band)][msksel],
                band=band,
                ebv=assign["EBV"][msksel],
            )
            # AR xlabel
            ax.set_xlabel(
                "22.5 - 2.5*log10(FLUX_{}) - {:.3f} * EBV".format(
                    band, get_ext_coeffs(band)
                )
            )
        # AR if no ls-dr9 flux, we try gaia_g
        elif ((np.isfinite(parent["GAIA_PHOT_G_MEAN_MAG"])) & (mskpsel)).sum() > 0:
            dohist = 1
            magp = parent["GAIA_PHOT_G_MEAN_MAG"][mskpsel]
            mag = assign["GAIA_PHOT_G_MEAN_MAG"][msksel]
            ax.set_xlabel("GAIA_PHOT_G_MEAN_MAG")
        if dohist == 1:
            plot_hist(ax, mag, magp, msk)
            _, ymax = ax.get_ylim()
            ax.set_ylim(0.8, 100 * ymax)
            ax.set_yscale("log")


def plot_colcol_tracer(
    ax,
    xbands,
    ybands,
    survey,
    parent,
    assign,
    msk,
    mskkey,
    xlim,
    ylim,
    gridsize=20,
    cm=mycmap("coolwarm", 10, 0, 1),
):
    """
    Plots a color-color diagram, with color-coding with the fiber assignment rate,
        and transparency-coding the density.

    Args:
        ax: pyplot object
        survey: survey (string; e.g. "sv1", "sv2", "sv3", "main")
        xbands: two-elements list, the x-axis color being xbands[1] - xbands[0] (list of strings)
        ybands: two-elements list, the y-axis color being ybands[1] - ybands[0] (list of strings)
        parent: dictionary for the parent target sample (output by get_parent_assign_quants())                                                                                   
        assign: dictionary for the assigned target sample (output by get_parent_assign_quants()) 
        msk: mask name of the plotted sample (string)
        mskkey: key to select the mask on (string)
        xlim: plt.xlim
        ylim: plt.ylim
        gridsize (optional, defaults to 20): plt.hexbin gridsize parameter (int)
        cmap (optional, defaults to mycmap("coolwarm", 10, 0, 1)): colormap scheme
    """
    # AR YAML mask dictionary
    yaml_masks, _, _, _, _ = get_dt_masks(survey)

    # AR selecting the relevant tracer
    if mskkey in list(parent.keys()):
        mskpsel = (parent[mskkey] & yaml_masks[mskkey][msk]) > 0
    else:
        mskpsel = np.zeros(len(parent["TARGETID"]), dtype=bool)

    # AR plotting if some parent objects with valid colors
    keep = mskpsel.copy()
    for band in [xbands[0], xbands[1], ybands[0], ybands[1]]:
        keep &= parent["FLUX_{}".format(band)] > 0

    if keep.sum() > 0:
        # AR
        msksel = (assign[mskkey] & yaml_masks[mskkey][msk]) > 0

        # AR using a dictionary
        tmpdict = {"parent": {}, "assign": {}}
        for sample_name, sample, sel in zip(
            ["parent", "assign"], [parent, assign], [mskpsel, msksel]
        ):
            for axis_name, bands in zip(["x", "y"], [xbands, ybands]):
                mag0 = flux2mag(
                    sample["FLUX_{}".format(bands[0])][sel],
                    bands[0],
                    sample["EBV"][sel],
                )
                mag1 = flux2mag(
                    sample["FLUX_{}".format(bands[1])][sel],
                    bands[1],
                    sample["EBV"][sel],
                )
                tmpdict[sample_name][axis_name] = mag0 - mag1

        # AR first getting the hexbin outputs
        hbp = ax.hexbin(
            tmpdict["parent"]["x"],
            tmpdict["parent"]["y"],
            C=None,
            gridsize=gridsize,
            extent=(xlim[1], xlim[0], ylim[0], ylim[1]),
            mincnt=0,
            visible=False,
        )
        hb = ax.hexbin(
            tmpdict["assign"]["x"],
            tmpdict["assign"]["y"],
            C=None,
            gridsize=gridsize,
            extent=(xlim[1], xlim[0], ylim[0], ylim[1]),
            mincnt=0,
            visible=False,
        )

        # AR restricting to pixels with some parent data
        keep = hbp.get_array() > 0
        tmpx = hb.get_offsets()[keep, 0]
        tmpy = hb.get_offsets()[keep, 1]
        tmpc = hb.get_array()[keep]
        tmpcp = hbp.get_array()[keep].astype(float)

        # AR fraction assigned, clipped to famin,famax
        fafrac = msksel.sum() / float(mskpsel.sum())
        famin, famax = get_qa_farange(fafrac)
        c = cm(np.clip(((tmpc / tmpcp) - famin) / (famax - famin), 0, 1))

        # AR transparency = f(nb of parent obj)
        tmpmin, tmpmax = (
            1,
            1.2 * tmpcp.sum() / float(len(hbp.get_array())),
        )
        c[:, 3] = np.clip((tmpcp - tmpmin) / (tmpmax - tmpmin), 0, 1)
        sc = ax.scatter(tmpx, tmpy, c=c, s=15,)
        sc.cmap = cm
        ax.set_xlabel("{} - {}".format(xbands[0].lower(), xbands[1].lower()))
        ax.set_ylabel("{} - {}".format(ybands[0].lower(), ybands[1].lower()))
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.grid(True)
        ax.text(
            0.5,
            0.93,
            msk,
            color="k",
            fontweight="bold",
            fontsize=10,
            ha="center",
            transform=ax.transAxes,
        )
        cbar = plt.colorbar(sc)
        cbar.set_label("fraction assigned")
        cbar.mappable.set_clim(famin, famax)


def plot_sky_fa(
    axs,
    img,
    survey,
    parent,
    assign,
    dras,
    ddecs,
    msk,
    mskkey,
    width_deg,
    gridsize=30,
    cm=mycmap("coolwarm", 10, 0, 1),
):
    """
    Plots the sky distribution of the parent sample, the assigned sample, and of the fiber assignment rate.

    Args:
        axs: list of 3 pyplot objects, respectively for the parent sample, the assigned sample, and the fiber assignment rate
        img: mpimg.imread(ls-dr9-cutout)
        survey: survey (string; e.g. "sv1", "sv2", "sv3", "main")
        parent: dictionary for the parent target sample (output by get_parent_assign_quants())                                                                                   
        assign: dictionary for the assigned target sample (output by get_parent_assign_quants()) 
        dras: dictionary with projected distance (degrees) along R.A. to the center of the tile (np.array of floats),                                                            
            for each of the following subsample: "parent", "assign", "sky", "bad", "wd", "std" (all assigned subsamples,
            except parent)
        ddecs: same as dras, for projected distances along Dec.
        msk: mask name of the plotted sample (string)
        mskkey: key to select the mask on (string)
        width_deg: width of the cutout in degrees (np.array of floats)
        gridsize (optional, defaults to 30): plt.hexbin gridsize parameter (int)
        cmap (optional, defaults to mycmap("coolwarm", 10, 0, 1)): colormap scheme
    """
    # AR YAML mask dictionary
    yaml_masks, _, _, _, _ = get_dt_masks(survey)

    # AR selecting the relevant tracer
    if mskkey in list(parent.keys()):
        mskpsel = (parent[mskkey] & yaml_masks[mskkey][msk]) > 0
    else:
        mskpsel = np.zeros(len(parent["TARGETID"]), dtype=bool)

    if mskpsel.sum() > 0:
        # AR assign sample tracer selection
        msksel = (assign[mskkey] & yaml_masks[mskkey][msk]) > 0

        # AR xlim, ylim
        xlim = (width_deg / 2, -width_deg / 2)
        ylim = (-width_deg / 2, width_deg / 2)

        # AR area of the plotting window in deg2
        plot_area = (xlim[0] - xlim[1]) * (ylim[1] - ylim[0])

        # AR parent
        plot_cutout(
            axs[0],
            img,
            width_deg,
            dras["parent"][mskpsel],
            ddecs["parent"][mskpsel],
            dopetal=True,
            txts=[
                msk,
                "parent : {:.0f}".format(mskpsel.sum() / tile_area) + r" deg$^{-2}$",
            ],
            xtxts=[0.5, 0.5],
            ytxts=[0.98, 0.1],
        )

        # AR assigned
        plot_cutout(
            axs[1],
            img,
            width_deg,
            dras["assign"][msksel],
            ddecs["assign"][msksel],
            dopetal=True,
            txts=[
                msk,
                "assigned : {:.0f}".format(msksel.sum() / tile_area) + r" deg$^{-2}$",
            ],
            xtxts=[0.5, 0.5],
            ytxts=[0.98, 0.1],
        )

        # AR fraction assigned, clipped to famin,famax
        fafrac = msksel.sum() / float(mskpsel.sum())
        famin, famax = get_qa_farange(fafrac)
        txts = [msk, r"mean = {:.2f}".format(fafrac)]
        xtxts = [0.5, 0.5]
        ytxts = [0.93, 0.03]

        # AR assigned fraction
        ax = axs[2]
        x = dras["parent"][mskpsel]
        y = ddecs["parent"][mskpsel]
        C = np.in1d(parent["TARGETID"][mskpsel], assign["TARGETID"][msksel])
        hb = ax.hexbin(
            x,
            y,
            C=C,
            gridsize=gridsize,
            extent=(xlim[1], xlim[0], ylim[0], ylim[1]),
            mincnt=1,
            alpha=0.5,
            vmin=famin,
            vmax=famax,
        )
        hb.cmap = cm
        ax.set_xlabel(r"$\Delta$RA [deg.]")
        ax.set_ylabel(r"$\Delta$DEC [deg.]")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.grid(True)
        for txt, xtxt, ytxt in zip(txts, xtxts, ytxts):
            ax.text(
                xtxt,
                ytxt,
                txt,
                color="k",
                fontweight="bold",
                fontsize=10,
                ha="center",
                transform=ax.transAxes,
            )
        cbar = plt.colorbar(hb)
        cbar.set_label("fraction assigned")
        cbar.mappable.set_clim(famin, famax)


def make_qa(
    outpng,
    survey,
    program,
    faflavor,
    targfns,
    fiberassignfn,
    tileid,
    tilera,
    tiledec,
    obscon,
    rundate,
    tmpoutdir=tempfile.mkdtemp(),
    width_deg=4,
    log=Logger.get(),
    step="",
    start=time(),
):
    """
    Make fba_launch QA plot.
    
    Args:
        outpng: written output PNG file (string)
        survey: "sv1", "sv2", "sv3" or "main" (string)
        program: "DARK", "BRIGHT", or "BACKUP" (string)
        faflavor: usually {survey}{program} in lower cases (string)
        fiberassignfn: path to the output fiberassign-TILEID.fits file (string)
        tileid: tile TILEID (int)
        tilera: tile center R.A. in degrees (float)
        tiledec: tile center Dec. in degrees (float)
        obscon: tile allowed observing conditions (string; e.g. "DARK|GRAY|BRIGHT|BACKUP") 
        rundate: used rundate (string)
        tmpoutdir (optional, defaults to a temporary directory): temporary directory (to download the cutout)
        width_deg (optional, defaults to 4): width of the cutout in degrees (np.array of floats)
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()
    """
    log.info("")
    log.info("")
    log.info("{:.1f}s\t{}\tTIMESTAMP={}".format(time() - start, step, Time.now().isot))
    # AR WD and STD used masks
    _, wd_mskkeys, wd_msks, std_mskkeys, std_msks = get_dt_masks(
        survey, log=log, step=step, start=start,
    )
    log.info("{:.1f}s\t{}\twd_mskkeys = {}".format(time() - start, step, wd_mskkeys))
    log.info("{:.1f}s\t{}\twd_msks = {}".format(time() - start, step, wd_msks))
    log.info("{:.1f}s\t{}\tstd_mskkeys = {}".format(time() - start, step, std_mskkeys))
    log.info("{:.1f}s\t{}\tstd_msks = {}".format(time() - start, step, std_msks))

    # AR plotted tracers
    # AR TBD: handle secondary?
    trmskkeys, trmsks = get_qa_tracers(
        survey, program, log=log, step="doplot", start=start,
    )

    # AR storing parent/assigned quantities
    parent, assign, dras, ddecs, petals, nassign = get_parent_assign_quants(
        survey, targfns, fiberassignfn, tilera, tiledec, log=log, step=step, start=start
    )

    # AR stats : assigned / parent
    print_assgn_parent_stats(survey, parent, assign, log=log, step=step, start=start)

    # AR start plotting
    fig = plt.figure(figsize=(30, 3 * (1 + len(trmsks))))
    gs = gridspec.GridSpec(1 + len(trmsks), 7, wspace=0.5, hspace=0.3)

    # AR overall infos
    ax = plt.subplot(gs[0, 0])
    ax.axis("off")
    # AR infos : general
    qa_print_infos(
        ax,
        survey,
        program,
        faflavor,
        tileid,
        tilera,
        tiledec,
        obscon,
        rundate,
        parent,
        assign,
        log=log,
        step=step,
        start=start,
    )

    # AR stats per petal
    ax = plt.subplot(gs[0, 1])
    ax.axis("off")
    qa_print_petal_infos(
        ax, petals, assign,
    )

    # AR cutout
    img = get_viewer_cutout(
        tileid,
        tilera,
        tiledec,
        tmpoutdir=tmpoutdir,
        width_deg=width_deg,
        pixscale=10,
        dr="dr9",
    )

    # AR SKY, BAD, WD, STD, TGT
    iys = [2, 3, 4, 5, 6]
    keys = ["sky", "bad", "wd", "std", "assign"]
    txts = ["SKY", "BAD", "WD", "STD", "TGT"]
    alphas = [0.25, 1.0, 1.0, 1.0, 0.025]
    for iy, key, txt, alpha in zip(iys, keys, txts, alphas):
        ax = fig.add_subplot(gs[0, iy])
        plot_cutout(
            ax,
            img,
            width_deg,
            dras[key],
            ddecs[key],
            dopetal=True,
            alpha=alpha,
            txts=[txt],
            xtxts=[0.2],
            ytxts=[0.98],
        )

    # AR looping on tracers
    ix = 1
    for msk, mskkey in zip(trmsks, trmskkeys):

        # AR parent and assign magnitude distributions
        plot_hist_tracer(plt.subplot(gs[ix, 1]), survey, parent, assign, msk, mskkey)

        # AR color-color diagram, with fiber assignment rate color-coded, and density transparency-coded
        gridsize = 20
        for iy, xbands, ybands, xlim, ylim in zip(
            [2, 3],
            [("R", "Z"), ("R", "Z")],
            [("G", "R"), ("R", "W1")],
            [(-0.5, 2.5), (-0.5, 2.5)],
            [(-0.5, 2.5), (-2, 5)],
        ):
            ax = plt.subplot(gs[ix, iy])
            plot_colcol_tracer(
                ax,
                xbands,
                ybands,
                survey,
                parent,
                assign,
                msk,
                mskkey,
                xlim,
                ylim,
                gridsize=20,
            )

        # AR position in tile
        axs = [plt.subplot(gs[ix, 4]), plt.subplot(gs[ix, 5]), plt.subplot(gs[ix, 6])]
        plot_sky_fa(
            axs, img, survey, parent, assign, dras, ddecs, msk, mskkey, width_deg
        )

        #
        ix += 1

    #  AR saving plot
    plt.savefig(
        outpng, bbox_inches="tight",
    )
    plt.close()


def rmv_nonsvn(mytmpouts, myouts, log=Logger.get(), step="", start=time()):
    """
    Remove fba_launch non-SVN products
    
    Args:
        mytmpouts: dictionary with the temporary files location (dictionary);
            concerns the following keys:
                "tiles", "sky", "gfa", "targ", "scnd", "too", "fba"
        myouts: dictionary with the fba_launch args.outdir location (dictionary);
            contains same keys as mytmpouts
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()
    """
    log.info("")
    log.info("")
    log.info("{:.1f}s\t{}\tTIMESTAMP={}".format(time() - start, step, Time.now().isot))
    for key in ["tiles", "sky", "gfa", "targ", "scnd", "too", "fba"]:
        if os.path.isfile(mytmpouts[key]):
            os.remove(mytmpouts[key])
            log.info(
                "{:.1f}s\t{}\tdeleting file {}".format(
                    time() - start, step, mytmpouts[key]
                )
            )


def mv_temp2final(mytmpouts, myouts, doclean, log=Logger.get(), step="", start=time()):
    """
    Moves the fba_launch outputs from the temporary location to the args.outdir location.   

    Args:
        mytmpouts: dictionary with the temporary files location (dictionary);
            contains the following keys:
                "fiberassign", "png", "log"
                if doclean="n", also: "tiles", "sky", "gfa", "targ", "scnd", "too", "fba"
        myouts: dictionary with the fba_launch args.outdir location (dictionary);
            contains same keys as mytmpouts
        doclean: remove non-SVN files? "y" or "n" (string);
            args.clean fba_launch argument
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()

    Notes:
        actually, the log is not moved here; it is moved in fba_launch, after the main()
    """
    log.info("")
    log.info("")
    log.info("{:.1f}s\t{}\tTIMESTAMP={}".format(time() - start, step, Time.now().isot))
    # AR optional files
    if doclean == "n":
        for key in ["tiles", "sky", "gfa", "targ", "scnd", "too", "fba"]:
            if os.path.isfile(mytmpouts[key]):
                _ = shutil.move(mytmpouts[key], myouts[key])
                log.info(
                    "{:.1f}s\t{}\tmoving file {} to {}".format(
                        time() - start, step, mytmpouts[key], myouts[key]
                    )
                )
            else:
                log.info(
                    "{:.1f}s\t{}\tno file {}".format(
                        time() - start, step, mytmpouts[key]
                    )
                )
    # AR SVN-checked files
    for key in ["fiberassign", "png"]:
        _ = shutil.move(mytmpouts[key], myouts[key])
        log.info(
            "{:.1f}s\t{}\tmoving {} to {}".format(
                time() - start, step, mytmpouts[key], myouts[key]
            )
        )
    # AR actually here not moving the log file
    # AR moving it after the main() in fba_launch
    key = "log"
    log.info(
        "{:.1f}s\t{}\tmoving {} to {}".format(
            time() - start, step, mytmpouts[key], myouts[key]
        )
    )
