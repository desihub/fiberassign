"""
fiberassign.fba_launch_io
=============
Utility functions for fba_launch
"""
from __future__ import absolute_import, division

# system
import os
import sys
import tempfile

# time
from time import time
from datetime import datetime, timezone

#
import numpy as np

# astropy
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

# desitarget
from desitarget.geomask import match_to

# fiberassign
from fiberassign.utils import Logger

# fiberassign.fba_launch_io
from fiberassign.fba_launch_io import (
    print_config_infos,
    get_desitarget_paths,
    create_tile,
    create_sky,
    create_mtl,
    create_too,
)


def fba_rerun_get_settings(
    fn, bugfix=True, log=Logger.get(), step="settings", start=time(),
):
    """
    Get the required settings from a fiberassign-TILEID.fits.gz file to rerun the assignment,
        including the intermediate files (TILEID-{tiles,sky,targ,scnd,too}.fits).

    Args:
        fn: full path to the fiberassign-TILEID.fits.gz file to be rerun (string)
        bugfix (optional, defaults to True): bugfix the arguments with fba_rerun_bugfix_settings()? (boolean)
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start (optional, defaults to time()): start time for log (in seconds; output of time.time())

    Returns:
        mydict: dictionary with the required various settings (dictionary);
                the keys are:
                    fa_ver, skybrver, survey, program, tileid, tilera, tiledec, obscon,
                    dr, gaiadr, dtver, pmcorr, pmtime_utc_str,
                    rundate, mtltime, fieldrot, standards_per_petal, sky_per_petal
                    dotoo, mjd_min, mjd_max;
                if fa_ver >= 2.4.0: also sky_per_slitblock;
                if fa_ver >= 3.0.0: also ha, margin-pos, margin-petal, margin-gfa.

    Notes:
        Only runs for a SV3 or Main BRIGHT/DARK tile; exists with error otherwise.
        20211119 : now also allows Main BACKUP tile
        20211119 : add lookup_sky_source and SKYHEALPIXS_DIR
    """
    #
    hdr = fits.getheader(fn, 0)
    # AR SV3/Main BRIGHT/DARK?
    faflavors = [
        "sv3bright", "sv3dark",
        "mainbright", "mainbright1b",
        "maindark", "maindark1b",
        "mainbackup",
    ]
    if hdr["FAFLAVOR"] not in faflavors:
        msg = "{:.1f}s\t{}\tFAFLAVOR={} not in {}; exiting".format(
            time() - start, step, hdr["FAFLAVOR"], faflavors
        )
        log.error(msg)
        raise ValueError(msg)

    mydict = {}

    # AR from header: fa_ver, obscon, fieldrot, tileid, tilera, tiledec, mtltime
    # AR mtltime: args.mtltime has been introduced in fiberassign/2.4.0 only
    for key in ["fa_ver", "obscon", "fieldrot", "tileid", "tilera", "tiledec", "mtltime"]:
        mydict[key] = hdr[key.upper()]

    # AR settings we store
    keys = [
        "survey",
        "program",
        "rundate",
        "dr",
        "gaiadr",
        "dtver",
        "pmcorr",
        "pmtime_utc_str",
        "standards_per_petal",
        "sky_per_petal",
    ]
    if mydict["fa_ver"] >= "2.4.0":
        keys += ["sky_per_slitblock"]
    if mydict["fa_ver"] >= "3.0.0":
        keys += ["ha"]
        keys += ["margin-pos", "margin-petal", "margin-gfa"]
    if mydict["fa_ver"] > "5.3.0":
        keys += ["lookup_sky_source"]
    # AR storing the FAARGS in a dictionary
    faargs = np.array(hdr["FAARGS"].split())
    for key in keys:
        if key[:7] == "margin-":
            ii = np.where(faargs == "--{}".format(key.replace("margin-", "margin_")))[0]
        else:
            ii = np.where(faargs == "--{}".format(key))[0]
        if len(ii) > 0:
            mydict[key] = faargs[ii[0] + 1]
        # AR args.pmtime_utc_str initially names args.pmtime
        elif (key == "pmtime_utc_str") & ("--pmtime" in faargs):
            ii = np.where(faargs == "--pmtime")[0]
            mydict["pmtime_utc_str"] = faargs[ii[0] + 1]
        # AR should not be other cases
        else:
            log.error(
                "{:.1f}s\t{}\tkey={} not present in FAARGS, and not in expected missing keys; exiting".format(
                    time() - start, step, key,
                )
            )
            sys.exit(1)

    # AR fieldrot, ha: protecting against negative values,
    # AR    which cause problem as e.g. in "fba_run --fieldrot -1.0"
    # AR    changing those to a string is not a problem
    # AR    though, needs to be in single quotes, as this is written
    # AR    in a bash string encapsulated in double quotes
    for key in ["fieldrot", "ha"]:
        if key in mydict:
            mydict[key] = "' {}'".format(mydict[key])

    # AR SKYBRICKS_DIR
    mydict["skybrver"] = "-"  # default value if not set
    keys = [cards[0] for cards in hdr.cards]
    vals = [
        hdr[key.replace("NAM", "VER")].split("/")[-1]
        for key in keys
        if hdr[key] == "SKYBRICKS_DIR"
    ]
    if len(vals) > 0:
        mydict["skybrver"] = vals[0]

    # AR SKYHEALPIXS_DIR
    mydict["skyhpver"] = "-"  # default value if not set
    keys = [cards[0] for cards in hdr.cards]
    vals = [
        hdr[key.replace("NAM", "VER")].split("/")[-1]
        for key in keys
        if hdr[key] == "SKYHEALPIXS_DIR"
    ]
    if len(vals) > 0:
        mydict["skyhpver"] = vals[0]

    # AR ToO
    # AR default: run ToOs, with MJD_BEGIN < mjd_now < MJD_END
    mydict["dotoo"] = True
    try:
        mjd_now = Time(
            datetime.strptime(mydict["pmtime_utc_str"], "%Y-%m-%dT%H:%M:%S%z")
        ).mjd
    except ValueError:
        mjd_now = Time(
            datetime.strptime(mydict["pmtime_utc_str"], "%Y-%m-%dT%H:%M:%S")
        ).mjd
    mydict["mjd_min"] = mjd_now
    mydict["mjd_max"] = mjd_now

    # AR bugfix the arguments?
    if bugfix:
        mydict = fba_rerun_bugfix_settings(mydict, log=log, step=step ,start=start)

    return mydict


def fba_rerun_bugfix_settings(
    mydict, log=Logger.get(), step="settings", start=time(),
):
    """
    Fix the arguments read in fba_rerun_in_settings() to allow reproducibility.

    Args:
        mydict: dictionary with the required various settings (dictionary);
                the keys are:
                    fa_ver, skybrver, survey, program, tileid, tilera, tiledec, obscon,
                    dr, gaiadr, dtver, pmcorr, pmtime_utc_str,
                    rundate, mtltime, fieldrot, standards_per_petal, sky_per_petal
                    dotoo, mjd_min, mjd_max;
                if fa_ver >= 2.4.0: also sky_per_slitblock;
                if fa_ver >= 3.0.0: also ha, margin-pos, margin-petal, margin-gfa.
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start (optional, defaults to time()): start time for log (in seconds; output of time.time())

    Returns:
        mydict: dictionary with the bug-fixed values (dictionary)

    Notes:
        Fixes the following bugs:
            - rundate fix for 19 SV3 tiles designed on 2021-04-10
            - mtltime fix for SV3 TILEID=315 and 441
            - dtver for 2021 pre-shutdown Main (NIGHT<=20210518): 1.0.0 -> 1.1.1
            - fa_ver:
                - fa_ver = 2.2.0.dev2811 and rundate <  2021-04-14 -> fa_ver = 2.2.0
                - fa_ver = 2.2.0.dev2811 and rundate >= 2021-04-14 -> fa_ver = 2.3.0
                - fa_ver = 2.3.0.dev2838 -> fa_ver = 2.3.0
            - SKYBRICKS_DIR:
                - fa_ver = 2.4.0 -> skybrver = "v2" (SKYBRICKS_DIR was not recorded yet in the header)
            - ToO MJD window:
                - fa_ver <= 5.1.1 -> mjd_min = mjd_now - 1, mjd_max = mjd_now + 30
            - ToO dotoo:
                - pmtime_utc_str < 2021-04-19T07:00:00+00 -> dotoo=False (ToOs were introduced on Apr. 19th, 2021)
        The SUBPRIORITY overwritting for SV3 is handled with the desitarget code version
            used to rerun_intermediate()
    """
    # AR case of 19 SV3 tiles designed on 2021-04-10:
    # AR    those have rundate = 2021-04-10T21:28:37
    # AR    there has been a fp update on 2021-04-10T20:00:39+00:00
    # AR    but apparently the desi-state*ecsv file at NERSC was not updated then
    # AR    so the original fiberassign run considered the previous fp state
    # AR    =>
    # AR    we manually modify the rundate, pmtime, mtltime
    if mydict["tileid"] in [
        4,
        30,
        58,
        84,
        112,
        139,
        166,
        193,
        220,
        248,
        275,
        302,
        329,
        356,
        364,
        383,
        391,
        410,
        418,
    ]:
        fixed_time = "2021-04-10T20:00:00"
        for key in ["pmtime_utc_str", "rundate", "mtltime"]:
            log.info(
                "{:.1f}s\t{}\tmodifying mydict[{}] from {} to {} (TILEID={} designed when desi-state*ecsv file not updated at NERSC".format(
                    time() - start, step, key, mydict[key], fixed_time, mydict["tileid"],
                )
            )
            mydict[key] = fixed_time

    # AR SV3 TILEID=315, 441
    # AR    for both, in the original files:
    # AR    - MTLTIME = 2021-04-22T18:55:39+00:00
    # AR    - the latest MTL TIMESTAMP in the original target files is: 2021-04-19T20:38:52 (315) and 2021-04-20T18:42:10 (441)
    # AR    however the MTL ledgers have been updated at 2021-04-22T17:09:xx (including TILEID=314, 440)
    # AR    it is likely the ledgers files were not svn-updated then.
    # AR    so we modify the MTLTIME to 2021-04-19T21:00:00 (315) and 2021-04-20T19:00:00 (441)
    # AR    (ie after the previous MTL update, but before the 314, 440 MTL updates)
    if mydict["tileid"] in [315, 441]:
        if mydict["tileid"] == 315:
            fixed_time = "2021-04-19T21:00:00+00:00"
        if mydict["tileid"] == 441:
            fixed_time = "2021-04-20T19:00:00+00:00"
        log.info(
            "{:.1f}s\t{}\tmodifying mydict[mtltime] from {} to {} (TILEID={} designed when MTL ledgers not updated at NERSC".format(
                time() - start, step, mydict["mtltime"], fixed_time, mydict["tileid"],
            )
        )
        mydict["mtltime"] = fixed_time

    # AR dtver: pre-shutdown main, replace 1.0.0 by 1.1.1
    # AR because main-ledgers are built from 1.1.1
    if (mydict["survey"] == "main") & (mydict["dtver"] == "1.0.0"):
        log.info(
            "{:.1f}s\t{}\tmodifying mydict[dtver]=1.0.0 to mydict[dtver]=1.1.1".format(
                time() - start, step
            )
        )
        mydict["dtver"] = "1.1.1"

    # AR fiberassign code version
    # AR SV3 handling, as SV3 was run with fiberassign/master for RUNDATE<2021-04-30
    # AR on 2021-04-13: https://github.com/desihub/fiberassign/pull/321
    # AR                this PR changes the assignment
    # AR                so we use 2.2.0 before, 2.3.0 after
    fixed_faver = None
    if (mydict["fa_ver"] == "2.2.0.dev2811") & (mydict["rundate"] < "2021-04-14"):
        fixed_faver = "2.2.0"
    elif (mydict["fa_ver"] == "2.2.0.dev2811") & (mydict["rundate"] >= "2021-04-14"):
        fixed_faver = "2.3.0"
    elif mydict["fa_ver"] == "2.3.0.dev2838":
        fixed_faver = "2.3.0"
    if fixed_faver is not None:
        log.info(
            "{:.1f}s\t{}\tmodifying mydict[fa_ver]={} to mydict[fa_ver]={}".format(
                time() - start, step, mydict["fa_ver"], fixed_faver,
            )
        )
        mydict["fa_ver"] = fixed_faver

    # AR SKYBRICKS_DIR: for fiberassign/2.4.0, SKYBRICKS_DIR was not recorded yet in the header
    if mydict["fa_ver"] == "2.4.0":
        fixed_skybrver = "v2"
        log.info(
            "{:.1f}s\t{}\tmodifying mydict[skybrver]={} to mydict[skybrver]={}".format(
                time() - start, step, mydict["skybrver"], fixed_skybrver,
            )
        )
        mydict["skybrver"] = fixed_skybrver

    # AR MJD window for ToO
    # AR - fa_ver <= 5.1.1: mjd_min = mjd_now - 1, mjd_max = mjd_now + 30
    if mydict["fa_ver"] <= "5.1.1":
        try:
            mjd_now = Time(
                datetime.strptime(mydict["pmtime_utc_str"], "%Y-%m-%dT%H:%M:%S%z")
            ).mjd
        except ValueError:
            mjd_now = Time(
                datetime.strptime(mydict["pmtime_utc_str"], "%Y-%m-%dT%H:%M:%S")
            ).mjd
        mjd_min, mjd_max = mjd_now - 1, mjd_now + 30
        log.info(
            "{:.1f}s\t{}\tsetting mjd_min={} and mjd_max={} as fa_ver = 5.1.1".format(
                time() - start, step, mjd_min, mjd_max
            )
        )
        mydict["mjd_min"], mydict["mjd_max"] = mjd_min, mjd_max

    # AR ToO were first run for SV3 on 20210419 (Pacific)
    if mydict["pmtime_utc_str"] < "2021-04-19T07:00:00+00":
        log.info(
            "{:.1f}s\t{}\tsetting dotoo=False, because mydict[pmtime_utc_str]={} < 2021-04-19T07:00:00+00".format(
                time() - start, step, mydict["pmtime_utc_str"],
            )
        )
        mydict["dotoo"] = False

    return mydict


def fba_rerun_intermediate(
    fn, outdir, bugfix=True, tmpoutdir=tempfile.mkdtemp(), log=Logger.get(), start=time(),
):
    """
    Re-create the intermediate files (TILEID-{tiles,sky,targ,scnd,too}.fits) necessary
        to rerun the fiber assignment of a fiberassign-TILEID.fits.gz file.

    Args:
        fn: full path to the fiberassign-TILEID.fits file to be rerun (string)
        outdir: output folder (files will be written in outdir/ABC/) (string)
        bugfix (optional, defaults to True): bugfix the arguments with fba_rerun_bugfix_settings()? (boolean)
        tmpoutdir (optional, defaults to a temporary directory): temporary directory where
                write_targets will write (creating some sub-directories)
        log (optional, defaults to Logger.get()): Logger object
        start (optional, defaults to time()): start time for log (in seconds; output of time.time())

    Notes:
        The SUBPRIORITY overwritting for SV3 is handled with the desitarget code version
            used to rerun_intermediate(); it should be <1.2.2; we suggest to use 1.0.0
    """
    #
    start = time()
    log.info("{:.1f}s\tstart\tTIMESTAMP={}".format(time() - start, Time.now().isot))
    log.info("")
    log.info("")
    # AR config infos
    print_config_infos(log=log, step="settings", start=start)
    # AR get settings
    mydict = fba_rerun_get_settings(
        fn, bugfix=bugfix, log=log, step="settings", start=start
    )
    log.info("{:.1f}s\tsettings\tmydict = {}".format(time() - start, mydict))
    # AR setting outdir/ABC + safe check
    outdir = os.path.join(
        os.path.normpath(outdir), "{:06d}".format(mydict["tileid"])[:3]
    )
    if outdir == os.path.dirname(fn):
        log.error(
            "{:.1f}s\tsettings\tnot safe to write in the original folder {}; exiting".format(
                time() - start, outdir
            )
        )
        sys.exit(1)
    # AR output files
    myouts = {}
    keys = ["tiles", "sky", "gfa", "targ", "scnd"]
    if mydict["dotoo"]:
        keys.append("too")
    for key in keys:
        myouts[key] = os.path.join(
            outdir, "{:06d}-{}.fits".format(mydict["tileid"], key)
        )
    # AR get desitarget paths
    mydirs = get_desitarget_paths(
        mydict["dtver"],
        mydict["survey"],
        mydict["program"],
        dr=mydict["dr"],
        gaiadr=mydict["gaiadr"],
        log=log,
        step="settings",
        start=start,
    )
    # AR tiles
    create_tile(
        mydict["tileid"],
        mydict["tilera"],
        mydict["tiledec"],
        myouts["tiles"],
        mydict["survey"],
        mydict["obscon"],
        log=log,
        step="dotile",
        start=start,
    )
    # AR sky
    create_sky(
        myouts["tiles"],
        mydirs["sky"],
        myouts["sky"],
        suppskydir=mydirs["skysupp"],
        tmpoutdir=tmpoutdir,
        log=log,
        step="dosky",
        start=start,
    )
    # AR targ
    create_mtl(
        myouts["tiles"],
        mydirs["mtl"],
        mydict["mtltime"],
        mydirs["targ"],
        mydict["survey"],
        mydict["gaiadr"].replace("gaia", ""),
        mydict["pmcorr"],
        myouts["targ"],
        pmtime_utc_str=mydict["pmtime_utc_str"],
        tmpoutdir=tmpoutdir,
        log=log,
        step="domtl",
        start=start,
    )
    # AR secondary targets
    targdirs = [mydirs["scnd"]]
    for key in sorted(list(mydirs.keys())):
        if (key[:4] == "scnd") & (key != "scnd") & (key != "scndmtl"):
            targdirs.append(mydirs[key])
    create_mtl(
        myouts["tiles"],
        mydirs["scndmtl"],
        mydict["mtltime"],
        targdirs,
        mydict["survey"],
        mydict["gaiadr"].replace("gaia", ""),
        mydict["pmcorr"],
        myouts["scnd"],
        pmtime_utc_str=mydict["pmtime_utc_str"],
        tmpoutdir=tmpoutdir,
        log=log,
        step="doscnd",
        start=start,
    )
    # AR ToO targets
    if mydict["dotoo"]:
        _ = create_too(
            myouts["tiles"],
            mydirs["too"],
            mydict["mjd_min"],
            mydict["mjd_max"],
            mydict["survey"],
            mydict["gaiadr"].replace("gaia", ""),
            mydict["pmcorr"],
            myouts["too"],
            pmtime_utc_str=mydict["pmtime_utc_str"],
            tmpoutdir=tmpoutdir,
            log=log,
            step="dotoo",
            start=start,
        )
    #
    log.info("")
    log.info("")
    log.info("{:.1f}s\tend\tTIMESTAMP={}".format(time() - start, Time.now().isot))


def get_fba_rerun_scriptname(infiberassignfn, outdir, rerun="fa"):
    """
    Get the fba_rerun_script() *sh scripts filenames.

    Args:
        infiberassignfn: full path to the fiberassign-TILEID.fits.gz file to be rerun (string)
        outdir: output folder (files will be written in outdir/ABC/) (string)
        rerun (optional, defaults to "fa"): type of rerun, "fa" or "intermediate" (string)

    Returns:
        *sh script filename (string)
    """
    if rerun not in ["intermediate", "fa"]:
        sys.exit("rerun={} not in ['intermediate', 'rerun']; exiting".format(rerun))
    tileid = fits.getheader(infiberassignfn, 0)["TILEID"]
    subdir = "{:06d}".format(tileid)[:3]
    if rerun == "fa":
        return os.path.join(outdir, subdir, "rerun-fa-{:06d}.sh".format(tileid))
    else:
        return os.path.join(
            outdir, subdir, "rerun-intermediate-{:06d}.sh".format(tileid)
        )


def fba_rerun_fbascript(
    infiberassignfn,
    outdir,
    run_intermediate,
    bugfix=True,
    outfba_type="none",
    outfiberassign_type="zip",
    faver_noswap=False,
    run_check=True,
    intermediate_dir_final="-",
    overwrite=False,
    log=Logger.get(),
    start=time(),
):
    """
    Write a bash script (outdir/ABC/rerun-fa-TILEID.sh) to rerun fiber assignment of a fiberassign-TILEID.fits.gz file.
    If run_intermediate=True, also write a bash script (outdir/ABC/rerun-intermediate-TILEID.sh) to rerun the intermediate
        products (TILEID-{tiles,sky,targ,scnd,too}.fits).

    Args:
        infiberassignfn: full path to the fiberassign-TILEID.fits.gz file to be rerun (string)
        outdir: output folder (files will be written in outdir/ABC/) (string)
        run_intermediate: generate the intermediate products? (outdir/ABC/TILEID-{tiles,sky,targ,scnd,too}.fits) (boolean)
            if run_intermediate=False, the intermediate products need to be present
        bugfix (optional, defaults to True): bugfix the arguments with fba_rerun_bugfix_settings()? (boolean)
        outfba (optional, defaults to "none"): "none"->no fba-TILEID.fits, "unzip"->fba-TILEID.fits, "zip"->fba-TILEID.fits.gz (string)
        outfiberassign_type (optional, defaults to "zip"): "none"->no fiberassign-TILEID.fits, "unzip"->fiberassign-TILEID.fits, "zip"->fiberassign-TILEID.fits.gz (string)
        faver_noswap (optional, defaults to False): if True, does not "module swap fiberassign/xxx" to the original fiberassign code version  (boolean)
        run_check (optional, defaults to True): run fba_rerun_check()? (boolean)
        intermediate_dir_final (optional, defaults to "-"): folder to move all intermediate products (*sh, *fits, *log) (string)
        overwrite (optional, defaults to False): overwrite any existing file? (boolean)
        log (optional, defaults to Logger.get()): Logger object
        start (optional, defaults to time()): start time for log (in seconds; output of time.time())

    Notes:
        The bash script requires that the desi (master) environment is already loaded,
            i.e. the following has been executed: "source /global/cfs/cdirs/desi/software/desi_environment.sh master".
        The code will use all the intermediate products present in outdir/ABC, and only those.
        The code will create the sub-folders ABC/ if they do not exist.
        If any of the following (according to arguments) files already exist:
            - outdir/ABC/rerun-fa-TILEID.{sh, log, diff} , outdir/ABC/{fba,fiberassign}-TILEID.fits{.gz},
            - intermediate_dir_final/ABC/rerun-intermediate-TILEID.{sh,log} , TILEID-{tiles,sky,targ,scnd,too}.fits
        the code will:
            - exit with an error if overwrite=False;
            - beforehand delete those if overwrite=True.
    """
    # AR assess arguments
    if not os.path.isfile(infiberassignfn):
        log.error("no {} file; exiting".format(infiberassignfn))
        sys.exit(1)
    if not os.path.isdir(outdir):
        log.error("no {} folder; exiting".format(outdir))
        sys.exit(1)
    if run_intermediate not in [True, False]:
        log.error("run_intermediate={} not a boolean; exiting".format(run_intermediate))
        sys.exit(1)
    if outfba_type not in ["none", "unzip", "zip"]:
        log.error("fba={} not in 'none', 'unzip', zip'; exiting".format(outfba_type))
        sys.exit(1)
    if outfiberassign_type not in ["none", "unzip", "zip"]:
        log.error("fiberassign={} not in 'none', 'unzip', zip'; exiting").format(
            outfiberassign_type
        )
        sys.exit(1)
    if run_check not in [True, False]:
        log.error("run_check={} not a boolean; exiting".format(run_check))
        sys.exit(1)
    if (outfiberassign_type == "none") & (run_check):
        log.error(
            "outfiberassign_type=none and run_check=True: run_check=True requires outfiberassign_type!=none; exiting"
        )
        sys.exit(1)
    if intermediate_dir_final != "-":
        if not os.path.isdir(intermediate_dir_final):
            log.error("no {} folder; exiting".format(intermediate_dir_final))
            sys.exit(1)

    # AR get settings, fiberassign version, and SKYBRICKS_DIR and SKYHEALPIXS version
    mydict = fba_rerun_get_settings(
        infiberassignfn, bugfix=bugfix, log=log, step="settings", start=start
    )

    # AR create sub-folders?
    subdir = "{:06d}".format(mydict["tileid"])[:3]
    mydirs = [os.path.join(outdir, subdir)]
    if intermediate_dir_final != "-":
        mydirs.append(os.path.join(intermediate_dir_final, subdir))
    for mydir in mydirs:
        if not os.path.isdir(mydir):
            log.info("create {}".format(mydir))
            os.system("mkdir {}".format(mydir))

    # AR output files
    outsh = get_fba_rerun_scriptname(infiberassignfn, outdir, rerun="fa")
    outlog = outsh.replace(".sh", ".log")
    outfba = os.path.join(outdir, subdir, "fba-{:06d}.fits".format(mydict["tileid"]))
    outfiberassign = os.path.join(
        outdir, subdir, "fiberassign-{:06d}.fits".format(mydict["tileid"])
    )
    fns = [
        outsh,
        outlog,
        outfba,
        "{}.gz".format(outfba),
        outfiberassign,
        "{}.gz".format(outfiberassign),
    ]
    if run_intermediate:
        outintsh = get_fba_rerun_scriptname(
            infiberassignfn, outdir, rerun="intermediate"
        )
        outintlog = outintsh.replace(".sh", ".log")
        fns += [outintsh, outintlog]
    if run_check:
        outdiff = outsh.replace(".sh", ".diff")
        fns.append(outdiff)
    for fn in fns:
        if os.path.isfile(fn):
            if overwrite:
                log.info("remove already existing {}".format(fn))
                os.remove(fn)
            else:
                log.error("{} already exists; exiting".format(fn))
                sys.exit(1)

    # AR files
    mydict["dir"] = os.path.join(outdir, subdir)
    mydict["footprint"] = os.path.join(
        outdir, subdir, "{:06d}-tiles.fits".format(mydict["tileid"])
    )
    mydict["sky"] = os.path.join(outdir, subdir, "{:06d}-sky.fits".format(mydict["tileid"]))
    mydict["targets"] = os.path.join(outdir, subdir, "{:06d}-targ.fits".format(mydict["tileid"]))
    if not run_intermediate:
        for fn in [mydict["footprint"], mydict["sky"], mydict["targets"]]:
            if not os.path.isfile(fn):
                log.error("{} not present; exiting".format(fn))
                sys.exit(1)
    # AR scnd and too
    scndfn = os.path.join(outdir, subdir, "{:06d}-scnd.fits".format(mydict["tileid"]))
    toofn = os.path.join(outdir, subdir, "{:06d}-too.fits".format(mydict["tileid"]))

    # AR intermediate files to potentially move if intermediate_dir_final
    if intermediate_dir_final != "-":
        int2mv_fns = [outintsh, outintlog]
        int2mv_fns += [mydict["footprint"], mydict["sky"], mydict["targets"]]
        int2mv_fns += [scndfn, toofn]
        for fn in int2mv_fns:
            tmpfn = os.path.join(intermediate_dir_final, subdir, os.path.basename(fn))
            if os.path.isfile(tmpfn):
                if overwrite:
                    log.info("remove already existing {}".format(tmpfn))
                    os.remove(tmpfn)
                else:
                    log.error("{} already exists; exiting".format(tmpfn))
                    sys.exit(1)

    # AR run intermediate products?
    if run_intermediate:
        f = open(outintsh, "w")
        f.write("#!/bin/bash\n")
        f.write("\n")
        if mydict["survey"] == "sv3":
            f.write("module swap desitarget/1.0.0\n")
            f.write("echo 'module swap desitarget/1.0.0' >> {}\n".format(outintlog))
        f.write("module list 2> {}\n".format(outintlog))
        f.write("\n")
        f.write(
            'python -c \'from fiberassign.fba_rerun_io import fba_rerun_intermediate; fba_rerun_intermediate("{}", "{}")\' >> {} 2>&1\n'.format(
                infiberassignfn, outdir, outintlog,
            )
        )
        # AR close + make it executable
        f.close()
        os.system("chmod +x {}".format(outintsh))

    # AR writing outsh file
    f = open(outsh, "w")
    f.write("#!/bin/bash\n")
    f.write("\n")
    f.write("# ===============\n")
    f.write("# setting + printing the environment\n")
    f.write("# ===============\n")
    if not faver_noswap:
        f.write("module swap fiberassign/{}\n".format(mydict["fa_ver"]))
        f.write("echo 'module swap fiberassign/{}' >> {}\n".format(mydict["fa_ver"], outlog))
    else:
        f.write("echo 'faver_noswap=True, so we do not swap to the original fiberassign code version' >> {}\n".format(outlog))
    f.write("module list 2>> {}\n".format(outlog))
    if mydict["skybrver"] == "-":
        f.write("unset SKYBRICKS_DIR\n")
    else:
        f.write(
            "export SKYBRICKS_DIR=$DESI_ROOT/target/skybricks/{}\n".format(mydict["skybrver"])
        )
    # AR similarly adding SKYHEALPIXS
    if mydict["skyhpver"] == "-":
        f.write("unset SKYHEALPIXS_DIR\n")
    else:
        f.write(
            "export SKYHEALPIXS_DIR=$DESI_ROOT/target/skyhealpixs/{}\n".format(mydict["skyhpver"])
        )


    # AR control informations
    f.write('echo "fba_run executable: `which fba_run`" >> {}\n'.format(outlog))
    f.write(
        "echo \"fiberassign path: `python -c 'import fiberassign; print(fiberassign.__path__)'`\" >> {}\n".format(
            outlog
        )
    )
    f.write(
        "echo \"fiberassign version: `python -c 'import fiberassign; print(fiberassign.__version__)'`\" >> {}\n".format(
            outlog
        )
    )
    f.write('echo "SKYBRICKS_DIR=$SKYBRICKS_DIR" >> {}\n'.format(outlog))
    f.write('echo "SKYHEALPIXS_DIR=$SKYHEALPIXS_DIR" >> {}\n'.format(outlog))
    f.write("\n")
    f.write("\n")

    # AR constructing the fba_run call
    # AR purposely keep --targets at the end,
    # AR so that we can complete with scnd and too, if they exist
    # AR because we do not know if they will exist now
    #
    # AR fba_run arguments: --targets needs to be the last one!
    fbarun_keys = ["rundate", "fieldrot", "standards_per_petal", "sky_per_petal", "dir", "footprint", "sky"]
    if mydict["fa_ver"] >= "2.4.0":
        fbarun_keys += ["sky_per_slitblock"]
    if mydict["fa_ver"] >= "3.0.0":
        fbarun_keys += ["ha"]
        fbarun_keys += ["margin-pos", "margin-petal", "margin-gfa"]
    fbarun_keys += ["targets"]
    fbarun_cmd = "fba_run --write_all_targets"
    for key in fbarun_keys:
        fbarun_cmd += " --{} {}".format(key, mydict[key])
    f.write("# ===============\n")
    f.write("# fba_run call\n")
    f.write("# ===============\n")
    f.write('CMD="{}"\n'.format(fbarun_cmd))
    f.write("for FN in {} {}\n".format(scndfn, toofn))
    f.write("do\n")
    f.write("\tif [[ -f $FN ]]\n")
    f.write("\tthen\n")
    f.write('\t\tCMD=`echo $CMD " " $FN`\n')
    f.write("\tfi\n")
    f.write("done\n")
    f.write("echo $CMD >> {}\n".format(outlog))
    f.write("eval $CMD >> {} 2>&1\n".format(outlog))
    f.write("\n")
    f.write("\n")

    # AR constructing the merge_results + gzipping, if requested
    # AR need to use a non-straightforward approach
    # AR (with lots of quotes in quotes...)
    # AR to handle the possible existence of scnd and too files
    if outfiberassign_type != "none":
        merge_cmd = "CMD=\"python -c 'from fiberassign.assign import merge_results; "
        merge_cmd += 'merge_results([\\"{}\\""'.format(mydict["targets"])
        f.write("# ===============\n")
        f.write("# merge_results call\n")
        f.write("# ===============\n")
        f.write("{}\n".format(merge_cmd))
        f.write("for FN in {} {}\n".format(scndfn, toofn))
        f.write("do\n")
        f.write("\tif [[ -f $FN ]]\n")
        f.write("\tthen\n")
        f.write('\t\tCMD=`echo $CMD ", \\""$FN"\\""`\n')
        f.write("\tfi\n")
        f.write("done\n")
        f.write('CMD=`echo $CMD "], "`\n')
        merge_cmd = '[\\"{}\\"], '.format(mydict["sky"])
        merge_cmd += "[{}], ".format(mydict["tileid"])
        merge_cmd += 'result_dir=\\"{}\\", '.format(mydict["dir"])
        merge_cmd += "columns=None, "
        merge_cmd += "copy_fba=False"
        merge_cmd += ")'"
        f.write('CMD=`echo $CMD "{}"`\n'.format(merge_cmd))
        f.write("echo $CMD >> {}\n".format(outlog))
        f.write("eval $CMD >> {} 2>&1\n".format(outlog))
        f.write("\n")
        f.write("\n")
        if outfiberassign_type == "zip":
            f.write("# ===============\n")
            f.write("# gzipping fiberassign-TILEID.fits\n")
            f.write("# ===============\n")
            f.write("gzip {}\n".format(outfiberassign))
            f.write("\n")
            f.write("\n")
    if outfba_type == "none":
        f.write("# ===============\n")
        f.write("# removing fba-TILEID.fits\n")
        f.write("# ===============\n")
        f.write("rm {}\n".format(outfba))
        f.write("\n")
        f.write("\n")
    if outfba_type == "zip":
        f.write("# ===============\n")
        f.write("# gzipping fba-TILEID.fits\n")
        f.write("# ===============\n")
        f.write("gzip {}\n".format(outfba))
        f.write("\n")
        f.write("\n")

    # AR run check?
    if run_check:
        if outfiberassign_type == "unzip":
            tmpfn = outfiberassign
        if outfiberassign_type == "zip":
            tmpfn = "{}.gz".format(outfiberassign)
        f.write("# ===============\n")
        f.write("# running check\n")
        f.write("# ===============\n")
        # AR swapping back to fiberassign/master
        f.write("module swap fiberassign/master")
        f.write("echo 'module swap fiberassign/master to perform fba_rerun_check()' >> {}\n".format(outlog))
        f.write(
            'python -c \'from fiberassign.fba_rerun_io import fba_rerun_check; fba_rerun_check("{}", "{}", "{}")\' >> {} 2>&1\n'.format(
                infiberassignfn, tmpfn, outdiff, outlog,
            )
        )

    # AR move intermediate products to a final folder?
    if intermediate_dir_final != "-":
        mydir = os.path.join(intermediate_dir_final, subdir)
        f.write("\n")
        f.write("\n")
        f.write("# ===============\n")
        f.write("# moving intermediate files\n")
        f.write("# ===============\n")
        f.write("for FN in {}\n".format(" ".join(int2mv_fns)))
        f.write("do\n")
        f.write('\tif [[ -f "$FN" ]]\n')
        f.write("\tthen\n")
        f.write('\t\techo "moving $FN to {}" >> {}\n'.format(mydir, outlog))
        f.write("\t\tmv $FN {}\n".format(mydir))
        f.write("\telse\n")
        f.write('\t\techo "no $FN to move" >> {}\n'.format(outlog))
        f.write("\tfi\n")
        f.write("done\n")

    # AR close + make it executable
    f.close()
    os.system("chmod +x {}".format(outsh))


def fba_rerun_check(
    origfn, rerunfn, difffn, log=Logger.get(), start=time(),
):
    """
    Compare the {TARGETID}-{FIBER} values for the FIBERASSIGN and POTENTIAL_ASSIGNMENTS extensions,
        and report differences.

    Args:
        origfn: path to the original fiberassign file (string)
        rerunfn: path to the rerun fiberassign file (string)
        difffn: path for the output diagnosis file
        log (optional, defaults to Logger.get()): Logger object
        start (optional, defaults to time()): start time for log (in seconds; output of time.time())

    Notes:
        Will generate difffn in any case; if no differences, then difffn will just have the header line.
    """
    # AR internal function to add information from the TARGETS extension to POTENTIAL_ASSIGNMENTS rows
    # AR assume there are no TARGETID duplicate
    def populate_potass(h, keys):
        #
        targs_tids = h["TARGETS"].data["TARGETID"]
        potass_tids = h["POTENTIAL_ASSIGNMENTS"].data["TARGETID"]
        #
        unq_potass_tids, ii_inv = np.unique(potass_tids, return_inverse=True)
        ii_targs = match_to(targs_tids, unq_potass_tids)
        if ii_targs.size != unq_potass_tids.size:
            log.error(
                "mismatch: ii_targs.size={} and unq_potass_tids.size={}; exiting".format(
                    ii_targs.size, unq_potass_tids.size
                )
            )
            sys.exit(1)
        #
        potass_d = Table()
        for key in keys:
            if key not in ["FIBER", "PRIORITY_INIT"]:
                unq_vals = h["TARGETS"].data[key][ii_targs]
                potass_d[key] = unq_vals[ii_inv]
        return potass_d

    # AR read files
    orig_h = fits.open(origfn)
    orig_h["TARGETS"].columns["RA"].name = "TARGET_RA"
    orig_h["TARGETS"].columns["DEC"].name = "TARGET_DEC"
    rerun_h = fits.open(rerunfn)
    rerun_h["TARGETS"].columns["RA"].name = "TARGET_RA"
    rerun_h["TARGETS"].columns["DEC"].name = "TARGET_DEC"

    # AR TILEID
    tileid = orig_h[0].header["TILEID"]

    # AR keys we report for mismatches
    dtkeys = ["DESI_TARGET", "MWS_TARGET", "BGS_TARGET", "SCND_TARGET"]
    if "SV3_DESI_TARGET" in orig_h["FIBERASSIGN"].columns.names:
        dtkeys = ["SV3_{}".format(key) for key in dtkeys]
    keys = [
        "TARGETID",
        "FIBER",
        "TARGET_RA",
        "TARGET_DEC",
        "FA_TYPE",
        "PRIORITY_INIT",
        "PRIORITY",
    ]
    keys += dtkeys

    # AR compare
    f = open(difffn, "w")
    f.write("# TILEID\tEXTENSION\tSTATUS\t{}\n".format("\t".join(keys)))
    for ext in ["FIBERASSIGN", "POTENTIAL_ASSIGNMENTS"]:
        # AR FIBERASSIGN
        if ext == "FIBERASSIGN":
            orig_d = orig_h["FIBERASSIGN"].data
            rerun_d = rerun_h["FIBERASSIGN"].data
        # AR POTENTIAL_ASSIGNMENTS
        # AR dummy entry for PRIORITY_INIT
        else:
            orig_d = populate_potass(orig_h, keys)
            orig_d["FIBER"] = orig_h["POTENTIAL_ASSIGNMENTS"].data["FIBER"]
            orig_d["PRIORITY_INIT"] = -99 + np.zeros(len(orig_d), dtype=int)
            rerun_d = populate_potass(rerun_h, keys)
            rerun_d["FIBER"] = rerun_h["POTENTIAL_ASSIGNMENTS"].data["FIBER"]
            rerun_d["PRIORITY_INIT"] = -99 + np.zeros(len(rerun_d), dtype=int)
        # AR {TARGETID}-{FIBER} ids
        orig_ids = np.array(
            [
                "{}-{}".format(tid, fib)
                for tid, fib in zip(orig_d["TARGETID"], orig_d["FIBER"])
            ]
        )
        rerun_ids = np.array(
            [
                "{}-{}".format(tid, fib)
                for tid, fib in zip(rerun_d["TARGETID"], rerun_d["FIBER"])
            ]
        )
        # AR missed in rerun
        ii = np.where(~np.in1d(orig_ids, rerun_ids))[0]
        for i in ii:
            f.write(
                "{:06}\t{}\tmiss\t{}\n".format(
                    tileid,
                    ext,
                    "\t".join(["{}".format(orig_d[key][i]) for key in keys]),
                )
            )
        # AR added in rerun
        ii = np.where(~np.in1d(rerun_ids, orig_ids))[0]
        for i in ii:
            f.write(
                "{:06}\t{}\tadd\t{}\n".format(
                    tileid,
                    ext,
                    "\t".join(["{}".format(rerun_d[key][i]) for key in keys]),
                )
            )
    f.close()
