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
#import desitarget
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


def fba_rerun_intermediate_get_settings(
    fn, log=Logger.get(), step="rerun", start=time(),
):
    """
    Get the required settings from a fiberassign-TILEID.fits.gz file to recreate the
        intermediate files (TILEID-{tiles,sky,targ,scnd,too}.fits) necessary
        to rerun the fiber assignment.

    Args:
        fn: full path to the fiberassign-TILEID.fits file to be rerun (string0
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()

    Returns:
        mydict: dictionary with the required fba_launch settings,
            stored in the fiberassign*fits.gz header FAARGS (dictionary)

    Notes:
        Only runs for a SV3 or Main BRIGHT/DARK tile; exists with error otherwise.
        Special cases handling:
        - RUNDATE fix for 19 SV3 tiles designed on 2021-04-10
        - MTLTIME fix for SV3 TILEID=315 and 441
        - DTVER: 1.0.0->1.1.1 fix for pre-shutdown Main (NIGHT<=20210518)
        - UTC formatting for early SV3 tiles
        The SUBPRIORITY overwritting for SV3 is handled with the desitarget code version
            used to rerun_intermediate()
    """
    #
    hdr = fits.getheader(fn, 0)
    # AR SV3/Main BRIGHT/DARK?
    faflavors = ["sv3bright", "sv3dark", "mainbright", "maindark"]
    if hdr["FAFLAVOR"] not in faflavors:
        sys.exit("FAFLAVOR={} not in {}; exiting".format(hdr["FAFLAVOR"], faflavors))


    mydict = {}


    # AR fa_ver, obscon: from header
    mydict["fa_ver"] = hdr["FA_VER"]
    mydict["obscon"] = hdr["OBSCON"]


    # AR settings we store
    keys = [
        "survey",
        "program",
        "rundate",
        "mtltime",
        "dr",
        "gaiadr",
        "dtver",
        "pmcorr",
        "pmtime_utc_str",
        "tileid",
        "tiledec",
        "tilera",
    ]
    # AR storing the FAARGS in a dictionary
    faargs = np.array(hdr["FAARGS"].split())
    for key in keys:
        #
        ii = np.where(faargs == "--{}".format(key))[0]
        if len(ii) > 0:
            mydict[key] = faargs[ii[0] + 1]
        # AR args.mtltime has been introduced in 2.4.0
        # AR we add it as an argument if not here
        # AR (meaning that the tile was run with a previous version)
        elif key == "mtltime":
            mydict["mtltime"] = hdr["MTLTIME"]
        # AR args.pmtime_utc_str initially names args.pmtime
        elif (key == "pmtime_utc_str") & ("--pmtime" in faargs):
            ii = np.where(faargs == "--pmtime")[0]
            mydict["pmtime_utc_str"] = faargs[ii[0] + 1]
        # AR should not be other cases
        else:
            log.error(
                "key={} not present in FAARGS, and not in expected missing keys; exiting".format(
                    key
                )
            )
            sys.exit(1)
    #
    mydict["tileid"] = int(mydict["tileid"])
    mydict["tilera"], mydict["tiledec"] = (
        float(mydict["tilera"]),
        float(mydict["tiledec"]),
    )


    # AR case of 19 SV3 tiles designed on 2021-04-10:
    # AR    those have rundate = 2021-04-10T21:28:37
    # AR    there has been a fp update on 2021-04-10T20:00:39+00:00
    # AR    but apparently the desi-state*ecsv file at NERSC was not updated then
    # AR    so the original fiberassign run considered the previous fp state
    # AR    =>
    # AR    we manually modify the rundate, pmtime, mtltime
    if hdr["TILEID"] in [
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
                    time() - start, step, key, mydict[key], fixed_time, hdr["TILEID"],
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
    if hdr["TILEID"] in [315, 441]:
        if hdr["TILEID"] == 315:
            fixed_time = "2021-04-19T21:00:00+00:00"
        if hdr["TILEID"] == 441:
            fixed_time = "2021-04-20T19:00:00+00:00"
        log.info(
            "{:.1f}s\t{}\tmodifying mydict[mtltime] from {} to {} (TILEID={} designed when MTL ledgers not updated at NERSC".format(
                time() - start, step, mydict["mtltime"], fixed_time, hdr["TILEID"],
            )
        )
        mydict["mtltime"] = fixed_time


    # AR assume UTC for pmtime_utc_str, rundate, mtltime
    # AR (UTC formatting introduced in 20210422, during SV3)
    for key in ["pmtime_utc_str", "rundate", "mtltime"]:
        try:
            _ = datetime.strptime(mydict[key], "%Y-%m-%dT%H:%M:%S")
            log.info(
                "{:.1f}s\t{}\t{}: modifying {} to {}+00:00 for TILEID={}".format(
                    time() - start, step, key, mydict[key], mydict[key], hdr["TILEID"]
                )
            )
            mydict[key] += "+00:00"
        except ValueError:
            log.warning(
                "{:.1f}s\t{}\tmydict[{}]={} is not %Y-%m-%dT%H:%M:%S; not changing it, even if args.assume_utc".format(
                    time() - start, step, key, mydict[key],
                )
            )


    # AR dtver: pre-shutdown main, replace 1.0.0 by 1.1.1
    # AR because of missing CALIB targets that were bumped by MWS secondary targets
    if (mydict["survey"] == "main") & (mydict["dtver"] == "1.0.0"):
        log.info(
            "{:.1f}s\tstep\tmodifying mydict[dtver]=1.0.0 to mydict[dtver]=1.1.1".format(
                time() - start, step
            )
        )
        mydict["dtver"] = "1.1.1"


    return mydict




def fba_rerun_intermediate(
    fn, outdir, tmpoutdir=tempfile.mkdtemp(), log=Logger.get(), start=time(),
):
    """
    Re-create the intermediate files (TILEID-{tiles,sky,targ,scnd,too}.fits) necessary
        to rerun the fiber assignment of a fiberassign-TILEID.fits.gz file.

    Args:
        fn: full path to the fiberassign-TILEID.fits file to be rerun (string)
        outdir: output folder (files will be written in outdir/ABC/) (string)
        tmpoutdir (optional, defaults to a temporary directory): temporary directory where
                write_targets will write (creating some sub-directories)
        log (optional, defaults to Logger.get()): Logger object
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()

    Notes:
        No fiberassign here; but FA_VER is used to define the MJD window for the ToO,
            which has been redefined during 2021 shutdown.
        The SUBPRIORITY overwritting for SV3 is handled with the desitarget code version
            used to rerun_intermediate(); it should be <1.2.2; we suggest to use 1.0.0
        ToO SV3: tiles run before 20210419 (Pacific) did not have ToO design.
    """
    #
    start = time()
    log.info("{:.1f}s\tstart\tTIMESTAMP={}".format(time() - start, Time.now().isot))
    log.info("")
    log.info("")
    # AR config infos
    print_config_infos(log=log, step="settings", start=start)
    # AR get settings
    mydict = fba_rerun_intermediate_get_settings(
        fn, log=log, step="settings", start=start
    )
    print(mydict)
    # AR setting outdir/ABC + safe check
    outdir = os.path.join(
        os.path.normpath(outdir), "{:06d}".format(mydict["tileid"])[:3]
    )
    if outdir == os.path.dirname(fn):
        log.error("not safe to write in the original folder {}; exiting".format(outdir))
        sys.exit(1)
    # AR MJD window for ToO
    # AR - if FA_VER <= 5.1.1: mjd_min = mjd_now - 1, mjd_max = mjd_now + 30
    # AR - if FA_VER > 5.1.1: mjd_min = mjd_max = mjd_now
    # AR in any case, mjd_now is defined as pmtime_utc_str
    mjd_now = Time(datetime.strptime(mydict["pmtime_utc_str"], "%Y-%m-%dT%H:%M:%S%z")).mjd
    if mydict["fa_ver"] <= "5.1.1":
        mjd_min, mjd_max = mjd_now - 1, mjd_now + 30
    else:
        mjd_min, mjd_max = mjd_now, mjd_now
    # AR ToO were first run for SV3 on 20210419 (Pacific)
    if mydict["pmtime_utc_str"] < "2021-04-19T07:00:00+00":
        dotoo = False
    else:
        dotoo = True
    # AR output files
    myouts = {}
    keys = ["tiles", "sky", "gfa", "targ", "scnd"]
    if dotoo:
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
    create_mtl(
        myouts["tiles"],
        mydirs["scndmtl"],
        mydict["mtltime"],
        mydirs["scnd"],
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
    if dotoo:
        _ = create_too(
            myouts["tiles"],
            mydirs["too"],
            mjd_min,
            mjd_max,
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


def fba_rerun_get_settings(
    fn, log=Logger.get(), step="rerun", start=time(),
):
    """
    Get the required fiberassign settings from a fiberassign-TILEID.fits.gz file to execute
        fba_run, to rerun the fiber assignment.

    Args:
        fn: full path to the fiberassign-TILEID.fits file to be rerun (string)
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time())

    Returns:
        tileid: tileid (int)
        mydict: dictionary with "rundate", "tileid", "tiledec", "tilera", "fieldrot", "ha" (dictionary)
        survey: survey (string)
        faver: fiberassign code version to use (string)
        skybrver: SKYBRICKS_DIR version ("-" if not set) (string)

    Notes:
        Special cases handling (see code for details):
        - fiberassign/2.2.0.dev2811 and RUNDATE <  2021-04-14 -> fiberassign/2.2.0
        - fiberassign/2.2.0.dev2811 and RUNDATE >= 2021-04-14 -> fiberassign/2.3.0
        - fiberassign2.3.0.dev2838 -> fiberassign/2.3.0
        - fiberassign/2.4.0: SKYBRICKS_DIR was not recorded yet in the header -> we set it to v2
        - RUNDATE fix for 19 SV3 tiles designed on 2021-04-10
    """
    #
    hdr = fits.getheader(fn, 0)
    # AR SV3/Main BRIGHT/DARK?
    faflavors = ["sv3bright", "sv3dark", "mainbright", "maindark"]
    if hdr["FAFLAVOR"] not in faflavors:
        log.error("FAFLAVOR={} not in {}; exiting".format(hdr["FAFLAVOR"], faflavors))
        sys.exit(1)

    # AR survey
    survey = hdr["SURVEY"]

    # AR fiberassign code version
    # AR SV3 handling, as SV3 was run with fiberassign/master for RUNDATE<2021-04-30
    # AR on 2021-04-13: https://github.com/desihub/fiberassign/pull/321
    # AR                this PR changes the assignment
    # AR                so we use 2.2.0 before, 2.3.0 after
    if (hdr["FA_VER"] == "2.2.0.dev2811") & (hdr["RUNDATE"] < "2021-04-14"):
        faver = "2.2.0"
    elif (hdr["FA_VER"] == "2.2.0.dev2811") & (hdr["RUNDATE"] >= "2021-04-14"):
        faver = "2.3.0"
    elif hdr["FA_VER"] == "2.3.0.dev2838":
        faver = "2.3.0"
    else:
        faver = hdr["FA_VER"]


    # AR SKYBRICKS_DIR
    skybrver = "-"  # default value if not set
    keys = [cards[0] for cards in hdr.cards]
    vals = [
        hdr[key.replace("NAM", "VER")].split("/")[-1]
        for key in keys
        if hdr[key] == "SKYBRICKS_DIR"
    ]
    if len(vals) > 0:
        skybrver = vals[0]
    # AR for fiberassign/2.4.0, SKYBRICKS_DIR was not recorded yet in the header
    if faver == "2.4.0":
        skybrver = "v2"


    mydict = {}


    # AR settings from the header that we store
    keys = [
        "rundate",
        "fieldrot",
    ]
    if faver >= "2.4.0":
        keys += ["sky_per_slitblock"]
    if faver >= "3.0.0":
        keys += ["ha"]
    if (faver >= "3.0.0") & (faver < "5.0.0"):
        keys += ["margin-pos", "margin-petal", "margin-gfa"]
    if faver >= "5.0.0":
        keys += ["margin_pos", "margin_petal", "margin_gfa"]


    # AR storing the FAARGS in a dictionary
    faargs = np.array(hdr["FAARGS"].split())
    for key in keys:
        #
        if key[:7] == "margin-":
            ii = np.where(faargs == "--{}".format(key.replace("margin-", "margin_")))[0]
        else:
            ii = np.where(faargs == "--{}".format(key))[0]
        if len(ii) > 0:
            mydict[key] = faargs[ii[0] + 1]
        # AR fieldrot: grab it from the header
        elif key == "fieldrot":
            mydict["fieldrot"] = hdr["FIELDROT"]
        # AR HA: if not provided, grab it from the header (though it should be 0)
        elif key == "ha":
            mydict["ha"] = hdr["FA_HA"]
        # AR sky_per_slitblock, margin_*: if not provided, set to 0
        elif key in ["sky_per_slitblock", "margin_pos", "margin_petal", "margin_gfa"]:
            mydict[key] = 0
        # AR should not be other cases
        else:
            log.error(
                "key={} not present in FAARGS, and not in expected missing keys; exiting".format(
                    key
                )
            )
            sys.exit(1)


    # AR case of 19 SV3 tiles designed on 2021-04-10:
    # AR    those have rundate = 2021-04-10T21:28:37
    # AR    there has been a fp update on 2021-04-10T20:00:39+00:00
    # AR    but apparently the desi-state*ecsv file at NERSC was not updated then
    # AR    so the original fiberassign run considered the previous fp state
    # AR    =>
    # AR    we manually modify the rundate, pmtime, mtltime
    if hdr["TILEID"] in [
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
        for key in ["rundate"]:
            log.info(
                "{:.1f}s\t{}\tmodifying mydict[{}] from {} to {} (TILEID={} designed when desi-state*ecsv file not updated at NERSC".format(
                    time() - start, step, key, mydict[key], fixed_time, hdr["TILEID"],
                )
            )
            mydict[key] = fixed_time


    return hdr["TILEID"], mydict, survey, faver, skybrver


def get_fba_rerun_scriptname(infiberassignfn, outdir, rerun="fa"):
    """
    Gets the fba_rerun_script() *sh scripts filenames.

    Args:
        infiberassignfn: full path to the fiberassign-TILEID.fits file to be rerun (string)
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
        return os.path.join(outdir, subdir, "rerun-intermediate-{:06d}.sh".format(tileid))


def fba_rerun_fbascript(
    infiberassignfn,
    outdir,
    run_intermediate,
    fba="none",
    fiberassign="zip",
    run_check=True,
    intermediate_dir_final="-",
    overwrite=False,
    log=Logger.get(),
    start=time(),
):
    """
    Writes a bash script (outdir/ABC/rerun-fa-TILEID.sh) to rerun fiber assignment of a fiberassign-TILEID.fits.gz file.
    If run_intermediate=True, also writes a bash script (outdir/ABC/rerun-intermediate-TILEID.sh) to rerun the intermediate
        products (TILEID-{tiles,sky,targ,scnd,too}.fits).

    Args:
        infiberassignfn: full path to the fiberassign-TILEID.fits file to be rerun (string)
        outdir: output folder (files will be written in outdir/ABC/) (string)
        run_intermediate: generate the intermediate products? (outdir/ABC/TILEID-{tiles,sky,targ,scnd,too}.fits) (boolean)
            if run_intermediate=False, the intermediate products need to be present
        fba (optional, defaults to "none"): "none"->no fba-TILEID.fits, "unzip"->fba-TILEID.fits, "zip"->fba-TILEID.fits.gz (string)
        fiberassign(optional, defaults to "zip"): "none"->no fiberassign-TILEID.fits, "unzip"->fiberassign-TILEID.fits, "zip"->fiberassign-TILEID.fits.gz (string)
        run_check (optional, defaults to True): run fba_rerun_check()? (boolean)
        intermediate_dir_final (optional, defaults to "-"): folder to move all intermediate products (*sh, *fits, *log) (string)
        overwrite (optional, defaults to False): overwrite any existing file? (boolean)
        log (optional, defaults to Logger.get()): Logger object
        start(optional, defaults to time()): start time for log (in seconds; output of time.time())

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
    if fba not in ["none", "unzip", "zip"]:
        log.error("fba={} not in 'none', 'unzip', zip'; exiting".format(fba))
        sys.exit(1)
    if fiberassign not in ["none", "unzip", "zip"]:
        log.error("fiberassign={} not in 'none', 'unzip', zip'; exiting").format(fiberassign)
        sys.exit(1)
    if run_check not in [True, False]:
        log.error("run_check={} not a boolean; exiting".format(run_check))
        sys.exit(1)
    if (fiberassign == "none") & (run_check):
        log.error("fiberassign=none and run_check=True: run_check=True requires fiberassign!=none; exiting")
        sys.exit(1)
    if intermediate_dir_final != "-":
        if not os.path.isdir(intermediate_dir_final):
            log.error("no {} folder; exiting".format(intermediate_dir_final))
            sys.exit(1)


    # AR get settings, fiberassign version, and SKYBRICKS_DIR version
    tileid, mydict, survey, faver, skybrver = fba_rerun_get_settings(
        infiberassignfn, log=log, step="settings", start=start
    )


    # AR create sub-folders?
    subdir = "{:06d}".format(tileid)[:3]
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
    outfba = os.path.join(outdir, subdir, "fba-{:06d}.fits".format(tileid))
    outfiberassign = os.path.join(
        outdir, subdir, "fiberassign-{:06d}.fits".format(tileid)
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
        outintsh = get_fba_rerun_scriptname(infiberassignfn, outdir, rerun="intermediate")
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
    mydict["footprint"] = os.path.join(outdir, subdir, "{:06d}-tiles.fits".format(tileid))
    mydict["sky"] = os.path.join(outdir, subdir, "{:06d}-sky.fits".format(tileid))
    mydict["targets"] = os.path.join(outdir, subdir, "{:06d}-targ.fits".format(tileid))
    if not run_intermediate:
        for fn in [mydict["footprint"], mydict["sky"], mydict["targets"]]:
            if not os.path.isfile(fn):
                log.error("{} not present; exiting".format(fn))
                sys.exit(1)
    # AR scnd and too
    scndfn = os.path.join(outdir, subdir, "{:06d}-scnd.fits".format(tileid))
    toofn = os.path.join(outdir, subdir, "{:06d}-too.fits".format(tileid))


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
        if survey == "sv3":
            f.write("module swap desitarget/1.0.0\n")
            f.write("echo 'module swap desitarget/1.0.0' >> {}\n".format(outintlog))
        f.write("module list 2> {}\n".format(outintlog))
        f.write("\n")
        f.write(
            "python -c 'from fiberassign.fba_rerun_io import fba_rerun_intermediate; fba_rerun_intermediate(\"{}\", \"{}\")' >> {} 2>&1\n".format(
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
    f.write("module swap fiberassign/{}\n".format(faver))
    f.write("echo 'module swap fiberassign/{}' >> {}\n".format(faver, outlog))
    f.write("module list 2>> {}\n".format(outlog))
    f.write("\n")
    if skybrver == "-":
        f.write("unset SKYBRICKS_DIR\n")
    else:
        f.write(
            "export SKYBRICKS_DIR=$DESI_ROOT/target/skybricks/{}\n".format(skybrver)
        )
    f.write("\n")


    # AR control informations
    f.write("echo \"fba_run executable: `which fba_run`\" >> {}\n".format(outlog))
    f.write("echo \"fiberassign path: `python -c 'import fiberassign; print(fiberassign.__path__)'`\" >> {}\n".format(outlog))
    f.write("echo \"fiberassign version: `python -c 'import fiberassign; print(fiberassign.__version__)'`\" >> {}\n".format(outlog))
    f.write("echo \"SKYBRICKS_DIR=$SKYBRICKS_DIR\" >> {}\n".format(outlog))
    f.write("\n")


    # AR constructing the fba_run call
    # AR purposely keep --targets at the end,
    # AR so that we can complete with scnd and too, if they exist
    # AR because we do not if they will exist now
    fbarun_cmd = "fba_run --write_all_targets"
    for key in list(mydict.keys()):
        if key != "targets":
            fbarun_cmd += " --{} {}".format(key, mydict[key])
    key = "targets"
    fbarun_cmd += " --{} {}".format(key, mydict[key])
    f.write("CMD=\"{}\"\n".format(fbarun_cmd))
    f.write("for FN in {} {}\n".format(scndfn, toofn))
    f.write("do\n")
    f.write("\tif [[ -f $FN ]]\n")
    f.write("\tthen\n")
    f.write("\t\tCMD=`echo $CMD \" \" $FN`\n")
    f.write("\tfi\n")
    f.write("done\n")
    f.write("echo $CMD >> {}\n".format(outlog))
    f.write("eval $CMD >> {} 2>&1\n".format(outlog))
    f.write("\n")


    # AR constructing the merge_results + gzipping, if requested
    # AR need to use a non-straightforward approach
    # AR (with lots of quotes in quotes...)
    # AR to handle the possible existence of scnd and too files
    if fiberassign != "none":
        merge_cmd = "CMD=\"python -c 'from fiberassign.assign import merge_results; "
        merge_cmd += "merge_results([\\\"{}\\\"\"".format(mydict["targets"])
        f.write("{}\n".format(merge_cmd))
        f.write("for FN in {} {}\n".format(scndfn, toofn))
        f.write("do\n")
        f.write("\tif [[ -f $FN ]]\n")
        f.write("\tthen\n")
        f.write("\t\tCMD=`echo $CMD \", \\\"\"$FN\"\\\"\"`\n")
        f.write("\tfi\n")
        f.write("done\n")
        f.write("CMD=`echo $CMD \"], \"`\n")
        merge_cmd = '[\\\"{}\\\"], '.format(mydict["sky"])
        merge_cmd += "[{}], ".format(tileid)
        merge_cmd += 'result_dir=\\\"{}\\\", '.format(mydict["dir"])
        merge_cmd += "columns=None, "
        merge_cmd += "copy_fba=False"
        merge_cmd += ")'"
        f.write("CMD=`echo $CMD \"{}\"`\n".format(merge_cmd))
        f.write("echo $CMD >> {}\n".format(outlog))
        f.write("eval $CMD >> {} 2>&1\n".format(outlog))
        f.write("\n")
        if fiberassign == "zip":
            f.write("gzip {}\n".format(outfiberassign))
            f.write("\n")
    if fba == "none":
        f.write("rm {}\n".format(outfba))
        f.write("\n")
    if fba == "zip":
        f.write("gzip {}\n".format(outfba))
        f.write("\n")

    # AR run check?
    if run_check:
        if fiberassign == "unzip":
            tmpfn = outfiberassign
        if fiberassign == "zip":
            tmpfn = "{}.gz".format(outfiberassign)
        # HACK for development
        f.write("export PYTHONPATH=/global/homes/r/raichoor/software_dev/fiberassign_fba_rerun4/py:$PYTHONPATH\n")
        # HACK
        f.write(
            "python -c 'from fiberassign.fba_rerun_io import fba_rerun_check; fba_rerun_check(\"{}\", \"{}\", \"{}\")' >> {} 2>&1\n".format(
                infiberassignfn, tmpfn, outdiff, outlog,
            )
        )


    # AR move intermediate products to a final folder?
    if intermediate_dir_final != "-":
        mydir = os.path.join(intermediate_dir_final, subdir)
        f.write("\n")
        f.write("for FN in {}\n".format(" ".join(int2mv_fns)))
        f.write("do\n")
        f.write("\tif [[ -f \"$FN\" ]]\n")
        f.write("\tthen\n")
        f.write("\t\techo \"moving $FN to {}\" >> {}\n".format(mydir, outlog))
        f.write("\t\tmv $FN {}\n".format(mydir))
        f.write("\telse\n")
        f.write("\t\techo \"no $FN to move\" >> {}\n".format(outlog))
        f.write("\tfi\n")
        f.write("done\n")


    # AR close + make it executable
    f.close()
    os.system("chmod +x {}".format(outsh))


def fba_rerun_check(
    origfn,
    rerunfn,
    difffn,
    log=Logger.get(),
    start=time(),
):
    """
    Compares the {TARGETID-FIBER} values for the FIBERASSIGN and POTENTIAL_ASSIGNMENTS extensions.

    Args:
        origfn: path to the original fiberassign file (string)
        rerunfn: path to the rerun fiberassign file (string)
        difffn: path for the output diagnosis file
        log (optional, defaults to Logger.get()): Logger object
        start(optional, defaults to time()): start time for log (in seconds; output of time.time())
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
            log.error("mismatch: ii_targs.size={} and unq_potass_tids.size={}; exiting".format(ii_targs.size, unq_potass_tids.size))
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
    keys = ["TARGETID", "FIBER", "TARGET_RA", "TARGET_DEC", "FA_TYPE", "PRIORITY_INIT", "PRIORITY"]
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
        # AR {TARGETID-FIBER} ids
        orig_ids = np.array(["{}-{}".format(tid, fib) for tid, fib in zip(orig_d["TARGETID"], orig_d["FIBER"])])
        rerun_ids = np.array(["{}-{}".format(tid, fib) for tid, fib in zip(rerun_d["TARGETID"], rerun_d["FIBER"])])
        # AR miss
        ii = np.where(~np.in1d(orig_ids, rerun_ids))[0]
        for i in ii:
            f.write("{:06}\t{}\tmiss\t{}\n".format(tileid, ext, "\t".join(["{}".format(orig_d[key][i]) for key in keys])))
        # AR added
        ii = np.where(~np.in1d(rerun_ids, orig_ids))[0]
        for i in ii:
            f.write("{:06}\t{}\tadd\t{}\n".format(tileid, ext, "\t".join(["{}".format(rerun_d[key][i]) for key in keys])))
    f.close()
