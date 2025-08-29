# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.utils
=======================

Utility functions.

"""
from __future__ import absolute_import, division, print_function

import os
import subprocess
import sys
from datetime import datetime
from time import time, sleep
import numpy as np
from astropy.table import Table
from astropy.time import Time
import yaml
from pkg_resources import resource_filename
from desiutil.log import get_logger
from ._internal import (Logger, Timer, GlobalTimers, Circle, Segments, Shape,
                        Environment)
log = get_logger()


# Multiprocessing environment setup

default_mp_proc = None
"""Default number of multiprocessing processes.
Set globally on first import.
"""

if "SLURM_CPUS_PER_TASK" in os.environ:
    default_mp_proc = int(os.environ["SLURM_CPUS_PER_TASK"])
else:
    import multiprocessing as _mp
    default_mp_proc = max(1, _mp.cpu_count() // 2)


def option_list(opts):
    """Convert key, value pairs into a list.

    This converts a dictionary into an options list that can be passed to
    ArgumentParser.parse_args().  The value for each dictionary key will be
    converted to a string.  Values that are True will be assumed to not have
    a string argument added to the options list.

    Args:
        opts (dict):  Dictionary of options.

    Returns:
        (list): The list of options.

    """
    optlist = []
    for key, val in opts.items():
        keystr = "--{}".format(key)
        if val is not None:
            if isinstance(val, bool):
                if val:
                    optlist.append(keystr)
            else:
                optlist.append(keystr)
                if isinstance(val, float):
                    optlist.append("{:.14e}".format(val))
                elif isinstance(val, (list, tuple)):
                    optlist.extend(val)
                else:
                    optlist.append("{}".format(val))
    return optlist


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


def get_mjd(utc_time_mjd):
    """
    Get the MJD value for an input date.

    Args:
        utc_time_mjd: date in the "YYYY-MM-DDThh:mm:ss+00:00" format (string)

    Returns:
        MJD value (float)
    """
    assert_isoformat_utc(utc_time_mjd)
    return Time(datetime.strptime(utc_time_mjd, "%Y-%m-%dT%H:%M:%S%z")).mjd


def get_date_cutoff(datetype, cutoff_case):
    """
    Returns the cutoff date for a particular case, from reading the data/cutoff-dates.yaml file.

    Args:
        datetype: one of the key of data/cutoff-dates.yaml (e.g., "rundate", "mtltime", "pmtime_utc_str") (string)
        cutoff_case: one of the config[datetype] keys in data/cutoff-dates.yaml, (e.g. "std_wd" for datetype="rundate") (string)

    Returns:
        date_cutoff: cutoff date for datetype and cutoff_case, in the "YYYY-MM-DDThh:mm:ss+00:00" formatting (string)

    """
    fn = resource_filename("fiberassign", "data/cutoff-dates.yaml")
    with open(fn, "r") as f:
        config = yaml.safe_load(f)
    f.close()
    date_cutoff = config[datetype][cutoff_case]
    return date_cutoff


def get_default_static_obsdate():
    """
    Returns the 'historical' default obsdate value.

    Args:
        None

    Returns:
        "2022-07-01"
    """
    return "2022-07-01"


def get_obsdate(rundate=None):
    """
    Returns the default obsdate: "2022-07-01" if rundate=None or rundate<obsdate_cutoff.

    Args:
        rundate (optional, defaults to None): rundate, in the "YYYY-MM-DDThh:mm:ss+00:00" format (string)

    Returns:
        obsdate: "YYYY-MM-DD" format (string)
        is_after_cutoff: is rundate after rundate_cutoff? (bool)
    """
    obsdate = get_default_static_obsdate()
    is_after_cutoff = False
    if rundate is None:
        log.info(
            "rundate={} -> (obsdate, is_after_cutoff)=({}, {})".format(
                rundate, obsdate, is_after_cutoff
            )
        )
    else:
        if not assert_isoformat_utc(rundate):
            msg = "rundate={} is not yyyy-mm-ddThh:mm:ss+00:00".format(rundate)
            log.info(msg)
            raise ValueError(msg)
        rundate_mjd = Time(datetime.strptime(rundate, "%Y-%m-%dT%H:%M:%S%z")).mjd
        rundate_cutoff = get_date_cutoff("rundate", "obsdate")
        rundate_mjd_cutoff = Time(datetime.strptime(rundate_cutoff, "%Y-%m-%dT%H:%M:%S%z")).mjd
        if rundate_mjd >= rundate_mjd_cutoff:
            is_after_cutoff = True
            yyyy = int(rundate[:4])
            obsdate = "{}{}".format(yyyy + 1, rundate[4:10])
            log.info(
                "rundate={} >= rundate_cutoff={} -> (obsdate, is_after_cutoff)=({}, {})".format(
                    rundate, rundate_cutoff, obsdate, is_after_cutoff)
            )
        else:
            log.info(
                "rundate={} < rundate_cutoff={} -> (obsdate, is_after_cutoff)=({}, {})".format(
                    rundate, rundate_cutoff, obsdate, is_after_cutoff)
            )
    return obsdate, is_after_cutoff


def get_fba_use_fabs(rundate):
    """
    Return the value to which set the $FIBERASSIGN_ALGORITHM_EPOCH environment variable,
        which drives some way the cpp computation is done.

    Args:
        rundate: rundate, in the "YYYY-MM-DDThh:mm:ss+00:00" formatting (string)

    Returns:
        fba_use_fabs: an integer (int)

    Notes:
        See this PR https://github.com/desihub/fiberassign/pull/470.
        In the fiberassign running function, then one will set:
            os.environ["FIBERASSIGN_USE_FABS"] = fba_use_fabs
        So far:
        - fba_use_fabs=0 means we reproduce the buggy behavior;
        - fba_use_fabs=1 means we use the expected behavior.
    """
    # AR get the cutoff dates, and corresponding values
    mydict = get_date_cutoff("rundate", "fba_use_fabs")
    cutoff_rundates = np.array([key for key in mydict])
    cutoff_mjds = np.array([get_mjd(_) for _ in cutoff_rundates])
    values = np.array([mydict[_] for _ in cutoff_rundates])
    log.info(
        "fba_use_fabs cutoff dates: {}".format(
            ", ".join(["{} = {}".format(rundate, value) for rundate, value in zip(cutoff_rundates, values)])
        )
    )
    # AR safe, order by increasing mjd
    ii = cutoff_mjds.argsort()
    cutoff_rundates, cutoff_mjds, values = cutoff_rundates[ii], cutoff_mjds[ii], values[ii]
    # AR now pick the earliest date after the input rundate
    # AR as we have set the first cutoff date as 2019-09-16,
    # AR    there will always have some index returned here
    # AR    though still checking if someone queries with an earlier rundate...
    input_mjd = get_mjd(rundate)
    if input_mjd < cutoff_mjds[0]:
        msg = "rundate={} is earlier than DESI! ({})".format(
            rundate, cutoff_rundates[0]
        )
        log.error(msg)
        raise ValueError(msg)

    i = np.where(cutoff_mjds <= get_mjd(rundate))[0][-1]
    log.info("pick fba_use_fabs = {} for rundate = {}".format(values[i], rundate))

    return values[i]


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
    cmd = ["svn", "info", os.path.expandvars(svn_dir)]
    try:
        svn_ver = (
            subprocess.check_output(cmd, stderr=subprocess.DEVNULL).strip().decode()
        )
        # search for "Last Changed Rev: " line and parse out revision number.  Recent versions
        # of svn have a --show-item argument that does this in a less fragile way,
        # but the svn installed at KPNO is old and doesn't support this option.
        searchstr = 'Last Changed Rev: '
        svn_ver = [line[len(searchstr):] for line in svn_ver.split('\n')
                   if searchstr in line][0]
    except subprocess.CalledProcessError:
        svn_ver = "unknown"

    return svn_ver


def get_last_line(fn):
    """
    Return the last line of a text file.

    Args:
        fn: file name (string)

    Returns:
        last_line: (string)

    Notes:
        Fails if fn has one line only; we do not protect for that case,
            as this function is intended to be used in get_program_latest_timestamp()
            to read *ecsv ledgers, which will always have more than one line,
            and we want the fastest function possible, to use in fiberassign on-the-fly.
        Copied from https://stackoverflow.com/questions/46258499/how-to-read-the-last-line-of-a-file-in-python.
    """
    with open(fn, "rb") as f:
        f.seek(-2, os.SEEK_END)
        while f.read(1) != b"\n":
            f.seek(-2, os.SEEK_CUR)
        last_line = f.readline().decode().strip()
    f.close()
    return last_line


def read_ecsv_keys(fn):
    """
    Returns the column content of an .ecsv file.

    Args:
        fn: filename with an .ecsv format (string)

    Returns:
        keys: list of the column names in fn (list)

    Notes:
        Gets the column names from the first line not starting with "#".
    """
    keys = []
    with open(fn) as f:
        for line in f:
            if line[0] == "#":
                continue
            if len(line.strip()) == 0:
                continue
            keys = line.split()
            break
    f.close()
    return keys


def get_fa_gitversion():
    """
    Returns `git describe --tags --dirty --always`, or "unknown" if not a git repo

    Notes:
        Copied (and slightly edited) from desitarget.io.gitversion()
    """
    #
    origdir = os.getcwd()
    os.chdir(os.path.dirname(__file__))
    try:
        p = subprocess.Popen(["git", "describe", "--tags", "--dirty", "--always"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except EnvironmentError:
        return "unknown"
    #
    os.chdir(origdir)
    out = p.communicate()[0]
    if p.returncode == 0:
        # - avoid py3 bytes and py3 unicode; get native str in both cases
        return str(out.rstrip().decode("ascii"))
    else:
        return "unknown"


def get_revs_dates(svndir, subdirs):
    """
    Returns the list of svn revisions (and dates) for a list of subfolders.

    Args:
        svndir: full path to a svn checkout (e.g. "$DESI_TARGET/fiberassign/tiles/trunk") (str)
        subdirs: comma-separated list of subfolders in svndir (e.g. "000,001,080") (str)

    Returns:
        revs: list of revisions (np array of int)
        dates: list of the dates of the revisions (YYYY-MM-DD) (np array of str)

    Notes:
        This function queries the server ("svn log"), so do not run in parallel!
    """
    revs, dates = [], []
    for subdir in subdirs.split(","):
        start = time()
        ps = subprocess.Popen(
            ["svn", "log", os.path.join(svndir, subdir)], stdout=subprocess.PIPE
        )
        output = subprocess.check_output(
            ["grep", "^r[1-9]"],
            stdin=ps.stdout,
        )
        output = output.strip().decode()
        ps.stdout.close()
        lines = output.split("\n")
        revs += [int(line.split("|")[0][1:]) for line in lines]
        dates += [line.split("|")[2].split()[0] for line in lines]
        dt = time() - start
        log.info("subdir={}\tnrev={}\t(dt={:.1f}s)".format(subdir, len(lines), dt))
        # AR protect from overloading the server
        sleep(0.5)
    revs, ii = np.unique(revs, return_index=True)
    dates = np.array(dates)[ii]
    return revs, dates


def get_rev_fiberassign_changes(svndir, rev, subdirs=None):
    """
    Returns a Table() listing the changed fiberassign files for a given svn revision.

    Args:
        svndir: full path to a svn checkout (e.g. "$DESI_TARGET/fiberassign/tiles/trunk") (str)
        rev: the revision (int)
        subdirs (optional, defaults to None):
                comma-separated list of subfolders in svndir (e.g. "000,001,080") (str)
                if provided, cuts the list of fiberassign files on subdirs

    Returns:
        d: a Table() with: TILEID, FILE, DATE, REVISION

    Notes:
        This function queries the server ("svn log"), so do not run in parallel!
    """
    # AR query svn log
    start = time()
    output = subprocess.check_output(
        ["svn", "log", "-r", str(rev), "-v", svndir], stderr=subprocess.DEVNULL
    )
    output = output.strip().decode()
    dt = time() - start
    # AR first get date
    line_date = [
        line
        for line in output.split("\n")
        if line[:2] in ["r{}".format(i) for i in range(1, 10)]
    ][0]
    date = line_date.split("|")[2].split()[0]
    # AR now cut on lines related to fiberassign changes
    lines = [
        line
        for line in output.split("\n")
        if "/tiles/trunk" in line and "/fiberassign-" in line
    ]
    changes = np.array([line.split()[0] for line in lines])
    fns = np.array([line.split()[1] for line in lines])
    #
    myd = {key: [] for key in ["TILEID", "FILE", "DATE", "REVISION", "CHANGE"]}
    dt = time() - start
    # AR cut on subdirs
    tmp_tileids = np.array([int(os.path.basename(fn)[12:18]) for fn in fns])
    tmp_subdirs = np.array(["{:06d}".format(tileid)[:3] for tileid in tmp_tileids])
    sel = np.in1d(tmp_subdirs, subdirs.split(","))
    fns, changes = fns[sel], changes[sel]
    log.info("rev={}\tdate={}\tnfile={}\t(dt={:.1f}s)".format(rev, date, len(fns), dt))
    # AR store
    # AR (forcing types to protect cases where len(fns)=0)
    d = Table()
    d.meta["SVNDIR"] = svndir
    d.meta["SUBDIRS"] = subdirs
    d["TILEID"] = np.array([int(os.path.basename(fn)[12:18]) for fn in fns], dtype=int)
    d["FILE"] = np.array(fns, dtype=str)
    d["DATE"] = np.array([date for fn in fns], dtype="<U10")
    d["REVISION"] = np.array([rev for fn in fns], dtype=int)
    d["CHANGE"] = np.array(changes, dtype=str)
    return d


def get_obstheta_corr(decs, has, clip_arcsec=600.):
    """
    Returns the computed correction to be applied to the field rotation
        computed in fiberassign.tiles.load_tiles().
    The correction should be applied as: obsthetas[deg] -= obstheta_corrs[deg].

    Args:
        decs: tile declinations (float or np array of floats)
        has: hour angles (float or np array of floats)
        clip_arcsec (optional, defaults to 600): abs(obstheta_corrs) is
            forced to be <obstheta_corrs (float)

    Returns:
        obstheta_corrs: correction (in deg) to be applied (float or np array of floats)

    Notes:
        See DocDB-8931 for details.
        During observations, PlateMaker computes the required field rotation,
            then asks the hexapod to be rotated by
            ROTOFFST = FIELDROT - PM_REQ_FIELDROT,
            where FIELDROT is the field rotation coming from fiberassign.tiles.load_tiles().
        When abs(ROTOFFST)>600 arcsec, the move is denied, and the exposure aborted.
        Correction computed here is a fit to ROTOFFST=f(DEC), plus a fit on the residuals.
    """
    assert clip_arcsec >= 0
    isoneval = isinstance(decs, float)
    if isoneval:
        decs, has = np.atleast_1d(decs), np.atleast_1d(has)
    # AR fitted function, ROTOFFST[arcsec] = f(DEC)
    # AR rescale decs into [0, 1]
    xs = (90. - decs) / 180.
    rotoffsts = 937.60578 - 697.06513 * xs ** -0.18835
    # AR fitted function to the residuals, residuals[arcsec] = f(HA)
    xs = (90. + has) / 180.
    residuals = -113.90162 + 222.18009 * xs
    sel = xs > 0.5
    residuals[sel] = -245.49007 + 485.35700 * xs[sel]
    # AR total correction
    obstheta_corrs = rotoffsts + residuals
    # AR clip
    obstheta_corrs = np.clip(obstheta_corrs, -clip_arcsec, clip_arcsec)
    # AR switch to degrees
    obstheta_corrs /= 3600
    if isoneval:
        decs, has = decs[0], has[0]
        return obstheta_corrs[0]
    else:
        return obstheta_corrs


def get_main_dtver(tileid, svndir=None):
    """
    Retrieve the desitarget catalog version for a main tile.

    Args:
        tileid: tileid (int)
        svndir (optional, defaults to $DESI_TARGET/fiberassign/tiles/trunk): svn folder (str)

    Returns:
        dtver: the desitarget catalog version (str)

    Notes:
        The information is obtained from the FAARGS keyword in the
            zero-th extension of the fiberassign-TILEID.fits.gz file
            in the svn folder.
        As of Aug. 2025, the possible outputs are:
            '-', '1.0.0', '1.1.1', '2.2.0', '3.0.0', '3.2.0'
        The '-' output is returned if no such file exists.

    """
    # AR default svn folder
    if svndir is None:
        svndir = os.path.join(os.getenv("DESI_TARGET"), "fiberassign", "tiles", "trunk")

    tileidpad = "{:06d}".format(tileid)
    fn = os.path.join(svndir, tileidpad[:3], "fiberassign-{}.fits.gz".format(tileidpad))
    if os.path.isfile(fn):
        faargs = fitsio.read_header(fn, 0)["FAARGS"].split()
        return [faargs[i+1] for i in range(len(faargs)-1) if faargs[i] == "--dtver"][0]
    else:
        return "-"
