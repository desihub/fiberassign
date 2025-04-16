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
