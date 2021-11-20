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

from ._internal import (Logger, Timer, GlobalTimers, Circle, Segments, Shape,
                        Environment)

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
