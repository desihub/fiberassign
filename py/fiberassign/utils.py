# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.utils
=======================

Utility functions.

"""
from __future__ import absolute_import, division, print_function

import os

from ._internal import (Logger, Timer, GlobalTimers, Circle, Segments, Shape)

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
