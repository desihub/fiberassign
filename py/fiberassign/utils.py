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
