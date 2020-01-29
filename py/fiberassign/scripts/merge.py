# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.scripts.merge
============================

High-level functions for merging target catalogs with output files.

"""
from __future__ import absolute_import, division, print_function

import os
import argparse

import numpy as np

from ..utils import Logger

from ..assign import (merge_results, result_tiles)


def parse_merge(optlist=None):
    """Parse merging options.

    This parses either sys.argv or a list of strings passed in.  If passing
    an option list, you can create that more easily using the
    :func:`option_list` function.

    Args:
        optlist (list, optional): Optional list of arguments to parse instead
            of using sys.argv.

    Returns:
        (namespace):  an ArgumentParser namespace.

    """
    parser = argparse.ArgumentParser()

    parser.add_argument("--targets", type=str, required=True, nargs="+",
                        help="Input file with targets of any type.  This "
                        "argument can be specified multiple times (for "
                        "example if standards / skies / science targets are "
                        "in different files).")

    parser.add_argument("--sky", type=str, required=False, nargs="+",
                        help="Input file with sky or 'bad sky' targets.  "
                        "This option exists in order to treat main-survey"
                        " sky target files as valid for other survey types."
                        "  If you are running a main survey assignment, you"
                        " can just pass the sky file to the --targets list.")

    parser.add_argument("--dir", type=str, required=True, default=None,
                        help="Directory containing fiberassign results.")

    parser.add_argument("--prefix", type=str, required=False,
                        default="fba-",
                        help="Prefix of each file (before the <tile>.fits).")

    parser.add_argument("--split", required=False, default=False,
                        action="store_true",
                        help="Results are in tile prefix directories.")

    parser.add_argument("--out", type=str, required=False, default=None,
                        help="Output directory for the merged files.  Default"
                        " is the directory containing the fiberassign output.")

    parser.add_argument("--out_prefix", type=str, required=False,
                        default="fiberassign-",
                        help="Prefix of each output file.")

    parser.add_argument("--out_split", required=False, default=False,
                        action="store_true",
                        help="Split output into tile prefix directories.")

    parser.add_argument("--columns", type=str, required=False, default=None,
                        help="Override the column names of target data to be "
                        "copied from the target files into the fiber "
                        "assignment files.  This should be a comma-separated "
                        "list.")

    parser.add_argument("--copy_raw", required=False, default=False,
                        action="store_true",
                        help="If true, copy the raw fiberassign HDUs to the "
                             "merged output.")

    args = None
    if optlist is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(optlist)

    if args.sky is None:
        args.sky = list()

    return args


def run_merge_init(args, comm=None):
    """Initialize merging inputs.

    This uses the previously parsed options to load the input files needed.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        (tuple):  The (Tiles, columns) needed to run the merging.

    """
    log = Logger.get()
    # Check directory
    if (comm is None) or (comm.rank == 0):
        if not os.path.isdir(args.dir):
            log.error("Results directory {} does not exist".format(args.dir))
            if comm is not None:
                comm.Abort()

    # Check columns
    columns = None
    if args.columns is not None:
        coltest = args.columns.split(",")
        columns = [x for x in coltest if x != "TARGETID"]

    tiles = None
    if (comm is None) or (comm.rank == 0):
        tiles = result_tiles(dir=args.dir, prefix=args.prefix)

    if comm is not None:
        tiles = comm.bcast(tiles, root=0)

    return (tiles, columns)


def run_merge(args):
    """Run output merging.

    This uses the previously parsed options to read input data and perform
    merging of the input catalogs.  This runs on one node and uses
    multiprocessing.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        None

    """
    tiles, columns = run_merge_init(args)
    merge_results(args.targets, args.sky, tiles, result_dir=args.dir,
                  result_prefix=args.prefix, result_split_dir=args.split,
                  out_dir=args.out, out_prefix=args.out_prefix,
                  out_split_dir=args.out_split, columns=columns,
                  copy_fba=args.copy_raw)
    return


def run_merge_mpi(args, comm):
    """Run output merging.

    This uses the previously parsed options to read input data and perform
    merging of the input catalogs.  This is designed to be run with one MPI
    process per node.  Each MPI process in the input communicator will then
    use multiprocessing for further parallelism.

    Args:
        args (namespace): The parsed arguments.
        comm (MPI.Comm):  The MPI communicator.

    Returns:
        None

    """
    log = Logger.get()
    tiles, columns = run_merge_init(args, comm)
    ptiles = np.array_split(tiles, comm.size)[comm.rank]
    if len(ptiles) > 0:
        log.info("proc {} doing {} tiles".format(comm.rank, len(ptiles)))
        merge_results(args.targets, args.sky, ptiles, result_dir=args.dir,
                      result_prefix=args.prefix, out_dir=args.out,
                      out_prefix=args.out_prefix, columns=columns,
                      copy_fba=args.copy_raw)
    return
