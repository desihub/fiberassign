# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.scripts.qa
==========================

High-level functions for running QA on outputs.

"""
from __future__ import absolute_import, division, print_function

import os
import argparse
from datetime import datetime

from ..hardware import load_hardware

from ..tiles import load_tiles

from ..qa import qa_tiles


def parse_qa(optlist=None):
    """Parse QA options.

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

    parser.add_argument("--dir", type=str, required=True, default=None,
                        help="Directory containing fiberassign outputs.")

    parser.add_argument("--prefix", type=str, required=False,
                        default="fiberassign_",
                        help="Prefix of each file (before the <tile>.fits).")

    parser.add_argument("--split", required=False, default=False,
                        action="store_true",
                        help="Results are in tile prefix directories.")

    parser.add_argument("--out", type=str, required=False, default=None,
                        help="Output file for QA.  Default is a file in"
                        " directory containing the fiberassign output.")

    parser.add_argument("--footprint", type=str, required=False, default=None,
                        help="Optional FITS file defining the footprint.  If"
                        " not specified, the default footprint from desimodel"
                        " is used.")

    parser.add_argument("--tiles", type=str, required=False, default=None,
                        help="Optional text file containing a subset of the"
                        " tile IDs to use in the footprint, one ID per line."
                        " Default uses all tiles in the footprint.")

    parser.add_argument("--rundate", type=str, required=False, default=None,
                        help="Optional date to simulate for this run of "
                        "fiber assignment, used to load the correct "
                        "focalplane properties and state from desimodel.  "
                        "Default uses the current date.  Format is "
                        "YYYY-MM-DDTHH:mm:ss in UTC time.")

    args = None
    if optlist is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(optlist)

    # Check directory
    if not os.path.isdir(args.dir):
        raise RuntimeError("Results directory {} does not exist"
                           .format(args.dir))

    if args.out is None:
        args.out = os.path.join(args.dir, "qa.json")

    return args


def run_qa_init(args):
    """Initialize QA inputs.

    This uses the previously parsed options to load the input files needed.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        (tuple):  The (Hardware, Tiles) needed to run QA.

    """
    # Read hardware properties
    hw = load_hardware(rundate=args.rundate)

    # Read tiles we are using
    tileselect = None
    if args.tiles is not None:
        tileselect = list()
        with open(args.tiles, "r") as f:
            for line in f:
                # Try to convert the first column to an integer.
                try:
                    tileselect.append(int(line.split()[0]))
                except ValueError:
                    pass
    tiles = load_tiles(tiles_file=args.footprint, select=tileselect)

    return (hw, tiles)


def run_qa(args):
    """Run QA.

    This uses the previously parsed options to read input data and run the
    QA.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        None

    """
    hw, tiles = run_qa_init(args)
    qa_tiles(hw, tiles, result_dir=args.dir, result_prefix=args.prefix,
             result_split_dir=args.split, qa_out=args.out)
    return
