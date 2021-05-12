# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.scripts.plot
===========================

High-level functions for plotting output files.

"""
from __future__ import absolute_import, division, print_function

import os
import argparse
from datetime import datetime

from ..hardware import load_hardware

from ..tiles import load_tiles

from ..vis import plot_tiles


def parse_plot(optlist=None):
    """Parse plotting options.

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
                        help="Directory containing fiberassign results.")

    parser.add_argument("--prefix", type=str, required=False,
                        default="fba-",
                        help="Prefix of each file (before the <tile>.fits).")

    parser.add_argument("--split", required=False, default=False,
                        action="store_true",
                        help="Results are in tile prefix directories.")

    parser.add_argument("--out", type=str, required=False, default=None,
                        help="Output directory for plots.  Default is the"
                        " directory containing the fiberassign output.")

    parser.add_argument("--out_prefix", type=str, required=False,
                        default=None,
                        help="Prefix of each output file.")

    parser.add_argument("--out_split", required=False, default=False,
                        action="store_true",
                        help="Split output into tile prefix directories.")

    parser.add_argument("--petals", type=str, required=False, default=None,
                        help="Comma-separated list of petals to plot "
                        "(default is all petals)")

    parser.add_argument("--real_shapes", required=False, default=False,
                        action="store_true",
                        help="Plot the actual positioner shapes.  This looks"
                        " better but takes much longer and makes bigger files."
                        "  Recommended only for plotting limited "
                        "tiles / petals.")

    parser.add_argument("--margin-pos", type=float, required=False, default=0.,
                        help="Add margin (in mm) around positioner keep-out polygons")
    parser.add_argument("--margin-petal", type=float, required=False, default=0.,
                        help="Add margin (in mm) around petal-boundary keep-out polygons")
    parser.add_argument("--margin-gfa", type=float, required=False, default=0.,
                        help="Add margin (in mm) around GFA keep-out polygons")

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
                        "YYYY-MM-DDTHH:mm:ss+-zz:zz.")

    parser.add_argument("--serial", required=False, default=False,
                        action="store_true",
                        help="Disable the use of multiprocessing.  Needed by "
                        "some unit tests to avoid issues with matplotlib.")

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
        args.out = args.dir

    if args.out_prefix is None:
        args.out_prefix = args.prefix

    # Set up margins dict
    args.margins = {}
    if args.margin_pos != 0:
        args.margins['theta'] = args.margin_pos
        args.margins['phi']   = args.margin_pos
    if args.margin_petal != 0:
        args.margins['petal'] = args.margin_petal
    if args.margin_gfa != 0:
        args.margins['gfa'] = args.margin_gfa

    return args


def run_plot_init(args):
    """Initialize plotting inputs.

    This uses the previously parsed options to load the input files needed.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        (tuple):  The (Hardware, Tiles, petals) needed to run plotting.

    """
    # Read hardware properties
    hw = load_hardware(rundate=args.rundate, add_margins=args.margins)

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

    petals = None
    if args.petals is not None:
        petals = [int(x) for x in args.petals.split(",")]

    return (hw, tiles, petals)


def run_plot(args):
    """Run plotting of fiberassign outputs.

    This uses the previously parsed options to read input data and make
    per-tile plots of the fiber assignment.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        None

    """
    hw, tiles, petals = run_plot_init(args)
    plot_tiles(hw, tiles, result_dir=args.dir, result_prefix=args.prefix,
               result_split_dir=args.split, plot_dir=args.out,
               plot_prefix=args.out_prefix, plot_split_dir=args.out_split,
               petals=petals, real_shapes=args.real_shapes)
    return
