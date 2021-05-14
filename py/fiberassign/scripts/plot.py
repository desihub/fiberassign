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

    parser.add_argument("--out_dir", type=str, required=False, default=None,
                        help="Output directory for all plots.")

    parser.add_argument("--petals", type=str, required=False, default=None,
                        help="Comma-separated list of petals to plot "
                        "(default is all petals)")

    parser.add_argument("--simple", required=False, default=False,
                        action="store_true",
                        help="Plot simple lines for the positioner arms and no "
                        "exclusion zones.  This is less informative, but faster.")

    parser.add_argument("--margin-pos", type=float, required=False, default=0.,
                        help="Add margin (in mm) around positioner keep-out polygons")
    parser.add_argument("--margin-petal", type=float, required=False, default=0.,
                        help="Add margin (in mm) around petal-boundary keep-out polygons")
    parser.add_argument("--margin-gfa", type=float, required=False, default=0.,
                        help="Add margin (in mm) around GFA keep-out polygons")

    parser.add_argument("--serial", required=False, default=False,
                        action="store_true",
                        help="Disable the use of multiprocessing.  Needed by "
                        "some unit tests to avoid issues with matplotlib.")

    parser.add_argument("files", nargs="+")

    args = None
    if optlist is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(optlist)

    # Set up margins dict
    args.margins = {}
    if args.margin_pos != 0:
        args.margins['theta'] = args.margin_pos
        args.margins['phi']   = args.margin_pos
    if args.margin_petal != 0:
        ## The petal polygon is "backwards", so have to give it a negative margin
        args.margins['petal'] = -args.margin_petal
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
    petals = None
    if args.petals is not None:
        petals = [int(x) for x in args.petals.split(",")]
    return petals


def run_plot(args):
    """Run plotting of fiberassign outputs.

    This uses the previously parsed options to read input data and make
    per-tile plots of the fiber assignment.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        None

    """
    petals = run_plot_init(args)
    plot_tiles(
        args.files,
        petals=petals,
        real_shapes=(not args.simple),
        serial=args.serial,
        out_dir=args.out_dir
    )
    return
