# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.scripts.fulltest
==============================

High-level functions for testing all commandline operations.

"""
from __future__ import absolute_import, division, print_function

import os
import argparse
import shutil

from .assign import (parse_assign, run_assign_full)

from .plot import (parse_plot, run_plot)

from .qa import (parse_qa, run_qa)

from .qa_plot import (parse_plot_qa, run_plot_qa)

from .merge import (parse_merge, run_merge)

from ..utils import option_list, Environment


def parse_fulltest(optlist=None):
    """Parse fulltest options.

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
                        "in different files).  By default, the "
                        "'--mask_column' (default DESI_TARGET)"
                        "column and bitfield values defined in desitarget "
                        "are used to determine the type of each target.  "
                        "Each filename may be optionally followed by comma "
                        "and then one of the strings 'science', 'standard', "
                        "'sky' or 'safe' to force all targets in that file "
                        "to be treated as a fixed target type.")

    parser.add_argument("--sky", type=str, required=False, nargs="+",
                        help="Input file with sky or 'bad sky' targets.  "
                        "This option exists in order to treat main-survey"
                        " sky target files as valid for other survey types."
                        "  If you are running a main survey assignment, you"
                        " can just pass the sky file to the --targets list.")

    parser.add_argument("--gfafile", type=str, required=False, default=None,
                        help="Optional GFA targets FITS file")

    parser.add_argument("--footprint", type=str, required=False, default=None,
                        help="Optional FITS file defining the footprint.  If"
                        " not specified, the default footprint from desimodel"
                        " is used.")

    parser.add_argument("--tiles", type=str, required=False, default=None,
                        help="Optional text file containing a subset of the"
                        " tile IDs to use in the footprint, one ID per line."
                        " Default uses all tiles in the footprint.")

    parser.add_argument("--out", type=str, required=False, default=None,
                        help="Top-level output directory.")

    parser.add_argument("--plotpetals", type=str, required=False, default=None,
                        help="Comma-separated list of petals to plot "
                        "(default is all petals)")

    parser.add_argument("--positioners", type=str, required=False,
                        default=None,
                        help="Optional FITS file describing the fiber "
                        "positioner locations.  Default uses the file from "
                        "desimodel.")

    parser.add_argument("--status", type=str, required=False, default=None,
                        help="Optional fiber status file in astropy ECSV "
                        "format.  Default treats all fibers as good.")

    parser.add_argument("--rundate", type=str, required=False, default=None,
                        help="Optional date to simulate for this run of "
                        "fiber assignment, used with the fiber status file "
                        "to determine which fibers currently have problems.  "
                        "Default uses the current date.  Format is "
                        "YYYY-MM-DDTHH:mm:ss in UTC time.")

    parser.add_argument("--standards_per_petal", type=int, required=False,
                        default=10, help="Required number of standards per"
                        " petal")

    parser.add_argument("--sky_per_petal", type=int, required=False,
                        default=40, help="Required number of sky targets per"
                        " petal")

    args = None
    if optlist is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(optlist)

    # Set output directory
    if args.out is None:
        if args.rundate is None:
            raise RuntimeError(
                "You must specify the output directory or the rundate")
        args.out = "out_fulltest_{}".format(args.rundate)

    return args


def run_fulltest(args):
    """Run a full test of all commandline tools.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        None

    """
    # output directories
    outroot = args.out
    rawout = os.path.join(outroot, "raw")
    rawplot = os.path.join(outroot, "raw_plots")
    mergedout = os.path.join(outroot, "merged")
    mergedplot = os.path.join(outroot, "merged_plots")

    # Assignment

    opts = dict()
    opts["targets"] = args.targets
    opts["sky"] = args.sky
    opts["gfafile"] = args.gfafile
    opts["footprint"] = args.footprint
    opts["tiles"] = args.tiles
    opts["positioners"] = args.positioners
    opts["status"] = args.status
    opts["rundate"] = args.rundate
    opts["standards_per_petal"] = args.standards_per_petal
    opts["sky_per_petal"] = args.sky_per_petal
    opts["dir"] = rawout
    opts["write_all_targets"] = True
    opts["overwrite"] = True

    options = option_list(opts)
    parsed = parse_assign(options)
    run_assign_full(parsed)

    # Set the number of OpenMP threads to one for the rest of the steps.
    # The steps below use multiprocessing instead for parallelism.

    env = Environment.get()
    env.set_threads(1)

    # QA on raw assignment outputs

    opts = dict()
    opts["footprint"] = args.footprint
    opts["dir"] = rawout
    opts["prefix"] = "fiberassign_"

    options = option_list(opts)
    parsed = parse_qa(options)
    run_qa(parsed)

    # QA plot on raw assignment outputs

    opts = dict()
    opts["qafile"] = os.path.join(rawout, "qa.json")

    options = option_list(opts)
    parsed = parse_plot_qa(options)
    run_plot_qa(parsed)

    # Tile plots on raw assignment outputs

    if os.path.exists(rawplot):
        shutil.rmtree(rawplot, ignore_errors=True)

    opts = dict()
    opts["footprint"] = args.footprint
    opts["dir"] = rawout
    opts["prefix"] = "fiberassign_"
    opts["out"] = rawplot
    opts["petals"] = args.plotpetals

    options = option_list(opts)
    parsed = parse_plot(options)
    run_plot(parsed)

    # Create merged outputs.  Do not propagate the raw HDUs, so that we
    # can test that the following steps work on such files.

    opts = dict()
    opts["targets"] = args.targets
    opts["sky"] = args.sky
    opts["dir"] = rawout
    opts["out"] = mergedout
    opts["skip_raw"] = True

    options = option_list(opts)
    parsed = parse_merge(options)
    run_merge(parsed)

    # QA on merged outputs

    opts = dict()
    opts["footprint"] = args.footprint
    opts["dir"] = mergedout
    opts["prefix"] = "tile-"

    options = option_list(opts)
    parsed = parse_qa(options)
    run_qa(parsed)

    # QA plot on merged outputs

    opts = dict()
    opts["qafile"] = os.path.join(mergedout, "qa.json")

    options = option_list(opts)
    parsed = parse_plot_qa(options)
    run_plot_qa(parsed)

    # Tile plots on merged outputs

    if os.path.exists(mergedplot):
        shutil.rmtree(mergedplot, ignore_errors=True)

    opts = dict()
    opts["footprint"] = args.footprint
    opts["dir"] = mergedout
    opts["prefix"] = "tile-"
    opts["out"] = mergedplot
    opts["petals"] = args.plotpetals

    options = option_list(opts)
    parsed = parse_plot(options)
    run_plot(parsed)

    return
