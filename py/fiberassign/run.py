# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.run
=====================

High-level functions for common types of assignment cases.

"""
from __future__ import absolute_import, division, print_function

import os
import sys
import argparse
from datetime import datetime
import json

import numpy as np

from desitarget.targetmask import desi_mask

from fiberassign.utils import GlobalTimers, Logger

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles

from fiberassign.gfa import get_gfa_targets

from fiberassign.targets import (str_to_target_type, TARGET_TYPE_SCIENCE,
                                 TARGET_TYPE_SKY, TARGET_TYPE_STANDARD,
                                 TARGET_TYPE_SAFE, Targets, TargetsAvailable,
                                 TargetTree, FibersAvailable,
                                 load_target_file,
                                 get_sciencemask, get_stdmask, get_skymask,
                                 get_safemask, get_excludemask
                                 )

from fiberassign.assign import (Assignment, write_assignment_fits,
                                result_path, merge_results, result_tiles)

from fiberassign.qa import qa_tiles

from fiberassign.vis import plot_tiles, plot_qa


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


def parse_survey(optlist=None):
    """Parse survey options.

    This parses either sys.argv or a list of strings passed in.  If passing
    an option list, you can create that more easily using the
    :func:`option_list` function.

    Args:
        optlist (list, optional): Optional list of arguments to parse instead
            of using sys.argv.

    Returns:
        (namespace):  an ArgumentParser namespace.

    """
    log = Logger.get()
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

    parser.add_argument("--dir", type=str, required=False, default=None,
                        help="Output directory.")

    parser.add_argument("--prefix", type=str, required=False,
                        default="fiberassign_",
                        help="Prefix of each file (before the <tile>.fits).")

    parser.add_argument("--split", required=False, default=False,
                        action="store_true",
                        help="Split output into tile prefix directories.")

    parser.add_argument("--standards_per_petal", type=int, required=False,
                        default=10, help="Required number of standards per"
                        " petal")

    parser.add_argument("--sky_per_petal", type=int, required=False,
                        default=40, help="Required number of sky targets per"
                        " petal")

    parser.add_argument("--write_all_targets", required=False, default=False,
                        action="store_true",
                        help="When writing target properties, write data "
                        "for all available targets, not just those which are "
                        "assigned.  This is convenient, but increases the "
                        "write time and the file size.")

    parser.add_argument("--overwrite", required=False, default=False,
                        action="store_true",
                        help="Overwrite any pre-existing output files")

    parser.add_argument("--mask_column", required=False, default="DESI_TARGET",
                        help="Default FITS column to use for applying target "
                             "masks")

    parser.add_argument("--sciencemask", required=False,
                        default=get_sciencemask(),
                        help="Default DESI_TARGET mask to use for science "
                             "targets")

    parser.add_argument("--stdmask", required=False,
                        default=get_stdmask(),
                        help="Default DESI_TARGET mask to use for stdstar "
                             "targets")

    parser.add_argument("--skymask", required=False,
                        default=get_skymask(),
                        help="Default DESI_TARGET mask to use for sky targets")

    parser.add_argument("--safemask", required=False,
                        default=get_safemask(),
                        help="Default DESI_TARGET mask to use for safe "
                        "backup targets")

    parser.add_argument("--excludemask", required=False,
                        default=get_excludemask(),
                        help="Default DESI_TARGET mask to exclude from "
                        "any assignments")

    parser.add_argument("--by_tile", required=False, default=False,
                        action="store_true",
                        help="Do assignment one tile at a time.  This disables"
                        " redistribution.")

    args = None
    if optlist is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(optlist)

    # Allow sciencemask, stdmask, etc. to be int or string
    if isinstance(args.sciencemask, str):
        args.sciencemask = desi_mask.mask(args.sciencemask.replace(",", "|"))

    if isinstance(args.stdmask, str):
        args.stdmask = desi_mask.mask(args.stdmask.replace(",", "|"))

    if isinstance(args.skymask, str):
        args.skymask = desi_mask.mask(args.skymask.replace(",", "|"))

    if isinstance(args.safemask, str):
        args.safemask = desi_mask.mask(args.safemask.replace(",", "|"))

    if isinstance(args.excludemask, str):
        args.excludemask = desi_mask.mask(args.excludemask.replace(",", "|"))

    log.info("sciencemask {}".format(
        " ".join(desi_mask.names(args.sciencemask))))
    log.info("stdmask     {}".format(" ".join(desi_mask.names(args.stdmask))))
    log.info("skymask     {}".format(" ".join(desi_mask.names(args.skymask))))
    log.info("safemask    {}".format(" ".join(desi_mask.names(args.safemask))))
    log.info("excludemask {}".format(
        " ".join(desi_mask.names(args.excludemask))))

    # Get run date
    if args.rundate is None:
        args.rundate = datetime.utcnow()
    else:
        args.rundate = datetime.strptime(args.rundate, "%Y-%m-%dT%H:%M:%S")

    # Set output directory
    if args.dir is None:
        rdstr = args.rundate.isoformat("T", "seconds")
        args.dir = "out_fiberassign_{}".format(rdstr)

    return args


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

    parser.add_argument("--dir", type=str, required=True, default=None,
                        help="Directory containing fiberassign results.")

    parser.add_argument("--prefix", type=str, required=False,
                        default="fiberassign_",
                        help="Prefix of each file (before the <tile>.fits).")

    parser.add_argument("--split", required=False, default=False,
                        action="store_true",
                        help="Results are in tile prefix directories.")

    parser.add_argument("--out", type=str, required=False, default=None,
                        help="Output directory for the merged files.  Default"
                        " is the directory containing the fiberassign output.")

    parser.add_argument("--out_prefix", type=str, required=False,
                        default="tile-",
                        help="Prefix of each output file.")

    parser.add_argument("--out_split", required=False, default=False,
                        action="store_true",
                        help="Split output into tile prefix directories.")

    parser.add_argument("--columns", type=str, required=False, default=None,
                        help="Override the column names of target data to be "
                        "copied from the target files into the fiber "
                        "assignment files.  This should be a comma-separated "
                        "list.")

    parser.add_argument("--skip_raw", required=False, default=False,
                        action="store_true",
                        help="Do not copy the raw fiberassign HDUs to the "
                             "merged output.")

    args = None
    if optlist is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(optlist)

    return args


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
                        "YYYY-MM-DD or YYYY-MM-DDTHH:mm:ss in UTC time.")

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

    # Get run date
    if args.rundate is None:
        args.rundate = datetime.now()
    else:
        args.rundate = datetime.fromisoformat(args.rundate)

    return args


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
                        default="fiberassign_",
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

    parser.add_argument("--footprint", type=str, required=False, default=None,
                        help="Optional FITS file defining the footprint.  If"
                        " not specified, the default footprint from desimodel"
                        " is used.")

    parser.add_argument("--tiles", type=str, required=False, default=None,
                        help="Optional text file containing a subset of the"
                        " tile IDs to use in the footprint, one ID per line."
                        " Default uses all tiles in the footprint.")

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
                        "YYYY-MM-DD or YYYY-MM-DDTHH:mm:ss in UTC time.")

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

    # Get run date
    if args.rundate is None:
        args.rundate = datetime.now()
    else:
        args.rundate = datetime.fromisoformat(args.rundate)

    return args


def parse_plot_qa(optlist=None):
    """Parse QA plotting options.

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

    parser.add_argument("--qafile", type=str, required=True, default=None,
                        help="Input QA file.")

    parser.add_argument("--outroot", type=str, required=False, default=None,
                        help="Output root file name.  Default uses input.")

    parser.add_argument("--labels", required=False, default=False,
                        action="store_true",
                        help="Plot tile IDs at center of circles.")

    args = None
    if optlist is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(optlist)

    # Check directory
    if not os.path.isfile(args.qafile):
        raise RuntimeError("Input file {} does not exist".format(args.qafile))

    if args.outroot is None:
        args.outroot = os.path.splitext(args.qafile)[0]

    return args


def run_survey_init(args):
    """Initialize survey inputs.

    This uses the previously parsed options to load the input files needed.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        (tuple):  The (Hardware, Tiles, Targets) needed to run assignment.

    """
    log = Logger.get()
    # Read hardware properties
    hw = load_hardware(fiberpos_file=args.positioners, rundate=args.rundate,
                       status_file=args.status)

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

    # Before doing significant calculations, check for pre-existing files
    if not args.overwrite:
        for tileid in tiles.id:
            outfile = result_path(tileid, dir=args.dir,
                                  prefix=args.prefix, split=args.split)
            if os.path.exists(outfile):
                outdir = os.path.split(outfile)[0]
                log.error("Output files already exist in {}".format(outdir))
                log.error("either remove them or use --overwrite")
                sys.exit(1)

    # Create empty target list
    tgs = Targets()

    # Append each input target file
    for tgarg in args.targets:
        tgprops = tgarg.split(",")
        tgfile = tgprops[0]
        typeforce = None
        if len(tgprops) > 1:
            # we are forcing the target type for this file
            typeforce = str_to_target_type(tgprops[1])
        load_target_file(tgs, tgfile, typeforce=typeforce,
                         typecol=args.mask_column,
                         sciencemask=args.sciencemask,
                         stdmask=args.stdmask,
                         skymask=args.skymask,
                         safemask=args.safemask,
                         excludemask=args.excludemask)
    return (hw, tiles, tgs)


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


def run_qa_init(args):
    """Initialize QA inputs.

    This uses the previously parsed options to load the input files needed.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        (tuple):  The (Hardware, Tiles) needed to run QA.

    """
    # Read hardware properties
    hw = load_hardware(fiberpos_file=args.positioners, rundate=args.rundate,
                       status_file=args.status)

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


def run_plot_init(args):
    """Initialize plotting inputs.

    This uses the previously parsed options to load the input files needed.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        (tuple):  The (Hardware, Tiles, petals) needed to run plotting.

    """
    # Read hardware properties
    hw = load_hardware(fiberpos_file=args.positioners, rundate=args.rundate,
                       status_file=args.status)

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


def run_survey_full(args):
    """Run fiber assignment over all tiles simultaneously.

    This uses the previously parsed options to read input data and run through
    the typical assignment sequence doing one step at a time over all tiles.
    It then writes to the outputs to disk.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        None

    """
    gt = GlobalTimers.get()
    gt.start("run_survey_full calculation")

    # Load data
    hw, tiles, tgs = run_survey_init(args)

    # Create a hierarchical triangle mesh lookup of the targets positions
    tree = TargetTree(tgs)

    # Compute the targets available to each fiber for each tile.
    tgsavail = TargetsAvailable(hw, tgs, tiles, tree)

    # Free the tree
    del tree

    # Compute the fibers on all tiles available for each target and sky
    favail = FibersAvailable(tgsavail)

    # Create assignment object
    asgn = Assignment(tgs, tgsavail, favail)

    # First-pass assignment of science targets
    asgn.assign_unused(TARGET_TYPE_SCIENCE)

    # Redistribute science targets across available petals
    asgn.redistribute_science()

    # Assign standards, 10 per petal
    asgn.assign_unused(TARGET_TYPE_STANDARD, args.standards_per_petal)
    asgn.assign_force(TARGET_TYPE_STANDARD, args.standards_per_petal)

    # Assign sky to unused fibers, up to 40 per petal
    asgn.assign_unused(TARGET_TYPE_SKY, args.sky_per_petal)
    asgn.assign_force(TARGET_TYPE_SKY, args.sky_per_petal)

    # If there are any unassigned fibers, try to place them somewhere.
    asgn.assign_unused(TARGET_TYPE_SCIENCE)
    asgn.assign_unused(TARGET_TYPE_SKY)

    # NOTE:  This was removed since we are treating BAD_SKY as science targets
    # with very low priority.
    #
    # # Assign safe location to unused fibers (no maximum).  There should
    # # always be at least one safe location (i.e. "BAD_SKY") for each fiber.
    # # So after this is run every fiber should be assigned to something.
    # asgn.assign_unused(TARGET_TYPE_SAFE)

    # Assign sky monitor fibers
    asgn.assign_unused(TARGET_TYPE_SKY, -1, "ETC")

    gt.stop("run_survey_full calculation")
    gt.start("run_survey_full write output")

    # Make sure that output directory exists
    if not os.path.isdir(args.dir):
        os.makedirs(args.dir)

    # Optionally get GFA targets
    gfa_targets = None
    if args.gfafile is not None:
        gfa_targets = get_gfa_targets(tiles, args.gfafile)

    # Write output
    write_assignment_fits(tiles, asgn, out_dir=args.dir,
                          out_prefix=args.prefix, split_dir=args.split,
                          all_targets=args.write_all_targets,
                          gfa_targets=gfa_targets, overwrite=args.overwrite)

    gt.stop("run_survey_full write output")

    gt.report()

    return


def run_survey_bytile(args):
    """Run fiber assignment tile-by-tile.

    This uses the previously parsed options to read input data and run through
    the typical assignment sequence on a single tile before moving on to the
    next.  It then writes to the outputs to disk.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        None

    """
    gt = GlobalTimers.get()
    gt.start("run_survey_bytile calculation")

    # Load data
    hw, tiles, tgs = run_survey_init(args)

    # Create a hierarchical triangle mesh lookup of the targets positions
    tree = TargetTree(tgs)

    # Compute the targets available to each fiber for each tile.
    tgsavail = TargetsAvailable(hw, tgs, tiles, tree)

    # Free the tree
    del tree

    # Compute the fibers on all tiles available for each target and sky
    favail = FibersAvailable(tgsavail)

    # Create assignment object
    asgn = Assignment(tgs, tgsavail, favail)

    # We are now going to loop over tiles and assign each one fully before
    # moving on to the next

    for tile_id in tiles.id:

        # First-pass assignment of science targets
        asgn.assign_unused(TARGET_TYPE_SCIENCE, -1, "POS", tile_id, tile_id)

        # Assign standards, 10 per petal
        asgn.assign_unused(TARGET_TYPE_STANDARD, args.standards_per_petal,
                           "POS", tile_id, tile_id)
        asgn.assign_force(TARGET_TYPE_STANDARD, args.standards_per_petal,
                          tile_id, tile_id)

        # Assign sky to unused fibers, up to 40 per petal
        asgn.assign_unused(TARGET_TYPE_SKY, args.sky_per_petal, "POS",
                           tile_id, tile_id)
        asgn.assign_force(TARGET_TYPE_SKY, args.sky_per_petal,
                          tile_id, tile_id)

        # NOTE:  This was removed since we are treating BAD_SKY as science
        # targets with very low priority.
        #
        # # Assign safe location to unused fibers (no maximum).  There should
        # # always be at least one safe location (i.e. "BAD_SKY") for each
        # # fiber.  So after this is run every fiber should be assigned to
        # # something.
        # asgn.assign_unused(TARGET_TYPE_SAFE)

        # Assign sky monitor fibers
        asgn.assign_unused(TARGET_TYPE_SKY, -1, "ETC", tile_id, tile_id)

    gt.stop("run_survey_bytile calculation")
    gt.start("run_survey_bytile write output")

    # Make sure that output directory exists
    if not os.path.isdir(args.dir):
        os.makedirs(args.dir)

    # Optionally get GFA targets
    gfa_targets = None
    if args.gfafile is not None:
        gfa_targets = get_gfa_targets(tiles, args.gfafile)

    # Write output
    write_assignment_fits(tiles, asgn, out_dir=args.dir,
                          out_prefix=args.prefix, split_dir=args.split,
                          all_targets=args.write_all_targets,
                          gfa_targets=gfa_targets, overwrite=args.overwrite)

    gt.stop("run_survey_bytile write output")

    gt.report()

    return


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
    merge_results(args.targets, tiles, result_dir=args.dir,
                  result_prefix=args.prefix, result_split_dir=args.split,
                  out_dir=args.out, out_prefix=args.out_prefix,
                  out_split_dir=args.out_split, columns=columns,
                  copy_fba=(not args.skip_raw))
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
        merge_results(args.targets, ptiles, result_dir=args.dir,
                      result_prefix=args.prefix, out_dir=args.out,
                      out_prefix=args.out_prefix, columns=columns,
                      copy_fba=(not args.skip_raw))
    return


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


def run_plot_qa(args):
    """Run QA plotting.

    This uses the previously parsed options to read input data and make a
    plot of the QA results.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        None

    """
    qadata = None
    with open(args.qafile, "r") as f:
        qadata = json.load(f)
    plot_qa(qadata, args.outroot, labels=args.labels)
    return
