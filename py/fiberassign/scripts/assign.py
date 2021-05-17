# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.scripts.assign
==============================

High-level functions for running assignment.

"""
from __future__ import absolute_import, division, print_function

import os
import sys
import argparse
import re

from ..utils import GlobalTimers, Logger

from ..hardware import load_hardware

from ..tiles import load_tiles

from ..gfa import get_gfa_targets

from ..targets import (str_to_target_type, TARGET_TYPE_SCIENCE,
                       TARGET_TYPE_SKY, TARGET_TYPE_SUPPSKY,
                       TARGET_TYPE_STANDARD,
                       TARGET_TYPE_SAFE, Targets, TargetsAvailable,
                       LocationsAvailable,
                       load_target_file, targets_in_tiles,
                       create_tagalong)

from ..assign import (Assignment, write_assignment_fits,
                      result_path, run)

from ..stucksky import stuck_on_sky

def parse_assign(optlist=None):
    """Parse assignment options.

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

    parser.add_argument("--rundate", type=str, required=False, default=None,
                        help="Optional date to simulate for this run of "
                        "fiber assignment, used to load the correct "
                        "focalplane properties and state from desimodel.  "
                        "Default uses the current date.  Format is "
                        "YYYY-MM-DDTHH:mm:ss+-zz:zz.")

    parser.add_argument("--obsdate", type=str, required=False, default="2022-07-01",
                        help="Plan field rotations for this date (YEARMMDD, "
                        "or ISO 8601 YEAR-MM-DD with or without time).")

    parser.add_argument("--ha", type=float, required=False, default=0.,
                        help="Design for the given Hour Angle in degrees.")

    parser.add_argument("--fieldrot", type=float, required=False, default=None,
                        help="Override obsdate and use this field rotation "
                        "for all tiles (degrees counter clockwise in CS5)")

    parser.add_argument("--dir", type=str, required=False, default=None,
                        help="Output directory.")

    parser.add_argument("--prefix", type=str, required=False,
                        default="fba-",
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

    parser.add_argument("--sky_per_slitblock", type=int, required=False,
                        default=1, help="Required number of sky targets per"
                        " fiber slitblock")

    parser.add_argument("--margin-pos", type=float, required=False, default=0.05,
                        help="Add margin (in mm) around positioner keep-out polygons")
    parser.add_argument("--margin-petal", type=float, required=False, default=0.4,
                        help="Add margin (in mm) around petal-boundary keep-out polygons")
    parser.add_argument("--margin-gfa", type=float, required=False, default=0.4,
                        help="Add margin (in mm) around GFA keep-out polygons")

    parser.add_argument("--write_all_targets", required=False, default=False,
                        action="store_true",
                        help="When writing target properties, write data "
                        "for all available targets, not just those which are "
                        "assigned.  This is convenient, but increases the "
                        "write time and the file size.")

    parser.add_argument("--overwrite", required=False, default=False,
                        action="store_true",
                        help="Overwrite any pre-existing output files")

    parser.add_argument("--mask_column", required=False, default=None,
                        help="FITS column to use for applying target "
                             "masks")

    parser.add_argument("--sciencemask", required=False,
                        default=None,
                        help="Default DESI_TARGET mask to use for science "
                             "targets")

    parser.add_argument("--stdmask", required=False,
                        default=None,
                        help="Default DESI_TARGET mask to use for stdstar "
                             "targets")

    parser.add_argument("--skymask", required=False,
                        default=None,
                        help="Default DESI_TARGET mask to use for sky targets")

    parser.add_argument("--safemask", required=False,
                        default=None,
                        help="Default DESI_TARGET mask to use for safe "
                        "backup targets")

    parser.add_argument("--excludemask", required=False,
                        default=None,
                        help="Default DESI_TARGET mask to exclude from "
                        "any assignments")

    parser.add_argument("--by_tile", required=False, default=False,
                        action="store_true",
                        help="Do assignment one tile at a time.  This disables"
                        " redistribution.")

    parser.add_argument("--no_redistribute", required=False, default=False,
                        action="store_true",
                        help="Disable redistribution of science targets.")

    parser.add_argument("--no_zero_obsremain", required=False, default=False,
                        action="store_true",
                        help="Disable oversubscription of science targets with leftover fibers.")

    args = None
    if optlist is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(optlist)

    if args.sky is None:
        args.sky = list()

    # If any of the masks are strings, determine the survey type from
    # the first target file to know which bitmask to use
    if isinstance(args.sciencemask, str) or \
       isinstance(args.stdmask, str) or \
       isinstance(args.skymask, str) or \
       isinstance(args.safemask, str) or \
       isinstance(args.excludemask, str):
        import fitsio
        from desitarget.targets import main_cmx_or_sv
        data = fitsio.read(args.targets[0], 1, rows=[0,1])
        filecols, filemasks, filesurvey = main_cmx_or_sv(data)
        desi_mask = filemasks[0]

        # convert str bit names -> int bit mask
        if isinstance(args.sciencemask, str):
            try:
                args.sciencemask = int(args.sciencemask)
            except ValueError:
                args.sciencemask = desi_mask.mask(args.sciencemask.replace(",","|"))

        if isinstance(args.stdmask, str):
            try:
                args.stdmask = int(args.stdmask)
            except ValueError:
                args.stdmask = desi_mask.mask(args.stdmask.replace(",", "|"))

        if isinstance(args.skymask, str):
            try:
                args.skymask = int(args.skymask)
            except ValueError:
                args.skymask = desi_mask.mask(args.skymask.replace(",", "|"))

        if isinstance(args.safemask, str):
            try:
                args.safemask = int(args.safemask)
            except ValueError:
                args.safemask = desi_mask.mask(args.safemask.replace(",", "|"))

        if isinstance(args.excludemask, str):
            try:
                args.excludemask = int(args.excludemask)
            except ValueError:
                args.excludemask = desi_mask.mask(args.excludemask.replace(",","|"))

    # convert YEARMMDD to YEAR-MM-DD to be ISO 8601 compatible
    if re.match('\d{8}', args.obsdate):
        year = args.obsdate[0:4]
        mm = args.obsdate[4:6]
        dd = args.obsdate[6:8]
        #- Note: ISO8601 does not require time portion
        args.obsdate = '{}-{}-{}'.format(year, mm, dd)

    # Set output directory
    if args.dir is None:
        if args.rundate is None:
            raise RuntimeError(
                "You must specify the output directory or the rundate")
        args.dir = "out_fiberassign_{}".format(args.rundate)

    # Set up margins dict
    args.margins = dict(pos=args.margin_pos,
                        petal=args.margin_petal,
                        gfa=args.margin_gfa)
    return args

def run_assign_init(args, plate_radec=True):
    """Initialize assignment inputs.

    This uses the previously parsed options to load the input files needed.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        (tuple):  The (Hardware, Tiles, Targets) needed to run assignment.

    """
    log = Logger.get()
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
    tiles = load_tiles(tiles_file=args.footprint, select=tileselect,
        obstime=args.obsdate, obstheta=args.fieldrot, obsha=args.ha)

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
    # Create structure for carrying along auxiliary target data not needed by C++.
    tagalong = create_tagalong(plate_radec=plate_radec)

    # Append each input target file.  These target files must all be of the
    # same survey type, and will set the Targets object to be of that survey.

    for tgarg in args.targets:
        tgprops = tgarg.split(",")
        tgfile = tgprops[0]
        typeforce = None
        if len(tgprops) > 1:
            # we are forcing the target type for this file
            typeforce = str_to_target_type(tgprops[1])
        load_target_file(tgs, tagalong, tgfile, typeforce=typeforce,
                         typecol=args.mask_column,
                         sciencemask=args.sciencemask,
                         stdmask=args.stdmask,
                         skymask=args.skymask,
                         safemask=args.safemask,
                         excludemask=args.excludemask)
    # Now load the sky target files.  These are main-survey files that we will
    # force to be treated as the survey type of the other target files.
    survey = tgs.survey()
    for tgarg in args.sky:
        load_target_file(tgs, tagalong, tgarg, survey=survey, typeforce=typeforce,
                         typecol=args.mask_column,
                         sciencemask=args.sciencemask,
                         stdmask=args.stdmask,
                         skymask=args.skymask,
                         safemask=args.safemask,
                         excludemask=args.excludemask)

    return (hw, tiles, tgs, tagalong)


def run_assign_full(args, plate_radec=True):
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
    gt.start("run_assign_full calculation")

    # Load data
    hw, tiles, tgs, tagalong = run_assign_init(args, plate_radec=plate_radec)

    # Create a hierarchical triangle mesh lookup of the targets positions
    gt.start("Compute targets locations in tile")
    tile_targetids, tile_x, tile_y = targets_in_tiles(hw, tgs, tiles, tagalong)
    gt.stop("Compute targets locations in tile")

    # Compute the targets available to each fiber for each tile.
    gt.start("Compute Targets Available")
    tgsavail = TargetsAvailable(hw, tiles, tile_targetids, tile_x, tile_y)
    gt.stop("Compute Targets Available")

    # Free the target locations
    del tile_targetids, tile_x, tile_y

    # Compute the fibers on all tiles available for each target and sky
    gt.start("Compute Locations Available")
    favail = LocationsAvailable(tgsavail)
    gt.stop("Compute Locations Available")

    # Find stuck positioners and compute whether they will land on acceptable
    # sky locations for each tile.
    gt.start("Compute Stuck locations on good sky")
    stucksky = stuck_on_sky(hw, tiles)
    if stucksky is None:
        # (the pybind code doesn't like None when a dict is expected...)
        stucksky = {}
    gt.stop("Compute Stuck locations on good sky")

    # Create assignment object
    gt.start("Construct Assignment")
    asgn = Assignment(tgs, tgsavail, favail, stucksky)
    gt.stop("Construct Assignment")

    run(
        asgn,
        args.standards_per_petal,
        args.sky_per_petal,
        args.sky_per_slitblock,
        redistribute=(not args.no_redistribute),
        use_zero_obsremain=(not args.no_zero_obsremain)
    )

    gt.stop("run_assign_full calculation")
    gt.start("run_assign_full write output")

    # Make sure that output directory exists
    if not os.path.isdir(args.dir):
        os.makedirs(args.dir)

    # Optionally get GFA targets
    gfa_targets = None
    if args.gfafile is not None:
        gfa_targets = get_gfa_targets(tiles, args.gfafile)

    # Write output
    write_assignment_fits(tiles, tagalong, asgn, out_dir=args.dir,
                          out_prefix=args.prefix, split_dir=args.split,
                          all_targets=args.write_all_targets,
                          gfa_targets=gfa_targets, overwrite=args.overwrite,
                          stucksky=stucksky)

    gt.stop("run_assign_full write output")

    gt.report()

    return


def run_assign_bytile(args):
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
    gt.start("run_assign_bytile calculation")

    # Load data
    hw, tiles, tgs, tagalong = run_assign_init(args)

    # Create a hierarchical triangle mesh lookup of the targets positions
    gt.start("Compute targets locations in tile")
    tile_targetids, tile_x, tile_y = targets_in_tiles(hw, tgs, tiles, tagalong)
    gt.stop("Compute targets locations in tile")

    # Compute the targets available to each fiber for each tile.
    gt.start("Compute Targets Available")
    tgsavail = TargetsAvailable(hw, tiles, tile_targetids, tile_x, tile_y)
    gt.stop("Compute Targets Available")

    # Free the target locations
    del tile_targetids, tile_x, tile_y

    # Compute the fibers on all tiles available for each target and sky
    gt.start("Compute Locations Available")
    favail = LocationsAvailable(tgsavail)
    gt.stop("Compute Locations Available")

    # Find stuck positioners and compute whether they will land on acceptable
    # sky locations for each tile.
    gt.start("Compute Stuck locations on good sky")
    stucksky = stuck_on_sky(hw, tiles)
    if stucksky is None:
        # (the pybind code doesn't like None when a dict is expected...)
        stucksky = {}
    gt.stop("Compute Stuck locations on good sky")

    # Create assignment object
    gt.start("Construct Assignment")
    asgn = Assignment(tgs, tgsavail, favail, stucksky)
    gt.stop("Construct Assignment")

    # We are now going to loop over tiles and assign each one fully before
    # moving on to the next

    for tile_id in tiles.id:
        run(
            asgn,
            args.standards_per_petal,
            args.sky_per_petal,
            args.sky_per_slitblock,
            start_tile=tile_id,
            stop_tile=tile_id,
            redistribute=(not args.no_redistribute),
            use_zero_obsremain=(not args.no_zero_obsremain)
        )

    gt.stop("run_assign_bytile calculation")
    gt.start("run_assign_bytile write output")

    # Make sure that output directory exists
    if not os.path.isdir(args.dir):
        os.makedirs(args.dir)

    # Optionally get GFA targets
    gfa_targets = None
    if args.gfafile is not None:
        gfa_targets = get_gfa_targets(tiles, args.gfafile)

    # Write output
    write_assignment_fits(tiles, tagalong, asgn, out_dir=args.dir,
                          out_prefix=args.prefix, split_dir=args.split,
                          all_targets=args.write_all_targets,
                          gfa_targets=gfa_targets, overwrite=args.overwrite,
                          stucksky=stucksky)

    gt.stop("run_assign_bytile write output")

    gt.report()

    return
