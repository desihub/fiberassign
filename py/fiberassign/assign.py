# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.assign
=====================

Functions and classes for computing and writing fiber assignment.

"""
from __future__ import absolute_import, division, print_function

import os

import re

import numpy as np

import multiprocessing as mp
from multiprocessing.sharedctypes import RawArray

from functools import partial

from collections import OrderedDict

import fitsio

from desiutil.depend import add_dependencies

import desimodel.focalplane

from desitarget.targetmask import desi_mask

from ._version import __version__

from .utils import Logger, Timer, default_mp_proc, GlobalTimers

from .targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY, TARGET_TYPE_SUPPSKY,
                      TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE, desi_target_type,
                      default_target_masks, default_survey_target_masks)

from .hardware import (FIBER_STATE_UNASSIGNED, FIBER_STATE_STUCK,
                       FIBER_STATE_BROKEN)

from ._internal import Assignment


# The columns and types that are written to FITS format.  The raw data has
# 3 HDUs, and these are designed to efficiently contain all information used
# internally to fiber assignment without duplicating information:
#
# FASSIGN
# - Per-location information, including assignment.  Sorted by location.
#
# FTARGETS
# - All available targets and their properties relevant to assignment.  Sorted
#   by target ID.
#
# FAVAIL
# - Target IDs available to each location.  Sorted by loc and then target ID.
#

results_assign_columns = OrderedDict([
    ("FIBER", "i4"),
    ("TARGETID", "i8"),
    ("LOCATION", "i4"),
    ("FIBERSTATUS", "i4"),
    ("LAMBDA_REF", "f4"),
    ("PETAL_LOC", "i2"),
    ("DEVICE_LOC", "i4"),
    ("DEVICE_TYPE", "a3"),
    ("TARGET_RA", "f8"),
    ("TARGET_DEC", "f8"),
    ("FA_TARGET", "i8"),
    ("FA_TYPE", "u1"),
    ("FIBERASSIGN_X", "f4"),
    ("FIBERASSIGN_Y", "f4"),
])

results_targets_columns = OrderedDict([
    ("TARGETID", "i8"),
    ("TARGET_RA", "f8"),
    ("TARGET_DEC", "f8"),
    ("FA_TARGET", "i8"),
    ("FA_TYPE", "u1"),
    ("PRIORITY", "i4"),
    ("SUBPRIORITY", "f8"),
    ("OBSCONDITIONS", "i4"),
])

results_avail_columns = OrderedDict([
    ("LOCATION", "i4"),
    ("FIBER", "i4"),
    ("TARGETID", "i8")
])


def result_tiles(dir=".", prefix="fba-"):
    # Find all the per-tile files and get the tile IDs
    tiles = list()
    for root, dirs, files in os.walk(dir):
        for f in files:
            mat = re.match(r"{}(\d+).fits".format(prefix), f)
            if mat is not None:
                # Matches the prefix
                tiles.append(int(mat.group(1)))
    return tiles


def result_path(tile_id, dir=".", prefix="fba-",
                ext="fits", create=False, split=False):
    tiledir = dir
    if split:
        tilegroup = tile_id // 1000
        tiledir = os.path.join(dir, "{:02d}".format(tilegroup))
    if create:
        os.makedirs(tiledir, exist_ok=True)
    path = os.path.join(tiledir,
                        "{}{:06d}.{}".format(prefix, tile_id, ext))
    return path


def write_assignment_fits_tile(asgn, fulltarget, overwrite, params):
    """Write a single tile assignment to a FITS file.

    Args:
        outroot (str):  full path of the output root file name.
        asgn (Assignment):  the assignment class instance.
        fulltarget (bool):  if True, dump the target information for all
            available targets, not just the ones that are assigned.
        overwrite (bool): overwrite output files or not
        params (tuple):  tuple containing the tile ID, RA, DEC, rotation,
            output path, and GFA targets

    Returns:
        None

    """
    tm = Timer()
    tm.start()
    tile_id, tile_ra, tile_dec, tile_obstheta, tile_obstime, tile_obsha, \
        tile_file, gfa_targets = params
    log = Logger.get()

    # Hardware properties
    hw = asgn.hardware()

    # Targets available for all tile / locs
    tgsavail = asgn.targets_avail()

    # Target properties
    tgs = asgn.targets()

    # Data for this tile
    tdata = asgn.tile_location_target(tile_id)
    avail = tgsavail.tile_data(tile_id)

    # The recarray dtypes
    assign_dtype = np.dtype([(x, y) for x, y in
                            results_assign_columns.items()])
    targets_dtype = np.dtype([(x, y) for x, y in
                             results_targets_columns.items()])
    avail_dtype = np.dtype([(x, y) for x, y in
                            results_avail_columns.items()])
    # tm.stop()
    # tm.report("  data pointers for tile {}".format(tile_id))
    # tm.clear()
    # tm.start()

    locs = np.array(hw.locations)
    nloc = len(locs)

    # Compute the total list of targets
    navail = np.sum([len(avail[x]) for x in avail.keys()])
    tgids_avail = np.unique(np.concatenate([np.array(avail[x], dtype=np.int64)
                                            for x in avail.keys()]))

    # tm.stop()
    # tm.report("  available targets tile {}".format(tile_id))

    if len(tgids_avail) > 0:
        # We have some available targets.
        # tm.clear()
        # tm.start()

        if os.path.isfile(tile_file):
            if overwrite:
                log.warning("Overwriting {}".format(tile_file))
                os.remove(tile_file)
            else:
                raise RuntimeError("output file {} already exists"
                                   .format(tile_file))
        # This tile has some available targets
        log.info("Writing tile {}".format(tile_id))
        tmp_file = tile_file + ".tmp"
        fd = fitsio.FITS(tmp_file, "rw")

        # tm.stop()
        # tm.report("  opening file for tile {}".format(tile_id))
        # tm.clear()
        # tm.start()

        # Unpack all our target properties from C++ objects into numpy arrays

        tgids = None
        if fulltarget:
            # We are dumping all targets
            tgids = tgids_avail
        else:
            # We are dumping only assigned targets
            tgids = np.array(sorted([y for x, y in tdata.items()]),
                             dtype=np.int64)

        ntarget = len(tgids)

        tg_ra = np.empty(ntarget, dtype=np.float64)
        tg_dec = np.empty(ntarget, dtype=np.float64)
        tg_bits = np.zeros(ntarget, dtype=np.int64)
        tg_x = np.empty(ntarget, dtype=np.float64)
        tg_y = np.empty(ntarget, dtype=np.float64)
        tg_type = np.empty(ntarget, dtype=np.uint8)
        tg_priority = np.empty(ntarget, dtype=np.int32)
        tg_subpriority = np.empty(ntarget, dtype=np.float64)
        tg_obscond = np.empty(ntarget, dtype=np.int32)
        tg_indx = dict()
        for indx, tg in enumerate(tgids):
            tg_indx[tg] = indx
            props = tgs.get(tg)
            tg_ra[indx] = props.ra
            tg_dec[indx] = props.dec
            tg_bits[indx] = props.bits
            tg_type[indx] = props.type
            tg_priority[indx] = props.priority
            tg_subpriority[indx] = props.subpriority
            tg_obscond[indx] = props.obscond

        # We compute the X / Y focalplane coordinates for ALL available
        # targets, not just the assigned ones.  This allows us to write out
        # the focalplane coordinates of all available targets as computed by
        # the code at the time it was run.
        #
        # NOTE:  The output format is explicitly CS5 coordinates, even though
        # we use curved focal surface internally.
        xy = hw.radec2xy_multi(
            tile_ra, tile_dec, tile_obstheta, tg_ra, tg_dec, True, 0
        )
        for indx, fxy in enumerate(xy):
            tg_x[indx] = fxy[0]
            tg_y[indx] = fxy[1]

        # tm.stop()
        # tm.report("  extract target props for {}".format(tile_id))
        # tm.clear()
        # tm.start()

        # Write a header-only HDU 0 with some basic keywords
        header = dict()
        header["TILEID"] = tile_id
        header["TILERA"] = tile_ra
        header["TILEDEC"] = tile_dec
        header["FIELDROT"] = tile_obstheta
        header["FA_PLAN"] = tile_obstime
        header["FA_HA"] = tile_obsha
        header["FA_RUN"] = hw.time()

        header["REQRA"] = tile_ra
        header["REQDEC"] = tile_dec
        header["FIELDNUM"] = 0
        header["FA_VER"] = __version__
        header["FA_SURV"] = tgs.survey()
        add_dependencies(
            header,
            module_names=[
                "numpy",
                "matplotlib",
                "astropy",
                "fitsio",
                "desiutil",
                "desimodel",
                "desitarget"
            ]
        )
        fd.write(None, header=header, extname="PRIMARY")

        # FIXME:  write "-1" for unassigned targets.  Write all other fiber
        # status and other properties to this table.

        log.debug("Write:  copying assignment data for tile {}"
                  .format(tile_id))

        fdata = np.zeros(nloc, dtype=assign_dtype)
        fdata["LOCATION"] = locs

        # For unassigned fibers, we give each location a unique negative
        # number based on the tile and loc.
        unassign_offset = tile_id * nloc
        assigned_tgids = np.array([tdata[x] if x in tdata.keys()
                                  else -(unassign_offset + x)
                                  for x in locs], dtype=np.int64)
        fdata["TARGETID"] = assigned_tgids

        # Rows containing assigned locations
        assigned_valid = np.where(assigned_tgids >= 0)[0]
        assigned_invalid = np.where(assigned_tgids < 0)[0]

        # Buffers for X/Y/RA/DEC
        assigned_tgx = np.full(nloc, 9999.9, dtype=np.float64)
        assigned_tgy = np.full(nloc, 9999.9, dtype=np.float64)
        assigned_tgra = np.full(nloc, 9999.9, dtype=np.float64)
        assigned_tgdec = np.full(nloc, 9999.9, dtype=np.float64)
        assigned_tgbits = np.zeros(nloc, dtype=np.int64)
        assigned_tgtype = np.zeros(nloc, dtype=np.uint8)

        if (len(assigned_invalid) > 0):
            # Fill our unassigned location X/Y coordinates with the central
            # positioner locations.  Then convert these to RA/DEC.
            # NOTE:  Positioner locations are in curved focal surface coordinates.
            empty_fibers = locs[assigned_invalid]
            fpos_xy_mm = dict(hw.loc_pos_curved_mm)
            empty_x = np.array(
                [fpos_xy_mm[f][0] for f in empty_fibers], dtype=np.float64)
            empty_y = np.array(
                [fpos_xy_mm[f][1] for f in empty_fibers], dtype=np.float64)
            radec = hw.xy2radec_multi(
                tile_ra, tile_dec, tile_obstheta, empty_x, empty_y, False, 0
            )
            assigned_tgra[assigned_invalid] = [x for x, y in radec]
            assigned_tgdec[assigned_invalid] = [y for x, y in radec]
            empty_xy = hw.radec2xy_multi(
                tile_ra, tile_dec, tile_obstheta, assigned_tgra[assigned_invalid],
                assigned_tgdec[assigned_invalid], True, 0
            )
            assigned_tgx[assigned_invalid] = [x for x, y in empty_xy]
            assigned_tgy[assigned_invalid] = [y for x, y in empty_xy]

        if (len(assigned_valid) > 0):
            # The target IDs assigned to fibers (note- NOT sorted)
            assigned_real = np.copy(assigned_tgids[assigned_valid])

            # Mapping of our assigned target rows into the target properties.
            target_rows = [tg_indx[x] for x in assigned_real]

            # Copy values out of the full target list.
            assigned_tgra[assigned_valid] = np.array(tg_ra[target_rows])
            assigned_tgdec[assigned_valid] = np.array(tg_dec[target_rows])
            assigned_tgx[assigned_valid] = np.array(tg_x[target_rows])
            assigned_tgy[assigned_valid] = np.array(tg_y[target_rows])
            assigned_tgtype[assigned_valid] = np.array(tg_type[target_rows])
            assigned_tgbits[assigned_valid] = np.array(tg_bits[target_rows])

        fdata["TARGET_RA"] = assigned_tgra
        fdata["TARGET_DEC"] = assigned_tgdec
        fdata["FIBERASSIGN_X"] = assigned_tgx
        fdata["FIBERASSIGN_Y"] = assigned_tgy
        fdata["FA_TARGET"] = assigned_tgbits
        fdata["FA_TYPE"] = assigned_tgtype

        fibers = dict(hw.loc_fiber)
        etcloc = sorted([x for x, y in fibers.items() if y < 0])
        fakefibers = {y: (5000 + x) for x, y in enumerate(etcloc)}
        fullfibers = fibers
        fullfibers.update(fakefibers)
        device = dict(hw.loc_device)
        petal = dict(hw.loc_petal)
        device_type = dict(hw.loc_device_type)

        fdata["FIBER"] = np.array(
            [fullfibers[x] for x in locs]).astype(np.int32)
        fdata["DEVICE_LOC"] = np.array(
            [device[x] for x in locs]).astype(np.int32)
        fdata["PETAL_LOC"] = np.array(
            [petal[x] for x in locs]).astype(np.int16)
        fdata["DEVICE_TYPE"] = np.array(
            [device_type[x] for x in locs]).astype(np.dtype("a3"))

        # This hard-coded value copied from the original code...
        lambda_ref = np.ones(nloc, dtype=np.float32) * 5400.0
        fdata["LAMBDA_REF"] = lambda_ref

        # Fiber status
        fstate = dict(hw.state)
        fstatus = np.zeros(nloc, dtype=np.int32)
        # Set unused bit
        fstatus |= [0 if x in tdata.keys() else FIBER_STATE_UNASSIGNED
                    for x in locs]
        # Set stuck / broken bits
        fstatus |= [2 if (fstate[x] & FIBER_STATE_STUCK) else 0
                    for x in locs]
        fstatus |= [4 if (fstate[x] & FIBER_STATE_BROKEN) else 0
                    for x in locs]
        fstatus[assigned_valid] |= \
            [8 if (tg_type[x] & TARGET_TYPE_SAFE) else 0
             for x in target_rows]
        fdata["FIBERSTATUS"] = fstatus

        # tm.stop()
        # tm.report("  copy and compute assign data tile {}".format(tile_id))
        # tm.clear()
        # tm.start()

        log.debug("Write:  writing assignment data for tile {}"
                  .format(tile_id))

        fd.write(fdata, header=header, extname="FASSIGN")
        del fdata

        # tm.stop()
        # tm.report("  write assign data tile {}".format(tile_id))
        # tm.clear()
        # tm.start()

        log.debug("Write:  writing target data for tile {}"
                  .format(tile_id))

        fdata = np.zeros(ntarget, dtype=targets_dtype)
        fdata["TARGETID"] = tgids
        fdata["TARGET_RA"] = tg_ra
        fdata["TARGET_DEC"] = tg_dec
        fdata["FA_TARGET"] = tg_bits
        fdata["FA_TYPE"] = tg_type
        fdata["PRIORITY"] = tg_priority
        fdata["SUBPRIORITY"] = tg_subpriority
        fdata["OBSCONDITIONS"] = tg_obscond

        # tm.stop()
        # tm.report("  copy targets data tile {}".format(tile_id))
        # tm.clear()
        # tm.start()

        fd.write(fdata, header=header, extname="FTARGETS")
        del fdata

        # tm.stop()
        # tm.report("  write targets data tile {}".format(tile_id))
        # tm.clear()
        # tm.start()

        log.debug("Write:  writing avail data for tile {}"
                  .format(tile_id))

        fdata = np.zeros(navail, dtype=avail_dtype)
        off = 0
        for lid in sorted(avail.keys()):
            for tg in avail[lid]:
                fdata[off] = (lid, fibers[lid], tg)
                off += 1

        # tm.stop()
        # tm.report("  copy avail data tile {}".format(tile_id))
        # tm.clear()
        # tm.start()

        fd.write(fdata, header=header, extname="FAVAIL")
        del fdata

        if gfa_targets is not None:
            try:
                # Astropy Table
                fd.write(gfa_targets.as_array(), extname="GFA_TARGETS")
            except AttributeError:
                # numpy structured array
                fd.write(gfa_targets, extname="GFA_TARGETS")

        fd.close()
        os.rename(tmp_file, tile_file)

        # tm.stop()
        # tm.report("  write avail data tile {}".format(tile_id))
    return


def write_assignment_fits(tiles, asgn, out_dir=".", out_prefix="fba-",
                          split_dir=False, all_targets=False,
                          gfa_targets=None, overwrite=False):
    """Write out assignment results in FITS format.

    For each tile, all available targets (not only the assigned targets) and
    their properties are written to the first HDU.  The second HDU contains
    one row for every location and target available to that location.

    Args:
        tiles (Tiles):  The Tiles object containing the properties of each
            tile.
        asgn (Assignment):  The assignment object.
        out_dir (str):  The output directory for writing per-tile files.
        out_prefix (str):  The output file name prefix.
        split_dir (bool):  Optionally split files by tile ID prefix.
        all_targets (bool):  Optionally dump the target properties of all
            available targets for each tile.  If False, only dump the target
            properties of assigned targets.
        gfa_targets (list of numpy arrays): Include these as GFA_TARGETS HDUs
        overwrite (bool): overwrite pre-existing output files

    Returns:
        None

    """
    tm = Timer()
    tm.start()

    # The tile IDs that were assigned
    tileids = asgn.tiles_assigned()
    tileorder = tiles.order
    tilera = tiles.ra
    tiledec = tiles.dec
    tiletheta = tiles.obstheta
    tiletime = tiles.obstime
    tileha = tiles.obshourang

    write_tile = partial(write_assignment_fits_tile,
                         asgn, all_targets, overwrite)

    for i, tid in enumerate(tileids):
        tra = tilera[tileorder[tid]]
        tdec = tiledec[tileorder[tid]]
        ttheta = tiletheta[tileorder[tid]]
        ttime = tiletime[tileorder[tid]]
        tha = tileha[tileorder[tid]]

        outfile = result_path(tid, dir=out_dir, prefix=out_prefix,
                              create=True, split=split_dir)
        if gfa_targets is None:
            params = (tid, tra, tdec, ttheta, ttime, tha, outfile, None)
        else:
            params = (
                tid, tra, tdec, ttheta, ttime, tha, outfile, gfa_targets[i]
            )

        write_tile(params)

    tm.stop()
    tm.report("Write output files")

    return


def write_assignment_ascii(tiles, asgn, out_dir=".", out_prefix="fba-",
                           split_dir=False):
    """Write out assignment results in ASCII format.

    For each tile, only the final assignment to each tile is written out.  For
    the full information, including available targets, use the FITS format.

    Args:
        tiles (Tiles):  The Tiles object containing the properties of each
            tile.
        asgn (Assignment):  The assignment object.
        out_dir (str):  The output directory for writing per-tile files.
        out_prefix (str):  The output file name prefix.
        split_dir (bool):  Optionally split files by tile ID prefix.

    Returns:
        None

    """
    log = Logger.get()
    # Go through the assignment, one tile at a time.  For each tile, get the
    # best assignment and potential targets.

    tileids = asgn.tiles_assigned()
    tileorder = tiles.order
    tilera = tiles.ra
    tiledec = tiles.dec
    # Target properties
    tgs = asgn.targets()

    for t in tileids:
        tdata = asgn.tile_location_target(t)
        nloc = len(tdata)
        tfile = result_path(t, dir=out_dir, prefix=out_prefix,
                            create=True, split=split_dir)
        if nloc > 0:
            log.debug("Writing tile {}".format(t))
            with open(tfile, "w") as f:
                f.write("# TILE_RA = {}\n".format(tilera[tileorder[t]]))
                f.write("# TILE_DEC = {}\n".format(tiledec[tileorder[t]]))
                f.write("#\n")
                f.write("# LOCATION  TARGETID  RA  DEC  PRIORITY  "
                        "SUBPRIORITY  OBSCONDITIONS  FBATYPE\n")
                for lid in sorted(tdata.keys()):
                    tgid = tdata[lid]
                    tg = tgs.get(tgid)
                    f.write("{:d} {:d} {:.6f} {:.6f}\n"
                            .format(lid, tgid, tg.ra, tg.dec, tg.priority,
                                    tg.subpriority, tg.obscond, tg.type))
    return


def avail_table_to_dict(avail_data):
    """Convert a recarray of available targets into a dictionary.

    Args:
        avail_data (array):  The available targets read from a FITS HDU
            (for example).

    Returns:
        (dict):  A dictionary of numpy arrays, one per location, containing the
            available target IDs for the location.

    """
    avail_target = avail_data["TARGETID"]
    avail_loc = avail_data["LOCATION"]
    avail = dict()
    for lid, tgid in zip(avail_loc, avail_target):
        if lid in avail:
            avail[lid].append(tgid)
        else:
            avail[lid] = list([tgid])
    avail = {f: np.array(av) for f, av in avail.items()}
    return avail

def gfa_table_to_dict(gfa_data):
    """Convert a recarray of gfa targets into a dictionary.

    Args:
        gfa_data (array):  The available targets read from a FITS HDU
            (for example).

    Returns:
        (dict):  A dictionary of numpy arrays, one per location, containing the
            available gfa target IDs for the location.

    """
    gfa_target = gfa_data["TARGETID"]
    gfa_loc = gfa_data["GFA_LOC"]
    gfa_gmag = gfa_data["GAIA_PHOT_G_MEAN_MAG"]
    gfa = dict()
    for lid, tgid, mag in zip(gfa_loc, gfa_target,gfa_gmag):
        print(zip(gfa_loc, gfa_target,gfa_gmag))
        if lid in gfa:
            gfa[lid].append(mag)
        else:
            gfa[lid] = list([mag])
    gfa = {f: np.array(av) for f, av in gfa.items()}
    return gfa

def read_assignment_fits_tile(params):
    """Read in results.

    This reads in only the result information that was originally written by
    the fiber assignment.  This function is used internally for reading data
    for use in plotting and for reading the original data for merging with
    input target catalogs.

    Args:
        params (tuple):  The tile ID and tile file path packed as a tuple for
            use in multiprocessing.

    Returns:
        (tuple):  The FITS header, assignment recarray, target property
            recarray, the available targets recarray, and the GFA targets

    """
    (tile_id, tile_file) = params
    log = Logger.get()
    # Output tile file
    if not os.path.isfile(tile_file):
        raise RuntimeError("input file {} does not exist".format(tile_file))

    log.debug("Reading tile data {}".format(tile_file))

    # Open the file
    fd = fitsio.FITS(tile_file, "r")

    header = None
    fiber_data = None
    targets_data = None
    avail_data = None
    if "FIBERASSIGN" in fd:
        # We have merged data.  Build new recarrays from the FIBERASSIGN,
        # TARGETS, and POTENTIAL_ASSIGNMENTS extensions.
        header = fd["FIBERASSIGN"].read_header()
        npos = fd["FIBERASSIGN"].get_nrows()
        fbassign = fd["FIBERASSIGN"].read()
        netc = fd["SKY_MONITOR"].get_nrows()
        fbsky = fd["SKY_MONITOR"].read()
        nfiber = npos + netc
        assign_dtype = np.dtype([(x, y) for x, y in
                                results_assign_columns.items()])
        fiber_data = np.zeros(nfiber, dtype=assign_dtype)

        survey = None
        if "FA_SURV" in header:
            survey = str(header["FA_SURV"]).rstrip()

        for col in results_assign_columns.keys():
            if col == "DEVICE_TYPE":
                fiber_data[col][0:npos] = "POS"
                fiber_data[col][npos:] = "ETC"
            elif col == "FIBER":
                fiber_data[col][0:npos] = fbassign[col]
                fiber_data[col][npos:] = -1
            else:
                fiber_data[col][0:npos] = fbassign[col]
                if col in fbsky.dtype.names:
                    fiber_data[col][npos:] = fbsky[col]

        full_targets_columns = [(x, y) for x, y in
                                results_targets_columns.items()]
        full_names = [x for x, y in results_targets_columns.items()]
        fbtargets = fd["TARGETS"].read()
        nrawtarget = fd["TARGETS"].get_nrows()

        # For some legacy files, the sky targets assigned to ETC fibers are
        # not listed in the TARGETS HDU.  We check that here.
        copy_etc = False
        ntarget = nrawtarget
        for etctg in fbsky["TARGETID"]:
            if (etctg >= 0) and (etctg not in fbtargets["TARGETID"]):
                ntarget += netc
                copy_etc = True

        for fld in list(fbtargets.dtype.names):
            subd = fbtargets.dtype[fld].subdtype
            if fld not in full_names:
                full_names.append(fld)
                if subd is None:
                    full_targets_columns.extend(
                        [(fld, fbtargets.dtype[fld].str)])
                else:
                    full_targets_columns.extend([(fld, subd[0], subd[1])])

        targets_dtype = np.dtype(full_targets_columns)

        targets_data = np.zeros(ntarget, dtype=targets_dtype)

        fsciencemask = None
        fstdmask = None
        fskymask = None
        fsuppskymask = None
        fexcludemask = None
        fcol = None
        if survey is not None:
            fsciencemask, fstdmask, fskymask, fsuppskymask, fsafemask, \
                fexcludemask = default_survey_target_masks(survey)
            fcol = "FA_TARGET"
        else:
            fsurvey, fcol, fsciencemask, fstdmask, fskymask, \
                fsuppskymask, fsafemask, \
                fexcludemask = default_target_masks(fbtargets)
            survey = fsurvey

        for col in full_names:
            if col == "FA_TYPE":
                targets_data[col][:nrawtarget] = [
                    desi_target_type(x, fsciencemask, fstdmask, fskymask,
                                     fsuppskymask, fsafemask, fexcludemask)
                    for x in fbtargets[fcol][:]]
            elif col == "FA_TARGET":
                targets_data[col][:nrawtarget] = fbtargets[fcol]
            elif col == "TARGET_RA":
                targets_data[col][:nrawtarget] = fbtargets["RA"][:]
            elif col == "TARGET_DEC":
                targets_data[col][:nrawtarget] = fbtargets["DEC"][:]
            else:
                targets_data[col][:nrawtarget] = fbtargets[col][:]
        if copy_etc:
            for col in full_names:
                if col == "FA_TYPE":
                    targets_data[col][nrawtarget:] = [
                        TARGET_TYPE_SKY for x in range(netc)]
                elif col == "FA_TARGET":
                    targets_data[col][nrawtarget:] = [
                        fskymask for x in range(netc)]
                elif col in fbsky.dtype.names:
                    targets_data[col][nrawtarget:] = fbsky[col][:]

        del fbassign
        del fbtargets
        if "POTENTIALTARGETID" in fd:
            log.warning("File {} is an old format with an incompatible "
                        "packing of the potential targets HDU."
                        .format(tile_file))
            avail_data = fd["POTENTIALTARGETID"].read()
        else:
            avail_data = fd["POTENTIAL_ASSIGNMENTS"].read()
    elif "FASSIGN" in fd:
        # We have raw outputs
        header = fd["FASSIGN"].read_header()
        fiber_data = fd["FASSIGN"].read()
        targets_data = fd["FTARGETS"].read()
        avail_data = fd["FAVAIL"].read()
    else:
        msg = "file {} does not contain FIBERASSIGN or FASSIGN HDUs".format(
            tile_file
        )
        log.error(msg)
        raise RuntimeError(msg)

    if "GFA_TARGETS" in fd:
        gfa_targets = fd["GFA_TARGETS"].read()
    else:
        gfa_targets = None

    return header, fiber_data, targets_data, avail_data, gfa_targets


merge_results_tile_tgbuffers = None
merge_results_tile_tgdtypes = None
merge_results_tile_tgshapes = None
merge_results_tile_skybuffers = None
merge_results_tile_skydtypes = None
merge_results_tile_skyshapes = None


def merge_results_tile_initialize(tgbufs, tgdtypes, tgshapes, skybufs,
                                  skydtypes, skyshapes):
    global merge_results_tile_tgbuffers
    global merge_results_tile_tgdtypes
    global merge_results_tile_tgshapes
    merge_results_tile_tgbuffers = tgbufs
    merge_results_tile_tgdtypes = tgdtypes
    merge_results_tile_tgshapes = tgshapes
    global merge_results_tile_skybuffers
    global merge_results_tile_skydtypes
    global merge_results_tile_skyshapes
    merge_results_tile_skybuffers = skybufs
    merge_results_tile_skydtypes = skydtypes
    merge_results_tile_skyshapes = skyshapes
    return


merged_fiberassign_swap = {
    "RA": "TARGET_RA",
    "DEC": "TARGET_DEC",
    "RA_IVAR": "TARGET_RA_IVAR",
    "DEC_IVAR": "TARGET_DEC_IVAR",
    "APFLUX_G": "FIBERFLUX_G",
    "APFLUX_R": "FIBERFLUX_R",
    "APFLUX_Z": "FIBERFLUX_Z",
    "APFLUX_IVAR_G": "FIBERFLUX_IVAR_G",
    "APFLUX_IVAR_R": "FIBERFLUX_IVAR_R",
    "APFLUX_IVAR_Z": "FIBERFLUX_IVAR_Z",
}

merged_fiberassign_req_columns = OrderedDict([
    ("TARGETID", "i8"),
    ("PETAL_LOC", "i2"),
    ("DEVICE_LOC", "i4"),
    ("LOCATION", "i4"),
    ("FIBER", "i4"),
    ("FIBERSTATUS", "i4"),
    ("TARGET_RA", "f8"),
    ("TARGET_DEC", "f8"),
    ("PMRA", "f4"),
    ("PMDEC", "f4"),
    ("PMRA_IVAR", "f4"),
    ("PMDEC_IVAR", "f4"),
    ("REF_EPOCH", "f4"),
    ("LAMBDA_REF", "f4"),
    ("FA_TARGET", "i8"),
    ("FA_TYPE", "u1"),
    ("OBJTYPE", "a3"),
    ("FIBERASSIGN_X", "f4"),
    ("FIBERASSIGN_Y", "f4"),
    ("NUMTARGET", "i2"),
    ("PRIORITY", "i4"),
    ("SUBPRIORITY", "f8"),
    ("OBSCONDITIONS", "i4"),
])

merged_skymon_columns = OrderedDict([
    ("FIBER", "i4"),
    ("LOCATION", "i4"),
    ("NUMTARGET", "i2"),
    ("TARGETID", "i8"),
    ("BRICKID", "i4"),
    ("BRICK_OBJID", "i4"),
    ("FA_TARGET", "i8"),
    ("FA_TYPE", "u1"),
    ("TARGET_RA", "f8"),
    ("TARGET_DEC", "f8"),
    ("FIBERASSIGN_X", "f4"),
    ("FIBERASSIGN_Y", "f4"),
    ("BRICKNAME", "a8"),
    ("FIBERSTATUS", "i4"),
    ("PETAL_LOC", "i2"),
    ("DEVICE_LOC", "i4"),
    ("PRIORITY", "i4"),
    ("SUBPRIORITY", "f8"),
    ("FIBERFLUX_G", "f4"),
    ("FIBERFLUX_R", "f4"),
    ("FIBERFLUX_Z", "f4"),
    ("FIBERFLUX_IVAR_G", "f4"),
    ("FIBERFLUX_IVAR_R", "f4"),
    ("FIBERFLUX_IVAR_Z", "f4"),
])

merged_potential_columns = OrderedDict([
    ("TARGETID", "i8"),
    ("FIBER", "i4"),
    ("LOCATION", "i4")
])


def merge_results_tile(out_dtype, copy_fba, params):
    """Merge results for one tile.

    This uses target catalog data which has been pre-staged to shared memory.
    This function should not be called directly, but only from the
    merge_results() function.

    Args:
        out_dtype (np.dtype):  The output recarray dtype for the merged target
            HDU.  This is the union of columns chosen from the input catalogs.
        copy_fba (bool):  If True, copy the original raw fiberassign HDUs onto
            the end of the output file.
        params (tuple):  The tile ID and input / output files.  Set by
            multiprocessing call.
    Returns:
        None.

    """
    (tile_id, infile, outfile) = params
    log = Logger.get()

    log.info("Reading raw tile data {}".format(infile))

    tm = Timer()
    tm.start()

    inhead, fiber_data, targets_data, avail_data, gfa_targets = \
        read_assignment_fits_tile((tile_id, infile))

    # tm.stop()
    # tm.report("  read input data {}".format(tile_id))
    # tm.clear()
    # tm.start()

    # Get the list of all targets we are considering for this tile.  This
    # will be either just the assigned targets or all available targets
    # depending on how the user wrote the file.

    tile_tgids = np.copy(targets_data["TARGETID"])
    tile_tgindx = {y: x for x, y in enumerate(tile_tgids)}

    # Extract just these targets from the full set of catalogs

    tile_targets_dtype = np.dtype(out_dtype.fields)
    tile_targets = np.zeros(len(tile_tgids), dtype=tile_targets_dtype)

    # Copy original data

    for field in tile_targets_dtype.names:
        if field in results_targets_columns:
            tile_targets[field] = targets_data[field]

    # tm.stop()
    # tm.report("  index input targets {}".format(tile_id))
    # tm.clear()
    # tm.start()

    # Loop over input target files and copy data.  Note: these are guaranteed
    # to be sorted by TARGETID, since that is done during staging to shared
    # memory.
    targetfiles = list(merge_results_tile_tgbuffers.keys())
    for tf in targetfiles:
        tgview = np.frombuffer(merge_results_tile_tgbuffers[tf],
                               dtype=merge_results_tile_tgdtypes[tf])\
                               .reshape(merge_results_tile_tgshapes[tf])
        # Some columns may not exist in all target files (e.g. PRIORITY),
        # So we select the valid columns for this file and only copy those.
        # The ordering of the targets in the catalog is not guaranteed to be
        # sorted, and the output table is sorted by fiber ID, not target.  So
        # must build an explicit mapping from target catalog rows to output
        # table rows.
        tgids = tgview["TARGETID"]
        inrows = np.where(np.isin(tgids, tile_tgids, assume_unique=True))[0]
        outrows = np.where(np.isin(tile_tgids, tgids, assume_unique=True))[0]
        tfcolsin = list()
        tfcolsout = list()
        for c in tgview.dtype.names:
            nm = c
            if c in merged_fiberassign_swap:
                nm = merged_fiberassign_swap[c]
            if nm in out_dtype.names:
                tfcolsin.append(c)
                tfcolsout.append(nm)

        # tm.stop()
        # tm.report("  indexing {} for {}".format(tf, tile_id))
        # tm.clear()
        # tm.start()

        if len(outrows) > 0:
            for irw, orw in zip(inrows, outrows):
                for c, nm in zip(tfcolsin, tfcolsout):
                    tile_targets[nm][orw] = tgview[c][irw]
        del tgids
        del inrows
        del outrows
        del tgview

        # tm.stop()
        # tm.report("  copy targets from {} for {}".format(tf, tile_id))
        # tm.clear()
        # tm.start()

    skytargetfiles = list(merge_results_tile_skybuffers.keys())
    for tf in skytargetfiles:
        skyview = np.frombuffer(merge_results_tile_skybuffers[tf],
                                dtype=merge_results_tile_skydtypes[tf])\
                               .reshape(merge_results_tile_skyshapes[tf])
        # Some columns may not exist in all target files (e.g. PRIORITY),
        # So we select the valid columns for this file and only copy those.
        # The ordering of the targets in the catalog is not guaranteed to be
        # sorted, and the output table is sorted by fiber ID, not target.  So
        # must build an explicit mapping from target catalog rows to output
        # table rows.
        skyids = skyview["TARGETID"]
        inrows = np.where(np.isin(skyids, tile_tgids, assume_unique=True))[0]
        outrows = np.where(np.isin(tile_tgids, skyids, assume_unique=True))[0]
        tfcolsin = list()
        tfcolsout = list()
        for c in skyview.dtype.names:
            nm = c
            if c in merged_fiberassign_swap:
                nm = merged_fiberassign_swap[c]
            if nm in out_dtype.names:
                tfcolsin.append(c)
                tfcolsout.append(nm)

        if len(outrows) > 0:
            for irw, orw in zip(inrows, outrows):
                for c, nm in zip(tfcolsin, tfcolsout):
                    tile_targets[nm][orw] = skyview[c][irw]
        del skyids
        del inrows
        del outrows
        del skyview

    # Now we have a reduced set of target data including only the targets
    # relevant for this tile, and with data merged from all the files.
    # Next, merge this with the assignment information.

    # Determine the rows of the assignment that are for science and sky
    # monitor positioners.
    science_rows = np.where(fiber_data["DEVICE_TYPE"].astype(str) == "POS")[0]
    sky_rows = np.where(fiber_data["DEVICE_TYPE"].astype(str) == "ETC")[0]

    # Construct output recarray
    outdata = np.zeros(len(fiber_data), dtype=out_dtype)

    # Build mapping from FTARGETS rows (sorted by target) to FASSIGN rows
    # (sorted by fiber).

    # Rows containing assigned fibers
    fassign_valid = np.where(fiber_data["TARGETID"] >= 0)[0]

    # Rows in target tables containing assigned targets These indices are
    # valid for both the original input FTARGETS data and also our per-tile
    # copy of the target catalog data.  These indices are essentially random
    # access into the target table.
    target_rows = np.array([tile_tgindx[x] for x in
                            fiber_data["TARGETID"][fassign_valid]])

    # tm.stop()
    # tm.report("  fiber / target index mapping {}".format(tile_id))
    # tm.clear()
    # tm.start()

    # Copy original data and also determine which of our required columns
    # will come from external catalogs
    external_cols = list()
    for field in out_dtype.names:
        if field in results_assign_columns:
            # Copy assignment and fiber property columns directly.
            outdata[field] = fiber_data[field]
        else:
            if field in results_targets_columns:
                # This is a column we are copying from our raw output.
                if (len(target_rows) > 0):
                    outdata[field][fassign_valid] = \
                        targets_data[field][target_rows]
            else:
                # This required column is coming from external catalogs
                external_cols.append(field)

    # tm.stop()
    # tm.report("  copy raw data to output {}".format(tile_id))
    # tm.clear()
    # tm.start()

    # Now copy external target properties for the remaining columns.
    if (len(target_rows) > 0):
        # # Looping over rows and then columns is faster, likely because there
        # # are so many columns and the data is stored row-major.
        # for irw, orw in zip(target_rows, fassign_valid):
        #     for c in external_cols:
        #         realname = c
        #         if c in merged_fiberassign_swap:
        #             realname = merged_fiberassign_swap[c]
        #         outdata[realname][orw] = tile_targets[c][irw]
        for c in external_cols:
            if c == "OBJTYPE":
                # FIXME:  This column is redundant to doing a simple bitwise
                # operation on the target bit field, and it is specific to
                # the main survey.  Personally I believe it should be
                # removed completely.  -TK
                objtype = np.zeros(len(fassign_valid), dtype="S3")
                if "DESI_TARGET" in out_dtype.names:
                    # This is a main survey file.
                    objtype[:] = "TGT"
                    is_sky = (tile_targets["DESI_TARGET"][target_rows]
                              & desi_mask.SKY) != 0
                    objtype[is_sky] = "SKY"
                    badmask = \
                        desi_mask.mask("BAD_SKY|NO_TARGET|IN_BRIGHT_OBJECT")
                    is_bad = (tile_targets["DESI_TARGET"][target_rows]
                              & badmask) != 0
                    objtype[is_bad] = "BAD"
                else:
                    # This is some other survey
                    objtype = ["NA" for x in range(len(fassign_valid))]
                outdata[c][fassign_valid] = objtype
            else:
                outdata[c][fassign_valid] = tile_targets[c][target_rows]

    # Special check for REF_EPOCH which isn't in MTL yet
    ismoving = (outdata['PMRA'] != 0) | (outdata['PMDEC'] != 0)
    if np.all(outdata['REF_EPOCH'] == 0.0) and np.any(ismoving):
        log.error('REF_EPOCH not set, using 2015.5')
        outdata['REF_EPOCH'][ismoving] = 2015.5

    # tm.stop()
    # tm.report("  copy external data to output {}".format(tile_id))
    # tm.clear()
    # tm.start()

    cols_to_keep = list()
    for c in outdata.dtype.names:
        if len(outdata[c].shape) < 2: #Don't propagate 2D target columns into FIBERASSIGN HDU
            cols_to_keep.append(c)
    outdata = outdata[cols_to_keep]

    cols_to_keep = list()
    for c in tile_targets.dtype.names:
        if len(tile_targets[c].shape) < 2: #Don't propagate 2D target columns into TARGETS HDU
            cols_to_keep.append(c)
    tile_targets = tile_targets[cols_to_keep]

    # Create the file
    if os.path.isfile(outfile):
        os.remove(outfile)
    fd = fitsio.FITS(outfile, "rw")

    # Write a heaader-only primary HDU
    fd.write(None, header=inhead, extname="PRIMARY")

    # Write the main FIBERASSIGN HDU- only the data for science positioners
    # Enforce sorting by FIBER
    log.info("Writing new data {}".format(outfile))
    ii = np.argsort(outdata['FIBER'][science_rows])
    fd.write(outdata[science_rows][ii], header=inhead, extname="FIBERASSIGN")

    # Now write out the sky monitor fibers.  We extract the rows and columns
    # from the already-computed recarray.

    skymon_dtype = np.dtype([(x, y) for x, y in merged_skymon_columns.items()])
    skymon = np.zeros(len(sky_rows), dtype=skymon_dtype)
    # Sky monitor fake FIBER column with values 0-19.  The fake FIBER value in
    # the raw data is already based on increasing LOCATION value.
    skymon_fiber = np.arange(len(sky_rows), dtype=np.int32)
    for field in skymon_dtype.names:
        if field == "FIBER":
            skymon["FIBER"] = skymon_fiber
        else:
            skymon[field] = outdata[field][sky_rows]
    fd.write(skymon, header=inhead, extname="SKY_MONITOR")

    # Copy GFA data if it exists
    if gfa_targets is not None:
        fd.write(gfa_targets, header=inhead, extname="GFA_TARGETS")

    # Write the per-tile catalog information also.  Sadly, this HDU is
    # expected to have the original column names.  We swap them back, which
    # is fine since we are going to delete this data after writing anyway.
    backswap = {y: x for x, y in merged_fiberassign_swap.items()}
    curnames = np.copy(tile_targets.dtype.names)
    newnames = list()
    for nm in curnames:
        if nm in backswap:
            newnames.append(backswap[nm])
        else:
            newnames.append(nm)
    tile_targets.dtype.names = newnames
    fd.write(tile_targets, header=inhead, extname="TARGETS")
    del tile_targets

    # Write the "POTENTIAL_ASSIGNMENTS" HDU
    potential_dtype = np.dtype([(x, y) for x, y
                                in merged_potential_columns.items()])
    potential = np.zeros(len(avail_data), dtype=potential_dtype)

    locfiber = {x: y for x, y in
                zip(fiber_data["LOCATION"], fiber_data["FIBER"])}
    potential["LOCATION"] = avail_data["LOCATION"]
    potential["TARGETID"] = avail_data["TARGETID"]
    potential["FIBER"] = [locfiber[x] for x in avail_data["LOCATION"]]
    fd.write(potential, header=inhead, extname="POTENTIAL_ASSIGNMENTS")

    # Now copy the original HDUs
    if copy_fba:
        fd.write(fiber_data, header=inhead, extname="FASSIGN")
        fd.write(targets_data, header=inhead, extname="FTARGETS")
        fd.write(avail_data, header=inhead, extname="FAVAIL")

    # Close the file
    fd.close()

    # tm.stop()
    # tm.report("  write data to file {}".format(tile_id))
    # tm.clear()
    # tm.start()

    # Try to encourage python to free some memory...
    del fd

    del avail_data
    del targets_data
    del fiber_data
    return


def merge_results(targetfiles, skyfiles, tiles, result_dir=".",
                  result_prefix="fba-", result_split_dir=False,
                  out_dir=None, out_prefix="fiberassign-", out_split_dir=False,
                  columns=None, copy_fba=True):
    """Merge target files and assignment output.

    Full target data is stored in shared memory and then multiple processes
    copy this data into the per-tile files.

    Args:
        targetfiles (list):  List of pathnames containing the original input
            target files.  The rows of the file MUST be sorted by TARGETID.
        skyfiles (list):  List of pathnames containing the original input
            sky target files.
        tiles (list):  List of tile IDs to process.
        result_dir (str):  Top-level directory of fiberassign results.
        result_prefix (str):  Prefix of each per-tile file name.
        result_split_dir (bool):  Results are in split tile directories.
        out_dir (str):  Top-level directory for merged outputs.
        out_prefix (str):  Prefix of each per-tile output file name.
        out_split_dir (bool):  Write outputs in split tile directories.
        columns (list):  List of column names to propagate from the input
            target files (default is all).
        copy_fba (bool):  If True, propagate the original raw HDUs at the end
            of each merged output file.

    Returns:
        None.

    """
    # Load the full set of target files into memory.  Also build a mapping of
    # target ID to row index.  We assume that the result columns have the same
    # dtype in any of the target files.  We take the first target file and
    # construct the output recarray dtype from the columns in that file.
    out_dtype = None
    dcols = [(x, y) for x, y in merged_fiberassign_req_columns.items()]
    dcolnames = [x for x in merged_fiberassign_req_columns.keys()]

    tgdata = dict()
    tgdtype = dict()
    tgshape = dict()
    tghead = dict()

    skydata = dict()
    skydtype = dict()
    skyshape = dict()
    skyhead = dict()

    survey = None

    for tf in targetfiles:
        tm = Timer()
        tm.start()
        fd = fitsio.FITS(tf)
        tghead[tf] = fd[1].read_header()
        # Allocate a shared memory buffer for the target data
        tglen = fd[1].get_nrows()
        tgshape[tf] = (tglen,)
        tgdtype[tf], tempoff, tempisvararray = fd[1].get_rec_dtype()
        tgbytes = tglen * tgdtype[tf].itemsize
        tgdata[tf] = RawArray("B", tgbytes)
        tgview = np.frombuffer(tgdata[tf],
                               dtype=tgdtype[tf]).reshape(tgshape[tf])
        # Read data directly into shared buffer
        tgview[:] = fd[1].read()
        if survey is None:
            (survey, col, sciencemask, stdmask, skymask, suppskymask,
             safemask, excludemask) = default_target_masks(tgview)

        # Sort rows by TARGETID if not already done
        tgviewids = tgview["TARGETID"]
        if not np.all(tgviewids[:-1] <= tgviewids[1:]):
            tgview.sort(order="TARGETID", kind="heapsort")

        tm.stop()
        tm.report("Read {} into shared memory".format(tf))

        # Add any missing columns to our output dtype record format.
        tfcols = list(tgview.dtype.names)
        if columns is not None:
            tfcols = [x for x in tfcols if x in columns]
        for col in tfcols:
            subd = tgview.dtype[col].subdtype
            colname = col
            if col in merged_fiberassign_swap:
                colname = merged_fiberassign_swap[col]
            if colname not in dcolnames:
                if subd is None:
                    dcols.extend([(colname, tgview.dtype[col].str)])
                else:
                    dcols.extend([(colname, subd[0], subd[1])])
                dcolnames.append(colname)

    for tf in skyfiles:
        tm = Timer()
        tm.start()
        fd = fitsio.FITS(tf)
        skyhead[tf] = fd[1].read_header()
        # Allocate a shared memory buffer for the target data
        skylen = fd[1].get_nrows()
        skyshape[tf] = (skylen,)
        skydtype[tf], tempoff, tempisvararray = fd[1].get_rec_dtype()
        skybytes = skylen * skydtype[tf].itemsize
        skydata[tf] = RawArray("B", skybytes)
        skyview = np.frombuffer(skydata[tf],
                                dtype=skydtype[tf]).reshape(skyshape[tf])
        # Read data directly into shared buffer
        skyview[:] = fd[1].read()

        # Sort rows by TARGETID if not already done
        skyviewids = skyview["TARGETID"]
        if not np.all(skyviewids[:-1] <= skyviewids[1:]):
            skyview.sort(order="TARGETID", kind="heapsort")

        tm.stop()
        tm.report("Read {} into shared memory".format(tf))

        # Add any missing columns to our output dtype record format.
        tfcols = list(skyview.dtype.names)
        if columns is not None:
            tfcols = [x for x in tfcols if x in columns]
        for col in tfcols:
            subd = skyview.dtype[col].subdtype
            colname = col
            if col in merged_fiberassign_swap:
                colname = merged_fiberassign_swap[col]
            if colname not in dcolnames:
                if subd is None:
                    dcols.extend([(colname, skyview.dtype[col].str)])
                else:
                    dcols.extend([(colname, subd[0], subd[1])])
                dcolnames.append(colname)

    out_dtype = np.dtype(dcols)

    # For each tile, find the target IDs used.  Construct the output recarray
    # and copy data into place.

    merge_tile = partial(merge_results_tile, out_dtype, copy_fba)

    if out_dir is None:
        out_dir = result_dir

    tile_map_list = [(x, result_path(x, dir=result_dir, prefix=result_prefix,
                                     split=result_split_dir),
                      result_path(x, dir=out_dir, prefix=out_prefix,
                                  create=True, split=out_split_dir))
                     for x in tiles]

    with mp.Pool(processes=default_mp_proc,
                 initializer=merge_results_tile_initialize,
                 initargs=(tgdata, tgdtype, tgshape, skydata,
                           skydtype, skyshape)) as pool:
        results = pool.map(merge_tile, tile_map_list)

    return

def run(
    asgn,
    std_per_petal=10,
    sky_per_petal=40,
    start_tile=-1,
    stop_tile=-1,
    redistribute=True
):
    """Run fiber assignment.

    Given an already-constructed Assignment class instance, run the assignment in a
    standard way.

    This is designed to be the main "driver" function used by higher-level code or
    commandline tools.  The purpose of this function is to ensure that all of those
    tools are assigning targets in the same way.

    Args:
        asgn (Assignment):  The assignment class
        std_per_petal (int):  The number of standards to assign per petal
        sky_per_petal (int):  The number of sky to assign per petal
        start_tile (int):  If specified, the first tile ID to assign.
        stop_tile (int):  If specified, the last tile ID to assign.
        redistribute (bool):  If True, attempt to shift science targets to unassigned
            fibers on later tiles in order to balance the number per petal.

    Returns:
        None

    """
    gt = GlobalTimers.get()

    # First-pass assignment of science targets
    gt.start("Assign unused fibers to science targets")
    asgn.assign_unused(TARGET_TYPE_SCIENCE, -1, "POS", start_tile, stop_tile)
    gt.stop("Assign unused fibers to science targets")

    # Redistribute science targets across available petals
    if redistribute:
        gt.start("Redistribute science targets")
        asgn.redistribute_science(start_tile, stop_tile)
        gt.stop("Redistribute science targets")

    # Assign standards, up to some limit
    gt.start("Assign unused fibers to standards")
    asgn.assign_unused(
        TARGET_TYPE_STANDARD, std_per_petal, "POS", start_tile, stop_tile
    )
    gt.stop("Assign unused fibers to standards")

    # Assign sky to unused fibers, up to some limit
    gt.start("Assign unused fibers to sky")
    asgn.assign_unused(
        TARGET_TYPE_SKY, sky_per_petal, "POS", start_tile, stop_tile
    )
    gt.stop("Assign unused fibers to sky")

    # Assign suppsky to unused fibers, up to some limit
    gt.start("Assign unused fibers to supp_sky")
    asgn.assign_unused(
        TARGET_TYPE_SUPPSKY, sky_per_petal, "POS", start_tile, stop_tile
    )
    gt.stop("Assign unused fibers to supp_sky")

    # Force assignment if needed
    gt.start("Force assignment of sufficient standards")
    asgn.assign_force(
        TARGET_TYPE_STANDARD, std_per_petal, start_tile, stop_tile
    )
    gt.stop("Force assignment of sufficient standards")

    gt.start("Force assignment of sufficient sky")
    asgn.assign_force(TARGET_TYPE_SKY, sky_per_petal, start_tile, stop_tile)
    gt.stop("Force assignment of sufficient sky")

    gt.start("Force assignment of sufficient supp_sky")
    asgn.assign_force(TARGET_TYPE_SUPPSKY, sky_per_petal, start_tile, stop_tile)
    gt.stop("Force assignment of sufficient supp_sky")

    # If there are any unassigned fibers, try to place them somewhere.
    # Assigning science again is a no-op, but...
    gt.start("Assign remaining unassigned fibers")
    asgn.assign_unused(TARGET_TYPE_SCIENCE, -1, "POS", start_tile, stop_tile)
    asgn.assign_unused(TARGET_TYPE_STANDARD, -1, "POS", start_tile, stop_tile)
    asgn.assign_unused(TARGET_TYPE_SKY, -1, "POS", start_tile, stop_tile)
    asgn.assign_unused(TARGET_TYPE_SUPPSKY, -1, "POS", start_tile, stop_tile)

    # Assign safe location to unused fibers (no maximum).  There should
    # always be at least one safe location (i.e. "BAD_SKY") for each fiber.
    # So after this is run every fiber should be assigned to something.
    asgn.assign_unused(TARGET_TYPE_SAFE, -1, "POS", start_tile, stop_tile)
    gt.stop("Assign remaining unassigned fibers")

    # Assign sky monitor fibers
    gt.start("Assign sky monitor fibers")
    asgn.assign_unused(TARGET_TYPE_SKY, -1, "ETC", start_tile, stop_tile)
    asgn.assign_unused(TARGET_TYPE_SUPPSKY, -1, "ETC", start_tile, stop_tile)
    asgn.assign_unused(TARGET_TYPE_SAFE, -1, "ETC", start_tile, stop_tile)
    gt.stop("Assign sky monitor fibers")

    return asgn
