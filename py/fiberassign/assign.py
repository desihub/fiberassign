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

from astropy import units
from astropy.coordinates import SkyCoord

import multiprocessing as mp
from multiprocessing.sharedctypes import RawArray

from functools import partial

from collections import OrderedDict
from types import SimpleNamespace

import fitsio

from desiutil.depend import add_dependencies, setdep
from desiutil.dust import SFDMap

import desimodel.focalplane

from desitarget.targetmask import desi_mask

from ._version import __version__

from .utils import Logger, Timer, default_mp_proc, GlobalTimers

from .targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY, TARGET_TYPE_SUPPSKY,
                      TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE, desi_target_type,
                      default_target_masks, default_survey_target_masks)

from .hardware import (FIBER_STATE_UNASSIGNED, FIBER_STATE_STUCK,
                       FIBER_STATE_BROKEN, FIBER_STATE_RESTRICT,
                       radec2xy, xy2radec, xy2cs5)

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
# AR these columns are for the fba-TILEID.fits files

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
    ("PLATE_RA", "f8"),
    ("PLATE_DEC", "f8"),
    #("PLATE_REF_EPOCH", "f8"),
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
    ("PLATE_RA", "f8"),
    ("PLATE_DEC", "f8"),
    #("PLATE_REF_EPOCH", "f8"),
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


def write_assignment_fits_tile(asgn, tagalong, fulltarget, overwrite, params):
    """Write a single tile assignment to a FITS file.

    Args:
        outroot (str):  full path of the output root file name.
        asgn (Assignment):  the assignment class instance.
        tagalong (TargetTagalong): a data structure that carries RA,Dec, and other information for targets from the targeting files, to be written to the fiberassign outputs.
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
        tile_file, gfa_targets, stuck_sky_tile, tile_xy_cs5 = params
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

        tg_indx = dict((t,i) for i,t in enumerate(tgids))
        tg_type, tg_priority, tg_subpriority = tgs.get_type_priority_subpriority(tgids)

        if tile_xy_cs5 is not None:
            (tg_obscond,) = tagalong.get_for_ids(tgids, ['OBSCOND'])

            xy = np.array([tile_xy_cs5[tid] for tid in tgids])
            tg_x = xy[:,0]
            tg_y = xy[:,1]
            del xy
        else:
            (tg_ra, tg_dec, tg_obscond) = tagalong.get_for_ids(
                tgids, ['RA', 'DEC', 'OBSCOND'])

            # We compute the X / Y focalplane coordinates for ALL available
            # targets, not just the assigned ones.  This allows us to write out
            # the focalplane coordinates of all available targets as computed by
            # the code at the time it was run.
            #
            # NOTE:  The output format is explicitly CS5 coordinates, even though
            # we use curved focal surface internally.
            tg_x,tg_y = radec2xy(
                hw, tile_ra, tile_dec, tile_obstime, tile_obstheta, tile_obsha,
                tg_ra, tg_dec, True)
            del tg_ra, tg_dec

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

        margins = hw.added_margins
        keys = list(margins.keys())
        keys.sort()
        for k in keys:
            header["FA_M_%s" % k[:3]] = margins[k]

        header["REQRA"] = tile_ra
        header["REQDEC"] = tile_dec
        header["FIELDNUM"] = 0
        header["FA_VER"] = __version__
        header["FA_SURV"] = tgs.survey()

        #- Add code dependency versions for default list of packages, then
        #- call again for ones that are most critical to record just in
        #- case they get dropped from the default list in the future.
        #- (it won't add a second copy if they are already there)
        add_dependencies(header)
        add_dependencies(
            header,
            module_names=[
                "numpy",
                "matplotlib",
                "astropy",
                "fitsio",
                "desiutil",
                "desimodel",
                "desitarget",
                "desimeter",
            ]
        )

        #- Keep SKYBRICKS_DIR used to lookup sky locations,
        #- shortening full path if possible
        # AR add SKYHEALPIX_DIR similarly, now looping over the variable names
        for skyname in ["SKYBRICKS_DIR", "SKYHEALPIXS_DIR"]:
            skydir = os.getenv(skyname, None)
            if (skydir is not None) and ("DESI_ROOT" in os.environ):
                if skydir.startswith(os.environ["DESI_ROOT"]):
                    skydir = skydir.replace(
                            os.environ["DESI_ROOT"], "$DESI_ROOT", 1)
            setdep(header, skyname, skydir)

        fd.write(None, header=header, extname="PRIMARY")

        # FIXME:  write "-1" for unassigned targets.  Write all other fiber
        # status and other properties to this table.

        log.debug("Write:  copying assignment data for tile {}"
                  .format(tile_id))

        fdata = np.zeros(nloc, dtype=assign_dtype)
        fdata["LOCATION"] = locs

        # For unassigned fibers, we give each location a unique negative
        # number based on the tile and loc.
        unassign_offset = tile_id * 10000
        assigned_tgids = np.array([tdata[x] if x in tdata.keys()
                                  else -(unassign_offset + x)
                                  for x in locs], dtype=np.int64)
        fdata["TARGETID"] = assigned_tgids
        # Fill in all the data values we've been carrying around!
        tagalong.set_data(assigned_tgids, fdata)

        # Rows containing assigned locations
        assigned_valid = np.where(assigned_tgids >= 0)[0]
        assigned_invalid = np.where(assigned_tgids < 0)[0]

        # Buffers for X/Y/RA/DEC
        assigned_tgx = np.full(nloc, 9999.9, dtype=np.float64)
        assigned_tgy = np.full(nloc, 9999.9, dtype=np.float64)
        assigned_tgtype = np.zeros(nloc, dtype=np.uint8)

        # For later, locs of stuck fibers that land on good sky
        stuck_sky_locs = set()

        if (len(assigned_invalid) > 0):
            # Fill our unassigned location X/Y coordinates with either their fixed
            # location (for positioners that are nonfunctional) or as "folded in"
            # as possible.  Then convert these to RA/DEC.
            # NOTE:  Positioner locations are in curved focal surface coordinates.
            empty_fibers = locs[assigned_invalid]

            empty_x = np.zeros(len(empty_fibers), dtype=np.float64)
            empty_y = np.zeros(len(empty_fibers), dtype=np.float64)
            empty_tgtype = np.zeros(len(empty_fibers), dtype=np.uint8)

            # pull relevant fields out of hw into python
            hwduck = SimpleNamespace()
            for attribute in [
                    'state', 'loc_pos_curved_mm', 'loc_theta_pos',
                    'loc_theta_offset', 'loc_phi_pos', 'loc_phi_offset',
                    'loc_theta_arm', 'loc_phi_arm', 'loc_theta_min',
                    'loc_phi_min', 'loc_theta_max', 'loc_phi_max']:
                setattr(hwduck, attribute, getattr(hw, attribute))

            for iloc, loc in enumerate(empty_fibers):
                # For stuck positioners on good sky, set the FA_TYPE SKY bit.
                if (hwduck.state[loc] & FIBER_STATE_STUCK) and stuck_sky_tile is not None:
                    # (note, the stuck_sky code does the check for STUCK, not BROKEN,
                    #  type POS, etc; that's in the stuck_sky_tile map.)
                    if stuck_sky_tile.get(loc, False):
                        empty_tgtype[iloc] = TARGET_TYPE_SKY
                        stuck_sky_locs.add(loc)

                if (
                    (hwduck.state[loc] & FIBER_STATE_STUCK) or
                    (hwduck.state[loc] & FIBER_STATE_BROKEN)
                ):
                    # This positioner is not moveable and therefore has a theta / phi
                    # position in the hardware model.  Used this fixed fiber location.
                    xy = hw.thetaphi_to_xy(
                        hwduck.loc_pos_curved_mm[loc],
                        hwduck.loc_theta_pos[loc] + hwduck.loc_theta_offset[loc],
                        hwduck.loc_phi_pos  [loc] + hwduck.loc_phi_offset  [loc],
                        hwduck.loc_theta_arm[loc],
                        hwduck.loc_phi_arm[loc],
                        hwduck.loc_theta_offset[loc],
                        hwduck.loc_phi_offset[loc],
                        hwduck.loc_theta_min[loc],
                        hwduck.loc_phi_min[loc],
                        hwduck.loc_theta_max[loc],
                        hwduck.loc_phi_max[loc],
                        ignore_range=True
                    )
                    empty_x[iloc] = xy[0]
                    empty_y[iloc] = xy[1]
                else:
                    # This is an unassigned, working positioner.
                    # Place it in a folded state.
                    theta,phi = get_parked_thetaphi(hwduck.loc_theta_offset[loc],
                                                    hwduck.loc_theta_min[loc],
                                                    hwduck.loc_theta_max[loc],
                                                    hwduck.loc_phi_offset[loc],
                                                    hwduck.loc_phi_min[loc],
                                                    hwduck.loc_phi_max[loc])
                    xy = hw.thetaphi_to_xy(
                        hwduck.loc_pos_curved_mm[loc],
                        theta,
                        phi,
                        hwduck.loc_theta_arm[loc],
                        hwduck.loc_phi_arm[loc],
                        hwduck.loc_theta_offset[loc],
                        hwduck.loc_phi_offset[loc],
                        hwduck.loc_theta_min[loc],
                        hwduck.loc_phi_min[loc],
                        hwduck.loc_theta_max[loc],
                        hwduck.loc_phi_max[loc]
                    )
                    empty_x[iloc] = xy[0]
                    empty_y[iloc] = xy[1]
                    if xy[0] is None:
                        print('WARNING: X,Y position for unassigned, parked positioner loc =', loc, ' is invalid.')

            ra,dec = xy2radec(
                hw, tile_ra, tile_dec, tile_obstime, tile_obstheta, tile_obsha,
                empty_x, empty_y, False, 0)
            # These rows will have had default values; set them now!
            fdata["PLATE_RA" ][assigned_invalid] = ra
            fdata["PLATE_DEC"][assigned_invalid] = dec
            fdata["TARGET_RA" ][assigned_invalid] = ra
            fdata["TARGET_DEC"][assigned_invalid] = dec

            ex, ey = xy2cs5(empty_x, empty_y)
            assigned_tgx[assigned_invalid] = ex
            assigned_tgy[assigned_invalid] = ey
            assigned_tgtype[assigned_invalid] |= empty_tgtype

        if (len(assigned_valid) > 0):
            # The target IDs assigned to fibers (note- NOT sorted)
            assigned_real = np.copy(assigned_tgids[assigned_valid])

            # Mapping of our assigned target rows into the target properties.
            target_rows = [tg_indx[x] for x in assigned_real]

            # Copy values out of the full target list.
            assigned_tgx[assigned_valid] = np.array(tg_x[target_rows])
            assigned_tgy[assigned_valid] = np.array(tg_y[target_rows])
            assigned_tgtype[assigned_valid] = np.array(tg_type[target_rows])

        fdata["FIBERASSIGN_X"] = assigned_tgx
        fdata["FIBERASSIGN_Y"] = assigned_tgy
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
        fstatus |= [0 if (x in tdata.keys()) or (x in stuck_sky_locs)
                    else FIBER_STATE_UNASSIGNED
                    for x in locs]
        # Set stuck / broken / restricted bits.
        fstatus |= [2 if (fstate[x] & FIBER_STATE_STUCK) else 0
                    for x in locs]
        fstatus |= [4 if (fstate[x] & FIBER_STATE_BROKEN) else 0
                    for x in locs]
        fstatus |= [8 if (fstate[x] & FIBER_STATE_RESTRICT) else 0
                    for x in locs]
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
        # Fill in all the data values we've been carrying around!
        tagalong.set_data(tgids, fdata)
        fdata["TARGETID"] = tgids
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
        # The "FAVAIL" (available targets) HDU is sorted first by LOCATION,
        # then by TARGETID.
        for lid in sorted(avail.keys()):
            # lid (location id) is a scalar, tg (target ids) is an array
            tg = avail[lid]
            fdata['LOCATION'][off:off+len(tg)] = lid
            fdata['FIBER']   [off:off+len(tg)] = fibers[lid]
            fdata['TARGETID'][off:off+len(tg)] = sorted(tg)
            off += len(tg)

        # tm.stop()
        # tm.report("  copy avail data tile {}".format(tile_id))
        # tm.clear()
        # tm.start()

        fd.write(fdata, header=header, extname="FAVAIL")
        del fdata

        if gfa_targets is not None:
            # Sort by TARGETID
            gfa_targets = gfa_targets[np.argsort(gfa_targets['TARGETID'])]
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


def write_assignment_fits(tiles, tagalong, asgn, out_dir=".", out_prefix="fba-",
                          split_dir=False, all_targets=False,
                          gfa_targets=None, overwrite=False, stucksky=None,
                          tile_xy_cs5=None):
    """Write out assignment results in FITS format.

    For each tile, all available targets (not only the assigned targets) and
    their properties are written to the first HDU.  The second HDU contains
    one row for every location and target available to that location.

    Args:
        tiles (Tiles):  The Tiles object containing the properties of each
            tile.
        tagalong (TargetTagalong): a data structure that carries RA,Dec, and other information for targets from the targeting files, to be written to the fiberassign outputs.
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
                         asgn, tagalong, all_targets, overwrite)

    for i, tid in enumerate(tileids):
        tra = tilera[tileorder[tid]]
        tdec = tiledec[tileorder[tid]]
        ttheta = tiletheta[tileorder[tid]]
        ttime = tiletime[tileorder[tid]]
        tha = tileha[tileorder[tid]]

        outfile = result_path(tid, dir=out_dir, prefix=out_prefix,
                              create=True, split=split_dir)
        gfa = None
        if gfa_targets is not None:
            gfa = gfa_targets[i]

        stuck = None
        if stucksky is not None:
            stuck = stucksky.get(tid,{})

        txy = None
        if tile_xy_cs5 is not None:
            txy = tile_xy_cs5.get(tid, None)

        params = (tid, tra, tdec, ttheta, ttime, tha, outfile, gfa, stuck, txy)
        write_tile(params)

    tm.stop()
    tm.report("Write output files")

    return


def write_assignment_ascii(tiles, asgn, tagalong, out_dir=".", out_prefix="fba-",
                           split_dir=False):
    """Write out assignment results in ASCII format.

    For each tile, only the final assignment to each tile is written out.  For
    the full information, including available targets, use the FITS format.

    Args:
        tiles (Tiles):  The Tiles object containing the properties of each
            tile.
        asgn (Assignment):  The assignment object.
        tagalong (TargetTagalong): a data structure that carries RA,Dec, and other information for targets from the targeting files, to be written to the fiberassign outputs.
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

                lids = sorted(tdata.keys())
                targetids = [tdata[lid] for lid in lids]
                (ra, dec, obscond) = tagalong.get_for_ids(targetids,
                                              ['RA', 'DEC', 'OBSCOND'])
                for i,lid in enumerate(lids):
                    tgid = tdata[lid]
                    tg = tgs.get(tgid)
                    f.write("{:d} {:d} {:.6f} {:.6f}\n"
                            .format(lid, tgid, ra[i], dec[i], tg.priority,
                                    tg.subpriority, obscond[i], tg.type))
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

def read_assignment_fits_tile(tile_file):
    """Read in results.

    This reads in only the result information that was originally written by
    the fiber assignment.  This function is used internally for reading data
    for use in plotting and for reading the original data for merging with
    input target catalogs.

    Args:
        tile_file: pathname

    Returns:
        (tuple):  The FITS header, assignment recarray, target property
            recarray, the available targets recarray, and the GFA targets

    """
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
                                merged_targets_columns.items()]
        full_names = [x for x, y in merged_targets_columns.items()]
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
        fgaia_stdmask = None
        fskymask = None
        fsuppskymask = None
        fexcludemask = None
        fcol = None
        if survey is not None:
            fsciencemask, fstdmask, fskymask, fsuppskymask, fsafemask, \
                fexcludemask, fgaia_stdmask = default_survey_target_masks(survey)
            fcol = "FA_TARGET"
        else:
            fsurvey, fcol, fsciencemask, fstdmask, fskymask, \
                fsuppskymask, fsafemask, \
                fexcludemask, fgaia_stdmask = default_target_masks(fbtargets)
            survey = fsurvey

        for col in full_names:
            if col == "FA_TYPE":
                # AR protect case where *MWS_TARGET does not exist
                mws_targets = 0 * fbtargets[fcol][:]
                if fcol.replace("DESI", "MWS") in fbtargets.dtype.names:
                    mws_targets = fbtargets[fcol.replace("DESI", "MWS")][:]
                targets_data[col][:nrawtarget] = [
                    desi_target_type(x, y, fsciencemask, fstdmask, fskymask,
                                     fsuppskymask, fsafemask, fexcludemask,
                                     fgaia_stdmask)
                    for x, y in zip(fbtargets[fcol][:], mws_targets)]
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

# minimal set of columns to read from the target file
# AR we also propagate any *_TARGET existing column
minimal_target_columns = OrderedDict([
    ('RELEASE', '>i2'),
    ('BRICKNAME', 'S8'),
    ('BRICKID', '>i4') ,
    ('BRICK_OBJID', '>i4'),
    ('MORPHTYPE', 'S4'),
    ('RA', '>f8'),
    ('DEC', '>f8'),
    ('EBV', '>f4'),
    ('FLUX_G', '>f4'),
    ('FLUX_R', '>f4'),
    ('FLUX_Z', '>f4'),
    ('FLUX_W1', '>f4'), # AR for QA plots
    ('FLUX_W2', '>f4'), # AR for QA plots
    ('FLUX_IVAR_G', '>f4'),
    ('FLUX_IVAR_R', '>f4'),
    ('FLUX_IVAR_Z', '>f4'),
    ('FLUX_IVAR_W1', '>f4'), # https://github.com/desihub/fiberassign/issues/300
    ('FLUX_IVAR_W2', '>f4'), # https://github.com/desihub/fiberassign/issues/300
    ('FIBERFLUX_G', '>f4'),
    ('FIBERFLUX_R', '>f4'),
    ('FIBERFLUX_Z', '>f4'),
    ('FIBERTOTFLUX_G', '>f4'),
    ('FIBERTOTFLUX_R', '>f4'),
    ('FIBERTOTFLUX_Z', '>f4'),
    ('REF_EPOCH', '>f4'),
    ('MASKBITS', '>i2'),
    ('SERSIC', '>f4'),
    ('SHAPE_R', '>f4'),
    ('SHAPE_E1', '>f4'),
    ('SHAPE_E2', '>f4'),
    ('REF_ID', '>i8'),
    ('REF_CAT', 'S2'),
    ('GAIA_PHOT_G_MEAN_MAG', '>f4'),
    ('GAIA_PHOT_BP_MEAN_MAG', '>f4'),
    ('GAIA_PHOT_RP_MEAN_MAG', '>f4'),
    ('PARALLAX', '>f4'),
    ('PMRA', '>f4'),
    ('PMDEC', '>f4'),
    ('PHOTSYS', 'S1'),
    ('TARGETID', '>i8'),
    #('DESI_TARGET', '>i8'),
    #('BGS_TARGET', '>i8'),
    #('MWS_TARGET', '>i8'),
    #('CMX_TARGET', '>i8'),
    #('SV1_DESI_TARGET', '>i8'),
    #('SV1_BGS_TARGET', '>i8'),
    #('SV1_MWS_TARGET', '>i8'),
    #('SCND_TARGET', '>i8'),
    ('SUBPRIORITY', '>f8'),
    ('OBSCONDITIONS', '>i8'),
    ('PRIORITY_INIT', '>i8'),
    ('NUMOBS_INIT', '>i8')
])

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

# AR commenting out NUMTARGET ([desi-survey 1032])
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
    ("REF_EPOCH", "f4"),
    ("LAMBDA_REF", "f4"),
    ("FA_TARGET", "i8"),
    ("FA_TYPE", "u1"),
    ("OBJTYPE", "a3"),
    ("FIBERASSIGN_X", "f4"),
    ("FIBERASSIGN_Y", "f4"),
    #("NUMTARGET", "i2"),
    ("PRIORITY", "i4"),
    ("SUBPRIORITY", "f8"),
    ("OBSCONDITIONS", "i4"),
])
# Columns that should appear at the end of the table.
merged_fiberassign_req_columns_at_end = OrderedDict([
    ("PLATE_RA", "f8"),
    ("PLATE_DEC", "f8"),
])

# AR commenting out NUMTARGET ([desi-survey 1032])
merged_skymon_columns = OrderedDict([
    ("FIBER", "i4"),
    ("LOCATION", "i4"),
    #("NUMTARGET", "i2"),
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
    #("FIBERFLUX_IVAR_G", "f4"),
    #("FIBERFLUX_IVAR_R", "f4"),
    #("FIBERFLUX_IVAR_Z", "f4"),
])

# AR TARGETS extension columns
# AR should be a subsample of merged_fiberassign_req_columns
# AR modulo the name swapping (merged_fiberassign_swap)
# AR we also propagate any *_TARGET existing column
merged_targets_columns = OrderedDict([
    ("TARGETID", "i8"),
    ("RA", "f8"),
    ("DEC", "f8"),
    ("FA_TARGET", "i8"),
    ("FA_TYPE", "u1"),
    ("PRIORITY", "i4"),
    ("SUBPRIORITY", "f8"),
    ("OBSCONDITIONS", "i4"),
    #("PLATE_RA", "f8"),
    #("PLATE_DEC", "f8"),
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

    Notes:
        20210928 : we populate EBV=0 rows in FIBERASSIGN with correct values.

    """
    (tile_id, infile, outfile) = params
    log = Logger.get()

    log.info("Reading raw tile data {}".format(infile))

    tm = Timer()
    tm.start()

    inhead, fiber_data, targets_data, avail_data, gfa_targets = \
        read_assignment_fits_tile(infile)

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

        if len(tgids) != len(set(tgids)):
            raise RuntimeError('TARGETID values in file %s are not unique.' % tf)

        inrows  = np.where(np.isin(tgids, tile_tgids, assume_unique=True))[0]
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

        # FIXME -- work around a bug that is fixed in #375.  Remove once that is merged.
        N = min(len(inrows), len(outrows))
        inrows = inrows[:N]
        outrows = outrows[:N]

        if len(outrows) > 0:
            for c, nm in zip(tfcolsin, tfcolsout):
                tile_targets[nm][outrows] = tgview[c][inrows]
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
        inrows  = np.where(np.isin(skyids, tile_tgids, assume_unique=True))[0]
        outrows = np.where(np.isin(tile_tgids, skyids, assume_unique=True))[0]
        # Demand a unique mapping
        assert(len(inrows) == len(outrows))
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
            for c, nm in zip(tfcolsin, tfcolsout):
                tile_targets[nm][outrows] = skyview[c][inrows]
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

    # Special handling for STUCK positioners that land on good SKY positions.
    # In the FA files these get FA_TYPE & 4 (SKY) and a negative TARGETID.
    stucksky_rows = np.where((fiber_data['TARGETID'] < 0) *
                             ((fiber_data['FA_TYPE'] & TARGET_TYPE_SKY) != 0))[0]
    if len(stucksky_rows):
        outdata['OBJTYPE'][stucksky_rows] = 'SKY'

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
                # AR relaxing the condition on the presence of any
                # AR *_TARGET, excluding FA_TARGET
                # AR hence accepting CMX_.., SV1_..., SV2_...
                #if "DESI_TARGET" in out_dtype.names: # AR commented out
                if "_TARGET" in [name[-7:] for name in out_dtype.names if name!="FA_TARGET"]:
                    # This is a main survey file.
                    objtype[:] = "TGT"
                    is_sky = (tile_targets["DESI_TARGET"][target_rows]
                              & desi_mask.SKY) != 0
                    is_sky |= (tile_targets["DESI_TARGET"][target_rows]
                              & desi_mask.SUPP_SKY) != 0
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

    # AR 20210928 : populate EBV=0 rows in FIBERASSIGN with correct values
    if "EBV" in outdata.dtype.names:
        sel = outdata["EBV"] == 0
        if sel.sum() > 0:
            ra_key, dec_key = merged_fiberassign_swap["RA"], merged_fiberassign_swap["DEC"]
            ras, decs = outdata[ra_key][sel], outdata[dec_key][sel]
            cs = SkyCoord(ra=ras * units.deg, dec=decs * units.deg, frame="icrs")
            ebvs = SFDMap(scaling=1).ebv(cs)
            outdata["EBV"][sel] = ebvs
            log.info("Populating {} rows of FIBERASSIGN with EBV values (min={:.2f}, max={:.2f})".format(
                    sel.sum(), ebvs.min(), ebvs.max()
                )
            )


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

    # Write a header-only primary HDU
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
    # AR hard-cutting on merged_targets_columns + *_TARGET columns
    tile_targets = tile_targets[list(merged_targets_columns.keys())]

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
    cols = merged_fiberassign_req_columns.copy()
    cols.update(merged_fiberassign_req_columns_at_end)
    dcols = [(x, y) for x, y in cols.items()]
    dcolnames = [x for x in cols.keys()]

    tgdata = dict()
    tgdtype = dict()
    tgshape = dict()
    tghead = dict()

    skydata = dict()
    skydtype = dict()
    skyshape = dict()
    skyhead = dict()

    survey = None

    # AR adding any possible *_TARGET column
    # AR future-proofing for SV2, SV3, or other
    # https://github.com/desihub/fiberassign/issues/296
    for tf in targetfiles + skyfiles:
        fd = fitsio.FITS(tf, "r")
        minimal_target_columns.update(
            OrderedDict([
                (key,fd[1].get_rec_dtype()[0][key].str)
                for key in fd[1].get_colnames()
                if key[-7:]=="_TARGET" and key!="FA_TARGET"]))
        fd.close()

    # minimal_target_columns to read
    minimal_dcolnames  = [x for x in minimal_target_columns.keys()]
    minimal_dcols = [(x, y) for x, y in minimal_target_columns.items()]

    for tf in targetfiles:
        tm = Timer()
        tm.start()
        fd = fitsio.FITS(tf)
        tghead[tf] = fd[1].read_header()
        # Allocate a shared memory buffer for the target data
        tglen = fd[1].get_nrows()
        tgshape[tf] = (tglen,)
        #tgdtype[tf], tempoff, tempisvararray = fd[1].get_rec_dtype()

        #select what subset of the 'minimal_dcolnames' are present in the data.
        file_tgdtype, tempoff, tempisvararray = fd[1].get_rec_dtype()
        file_dcolnames  = [x for x in file_tgdtype.names]
        dcols_to_read = []
        for i in range(len(minimal_dcolnames)):
            if minimal_dcolnames[i] in file_dcolnames:
                dcols_to_read.append(minimal_dcols[i])
        some_dt = np.dtype(dcols_to_read)
        some_columns = list(some_dt.fields.keys())

        #print(file_tgdtype)
        tgdtype[tf] = some_dt
        tgbytes = tglen * tgdtype[tf].itemsize
        tgdata[tf] = RawArray("B", tgbytes)
        tgview = np.frombuffer(tgdata[tf],
                               dtype=tgdtype[tf]).reshape(tgshape[tf])
        # Read data directly into shared buffer
        tgview[:] = fd[1].read(columns=some_columns)[some_columns]
        #if survey is None: # AR commented out, not used apparently
        #    (survey, col, sciencemask, stdmask, skymask, suppskymask, # AR commented out, not used apparently
        #     safemask, excludemask) = default_target_masks(tgview) # AR commented out, not used apparently

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

        # AR adding here the minimal set of columns to be read
        # AR just copying what is done for the targets,
        # AR with replacing "tg" by "sky"
        #select what subset of the 'minimal_dcolnames' are present in the data.
        file_skydtype, tempoff, tempisvararray = fd[1].get_rec_dtype()
        file_dcolnames  = [x for x in file_skydtype.names]
        dcols_to_read = []
        for i in range(len(minimal_dcolnames)):
            if minimal_dcolnames[i] in file_dcolnames:
                dcols_to_read.append(minimal_dcols[i])
        some_dt = np.dtype(dcols_to_read)
        some_columns = list(some_dt.fields.keys())

        #print(file_skydtype)
        skydtype[tf] = some_dt
        skybytes = skylen * skydtype[tf].itemsize
        skydata[tf] = RawArray("B", skybytes)

        #skydtype[tf], tempoff, tempisvararray = fd[1].get_rec_dtype() # AR commented out
        #skybytes = skylen * skydtype[tf].itemsize # AR commented out
        #skydata[tf] = RawArray("B", skybytes) # AR commented out
        skyview = np.frombuffer(skydata[tf],
                                dtype=skydtype[tf]).reshape(skyshape[tf])
        # Read data directly into shared buffer
        #skyview[:] = fd[1].read() # AR commented out
        skyview[:] = fd[1].read(columns=some_columns)[some_columns]


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


    end_keys = [k for k,v in merged_fiberassign_req_columns_at_end.items()]
    dcols = ([c for c in dcols if c[0] not in end_keys] +
             [c for c in dcols if c[0] in end_keys])
    out_dtype = np.dtype(dcols)

    # AR adding any *_TARGET columns to the TARGETS columns
    merged_targets_columns.update(
            OrderedDict([
                (name,out_dtype[name]) for name in out_dtype.names
                if name[-7:]=="_TARGET"]))

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

def get_parked_thetaphi(theta_offset, theta_min, theta_max,
                        phi_offset,   phi_min,   phi_max):
    # Arbitrarily,
    theta = theta_offset + theta_min

    # Put phi arm at 150 degrees (https://github.com/desihub/fiberassign/issues/194)
    phi_target = 150.
    phi_target = np.clip(np.deg2rad(phi_target), phi_min, phi_max)
    phi = phi_offset + phi_target

    return theta, phi

def run(
    asgn,
    std_per_petal=10,
    sky_per_petal=40,
    sky_per_slitblock=0,
    start_tile=-1,
    stop_tile=-1,
    redistribute=True,
    use_zero_obsremain=True
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
        sky_per_slitblock (int):  The number of sky to assign per slitblock
        start_tile (int):  If specified, the first tile ID to assign.
        stop_tile (int):  If specified, the last tile ID to assign.
        redistribute (bool):  If True, attempt to shift science targets to unassigned
            fibers on later tiles in order to balance the number per petal.

    Returns:
        None

    """
    gt = GlobalTimers.get()

    log = Logger.get()

    def print_counts(when=None):
        counts = asgn.get_counts(start_tile, stop_tile)
        tiles = list(counts.keys())
        tiles.sort()
        for tile in tiles:
            msg = 'Tile %i: ' % tile
            if when is not None:
                msg += when
            tilecounts = counts[tile]
            keys = [('SCIENCE',True), ('SCIENCE not STANDARD',False), ('STANDARD',True),
                    ('SKY',True), ('SUPPSKY',False), ('SAFE',False)]
            ss = []
            for k,always in keys:
                n = tilecounts.get(k, None)
                if n is None:
                    log.warning('Key', k, 'missing from Assignment.get_counts return value')
                else:
                    if n>0 or always:
                        ss.append('%s: %i' % (k,n))
            log.info(msg + ', '.join(ss))

    print_counts('Start: ')

    # First-pass assignment of science targets
    gt.start("Assign unused fibers to science targets")
    asgn.assign_unused(TARGET_TYPE_SCIENCE, -1, -1, "POS", start_tile, stop_tile)
    gt.stop("Assign unused fibers to science targets")
    print_counts('After assigning unused fibers to science targets: ')

    # Redistribute science targets across available petals
    if redistribute:
        gt.start("Redistribute science targets")
        asgn.redistribute_science(start_tile, stop_tile)
        gt.stop("Redistribute science targets")
        print_counts('After redistributing science targets: ')

    # Assign standards, up to some limit
    gt.start("Assign unused fibers to standards")
    asgn.assign_unused(
        TARGET_TYPE_STANDARD, std_per_petal, -1, "POS", start_tile, stop_tile
    )
    gt.stop("Assign unused fibers to standards")
    print_counts('After assigning standards: ')

    def do_assign_unused_sky(ttype, supp=False):
        tag = 'supp' if supp else ''
        if sky_per_petal > 0 and sky_per_slitblock > 0:
            # Assign using the slitblock requirement first, because it is
            # more specific
            asgn.assign_unused(
                ttype, -1, sky_per_slitblock, "POS",
                start_tile, stop_tile
            )
            print_counts('After assigning %ssky per-slitblock: ' % tag)

            # Then assign using the petal requirement, because it may(should) require
            # more fibers overall.
            asgn.assign_unused(
                ttype, sky_per_petal, -1, "POS",
                start_tile, stop_tile
            )
            print_counts('After assigning %ssky per-petal: ' % tag)
        else:
            asgn.assign_unused(
                ttype, sky_per_petal, sky_per_slitblock, "POS",
                start_tile, stop_tile
            )
            print_counts('After assigning %ssky: ' % tag)

    # Assign sky to unused fibers, up to some limit
    gt.start("Assign unused fibers to sky")
    do_assign_unused_sky(TARGET_TYPE_SKY)
    gt.stop("Assign unused fibers to sky")

    # Assign suppsky to unused fibers, up to some limit
    gt.start("Assign unused fibers to supp_sky")
    do_assign_unused_sky(TARGET_TYPE_SUPPSKY, supp=True)
    gt.stop("Assign unused fibers to supp_sky")

    # Force assignment if needed
    gt.start("Force assignment of sufficient standards")
    asgn.assign_force(
        TARGET_TYPE_STANDARD, std_per_petal, -1, start_tile, stop_tile
    )
    gt.stop("Force assignment of sufficient standards")
    print_counts('After force-assigning standards: ')

    def do_assign_forced_sky(ttype, supp=False):
        tag = 'supp' if supp else ''
        # This function really feels redundant with do_assign_unused_sky, but
        # when I tried to make a single function to do both calls, I had to call
        # f(*(preargs + pos_arg + postargs)) and it looked too mysterious.
        if sky_per_petal > 0 and sky_per_slitblock > 0:
            # Slitblock first
            asgn.assign_force(
                ttype, -1, sky_per_slitblock, start_tile, stop_tile)
            print_counts('After force-assigning %ssky per-slitblock: ' % tag)
            # Then petal
            asgn.assign_force(
                ttype, sky_per_petal, -1, start_tile, stop_tile)
            print_counts('After force-assigning %ssky per-petal: ' % tag)
        else:
            asgn.assign_force(
                ttype, sky_per_petal, sky_per_slitblock, start_tile, stop_tile)
            print_counts('After force-assigning %ssky: ' % tag)

    gt.start("Force assignment of sufficient sky")
    do_assign_forced_sky(TARGET_TYPE_SKY)
    gt.stop("Force assignment of sufficient sky")

    gt.start("Force assignment of sufficient supp_sky")
    do_assign_forced_sky(TARGET_TYPE_SUPPSKY, supp=True)
    gt.stop("Force assignment of sufficient supp_sky")

    # If there are any unassigned fibers, try to place them somewhere.
    # When assigning science targets to these unused fibers, also consider targets
    # with no remaining observations.  Getting extra observations of science
    # targets is preferred over additional standards and sky.  See desi-survey email
    # list archive message 1865 and preceding discussion thread.
    gt.start("Assign remaining unassigned fibers")
    asgn.assign_unused(
        TARGET_TYPE_SCIENCE,
        -1,
        -1,
        "POS",
        start_tile,
        stop_tile,
        use_zero_obsremain=use_zero_obsremain
    )
    print_counts('After assigning reobservations of science targets: ')

    asgn.assign_unused(TARGET_TYPE_STANDARD, -1, -1, "POS", start_tile, stop_tile)
    asgn.assign_unused(TARGET_TYPE_SKY, -1, -1, "POS", start_tile, stop_tile)
    asgn.assign_unused(TARGET_TYPE_SUPPSKY, -1, -1, "POS", start_tile, stop_tile)

    # Assign safe location to unused fibers (no maximum).  There should
    # always be at least one safe location (i.e. "BAD_SKY") for each fiber.
    # So after this is run every fiber should be assigned to something.
    asgn.assign_unused(TARGET_TYPE_SAFE, -1, -1, "POS", start_tile, stop_tile)
    gt.stop("Assign remaining unassigned fibers")
    print_counts('Final assignments: ')

    # Assign sky monitor fibers
    gt.start("Assign sky monitor fibers")
    asgn.assign_unused(TARGET_TYPE_SKY, -1, -1, "ETC", start_tile, stop_tile)
    asgn.assign_unused(TARGET_TYPE_SUPPSKY, -1, -1, "ETC", start_tile, stop_tile)
    asgn.assign_unused(TARGET_TYPE_SAFE, -1, -1, "ETC", start_tile, stop_tile)
    gt.stop("Assign sky monitor fibers")

    return asgn
