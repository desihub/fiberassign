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

from ._version import __version__

from .utils import Logger, Timer, default_mp_proc

from .tiles import load_tiles

from .targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                      TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                      TargetTree, TargetsAvailable, FibersAvailable,
                      load_target_file, desi_target_type)

from .hardware import load_hardware

from ._internal import Assignment


# The columns and types that are written to FITS format.
assign_result_columns = OrderedDict([
    ("TARGETID", "i8"),
    ("FIBER", "i4"),
    ("RA", "f8"),
    ("DEC", "f8"),
    ("FBATYPE", "u1"),
    ("PRIORITY", "i4"),
    ("SUBPRIORITY", "f8"),
    ("OBSCONDITIONS", "i4"),
    ("NUMOBS_MORE", "i4")
])


def write_assignment_fits_tile(outroot, asgn, params):
    """Write a single tile assignment to a FITS file.

    Args:
        outroot (str):  full path of the output root file name.
        asgn (Assignment):  the assignment class instance.
        params (tuple):  tuple containing the tile ID, RA, and DEC.

    """
    # tm = Timer()
    # tm.start()
    tile_id, tile_ra, tile_dec = params
    log = Logger.get()

    # Targets available for all tile / fibers
    tgsavail = asgn.targets_avail()

    # Target properties
    tgs = asgn.targets()

    # Data for this tile
    tdata = asgn.tile_fiber_target(tile_id)
    avail = tgsavail.tile_data(tile_id)

    # The recarray dtype for the assignment and available targets
    assign_dtype = np.dtype([(x, y) for x, y in assign_result_columns.items()])
    avail_dtype = np.dtype([("FIBER", "i4"), ("TARGETID", "i8")])
    # tm.stop()
    # tm.report("  data pointers for tile {}".format(tile_id))

    log.debug("Write:  indexing tile {}".format(tile_id))

    # Reverse lookup

    # tm.clear()
    # tm.start()
    tgfiber = {y: x for x, y in tdata.items()}
    # Compute the total list of targets
    navail = np.sum([len(avail[x]) for x in avail.keys()])
    tgids = np.unique(np.concatenate([np.array(avail[x], dtype=np.int64)
                                      for x in avail.keys()]))
    tgfids = [tgfiber[x] if x in tgfiber.keys() else -1 for x in tgids]

    # tm.stop()
    # tm.report("  indexing tile {}".format(tile_id))

    # Output tile file
    tfile = "{}_{:06d}.fits".format(outroot, tile_id)

    if len(tgids) > 0:
        # tm.clear()
        # tm.start()

        if os.path.isfile(tfile):
            raise RuntimeError("output file {} already exists"
                               .format(tfile))
        # This tile has some available targets
        log.info("Writing tile {}".format(tile_id))
        # Create the file
        fd = fitsio.FITS(tfile, "rw")
        # Construct the output recarray for the assignment and write
        log.debug("Write:  copying assignment data for tile {}"
                  .format(tile_id))

        # tm.stop()
        # tm.report("  opening file for tile {}".format(tile_id))
        #
        # tm.clear()
        # tm.start()
        fdata = np.zeros(len(tgids), dtype=assign_dtype)

        fdata[:] = [(x.id, fid, x.ra, x.dec, x.type, x.priority,
                    x.subpriority, x.obscond, x.obs_remain) for x, fid
                    in zip([tgs.get(y) for y in tgids], tgfids)]
        #
        # off = 0
        # for tg, fid in zip(tgids, tgfids):
        #     props = tgs.get(tg)
        #     obsrem = 0
        #     if props.obs_remain > 0:
        #         obsrem = props.obs_remain
        #     fdata[off] = (tg, fid, props.ra, props.dec, props.type,
        #                   props.priority, props.subpriority, props.obscond,
        #                   obsrem)
        #     off += 1
        # tm.stop()
        # tm.report("  copy data tile {}".format(tile_id))

        header = dict()
        header["TILEID"] = tile_id
        header["TILERA"] = tile_ra
        header["TILEDEC"] = tile_dec
        header["FIELDNUM"] = 0
        header["REQRA"] = tile_ra
        header["REQDEC"] = tile_dec
        header["FBAVER"] = __version__
        log.debug("Write:  FITS write tile {}".format(tile_id))

        # tm.clear()
        # tm.start()
        fd.write(fdata, header=header, extname='FIBERASSIGN_RAW')
        del fdata
        # tm.stop()
        # tm.report("  write / del assignment for tile {}".format(tile_id))

        # Construct recarray for available targets and write
        log.debug("Write:  available targets for tile {}".format(tile_id))

        # tm.clear()
        # tm.start()
        fdata = np.zeros(navail, dtype=avail_dtype)
        off = 0
        for fid in sorted(avail.keys()):
            for tg in avail[fid]:
                fdata[off] = (fid, tg)
                off += 1
        fd.write(fdata, header=header, extname='POTENTIAL_ASSIGNMENTS')
        del fdata
        # tm.stop()
        # tm.report("  write / del avail tile {}".format(tile_id))
    return


def write_assignment_fits(tiles, asgn, out_dir=".", out_prefix="fiberassign"):
    """Write out assignment results in FITS format.

    For each tile, all available targets (not only the assigned targets) and
    their properties are written to the first HDU.  The second HDU contains
    one row for every fiber and target available to that fiber.

    Args:
        tiles (Tiles):  The Tiles object containing the properties of each
            tile.
        asgn (Assignment):  The assignment object.
        out_dir (str):  The output directory for writing per-tile files.
        out_prefix (str):  The output file name prefix.

    """
    tm = Timer()
    tm.start()

    # The tile IDs that were assigned
    tileids = asgn.tiles_assigned()
    tileorder = tiles.order
    tilera = tiles.ra
    tiledec = tiles.dec

    outroot = os.path.join(out_dir, out_prefix)

    write_tile = partial(write_assignment_fits_tile, outroot, asgn)

    tile_map_list = [(tid, tilera[tileorder[tid]],
                      tiledec[tileorder[tid]]) for tid in tileids]

    for params in tile_map_list:
        write_tile(params)

    tm.stop()
    tm.report("Write output files")

    return


def write_assignment_ascii(tiles, asgn, out_dir=".", out_prefix="fiberassign"):
    """Write out assignment results in ASCII format.

    For each tile, only the final assignment to each tile is written out.  For
    the full information, including available targets, use the FITS format.

    Args:
        tiles (Tiles):  The Tiles object containing the properties of each
            tile.
        asgn (Assignment):  The assignment object.
        out_dir (str):  The output directory for writing per-tile files.
        out_prefix (str):  The output file name prefix.

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

    outroot = os.path.join(out_dir, out_prefix)

    for t in tileids:
        tdata = asgn.tile_fiber_target(t)
        nfiber = len(tdata)
        tfile = "{}_{}.txt".format(outroot, t)
        if nfiber > 0:
            log.debug("Writing tile {}".format(t))
            with open(tfile, "w") as f:
                f.write("# TILE_RA = {}\n".format(tilera[tileorder[t]]))
                f.write("# TILE_DEC = {}\n".format(tiledec[tileorder[t]]))
                f.write("#\n")
                f.write("# FIBER  TARGETID  RA  DEC  PRIORITY  "
                        "SUBPRIORITY  OBSCONDITIONS  NUMOBS_MORE  FBATYPE\n")
                for fid in sorted(tdata.keys()):
                    tgid = tdata[fid]
                    tg = tgs.get(tgid)
                    f.write("{:d} {:d} {:.6f} {:.6f}\n"
                            .format(fid, tgid, tg.ra, tg.dec, tg.priority,
                                    tg.subpriority, tg.obscond, tg.obs_remain,
                                    tg.type))
    return


def read_assignment_fits_tile_old(outroot, params):
    """Read in results from the old format.
    """
    (tile_id,) = params
    log = Logger.get()
    # Output tile file
    tfile = "{}_{:05d}.fits".format(outroot, tile_id)
    if not os.path.isfile(tfile):
        raise RuntimeError("input file {} does not exist".format(tfile))
    log.debug("Reading tile data {}".format(tfile))
    # Open the file
    fd = fitsio.FITS(tfile, "r")
    header = fd[1].read_header()
    rawdata = fd[1].read()
    # Repack this data into our expected dtype
    # The recarray dtype for the assignment and available targets
    assign_dtype = np.dtype([(x, y) for x, y in assign_result_columns.items()])
    tgdata = np.empty((fd[1].get_nrows(),), dtype=assign_dtype)
    tgdata["TARGETID"] = rawdata["TARGETID"]
    tgdata["FIBER"] = rawdata["FIBER"]
    if "TARGET_RA" in rawdata.dtype.names:
        tgdata["RA"] = rawdata["TARGET_RA"]
        tgdata["DEC"] = rawdata["TARGET_DEC"]
    else:
        tgdata["RA"] = rawdata["RA"]
        tgdata["DEC"] = rawdata["DEC"]
    tgdata["PRIORITY"] = rawdata["PRIORITY"]
    if "SUBPRIORITY" in rawdata.dtype.names:
        tgdata["SUBPRIORITY"] = rawdata["SUBPRIORITY"]
    else:
        tgdata["SUBPRIORITY"] = 0.0
    if "OBSCONDITIONS" in rawdata.dtype.names:
        tgdata["OBSCONDITIONS"] = rawdata["OBSCONDITIONS"]
    else:
        tgdata["OBSCONDITIONS"] = 0
    if "NUMOBS_MORE" in rawdata.dtype.names:
        tgdata["NUMOBS_MORE"] = rawdata["NUMOBS_MORE"]
    else:
        tgdata["NUMOBS_MORE"] = 0
    tgdata["FBATYPE"] = [desi_target_type(x) for x in rawdata["DESI_TARGET"]]
    del rawdata
    adata = fd[2].read()
    avail = dict()
    if "POTENTIALTARGETID" not in adata.dtype.names:
        # safe to read this.
        for row in adata:
            fid = row["FIBER"]
            tgid = row["TARGETID"]
            if fid in avail:
                avail[fid].append(tgid)
            else:
                avail[fid] = list([tgid])
        avail = {f: np.array(av) for f, av in avail.items()}
    return header, tgdata, avail


def read_assignment_fits_tile(outroot, params):
    """Read in results.
    """
    (tile_id,) = params
    log = Logger.get()
    # Output tile file
    tfile = "{}_{:06d}.fits".format(outroot, tile_id)
    if not os.path.isfile(tfile):
        # Try old naming
        tfile = "{}_{:05d}.fits".format(outroot, tile_id)
        if not os.path.isfile(tfile):
            raise RuntimeError("input file {} does not exist".format(tfile))
    log.debug("Reading tile data {}".format(tfile))
    # Open the file
    fd = fitsio.FITS(tfile, "r")
    header = fd[1].read_header()
    tgdata = fd[1].read()
    adata = fd[2].read()
    avail = dict()
    if "FBATYPE" not in tgdata.dtype.names:
        for row in adata:
            fid = row["LOCATION"]
            tgid = row["TARGETID"]
            if fid in avail:
                avail[fid].append(tgid)
            else:
                avail[fid] = list([tgid])
    else:
        for row in adata:
            fid = row["FIBER"]
            tgid = row["TARGETID"]
            if fid in avail:
                avail[fid].append(tgid)
            else:
                avail[fid] = list([tgid])
    avail = {f: np.array(av) for f, av in avail.items()}
    return header, tgdata, avail


merge_results_tile_tgbuffers = None
merge_results_tile_tgdtypes = None
merge_results_tile_tgshapes = None


def merge_results_tile_initialize(bufs, dtypes, shapes):
    global merge_results_tile_tgbuffers
    global merge_results_tile_tgdtypes
    global merge_results_tile_tgshapes
    merge_results_tile_tgbuffers = bufs
    merge_results_tile_tgdtypes = dtypes
    merge_results_tile_tgshapes = shapes
    return


def merge_results_tile(outroot, out_dtype, params):
    """Merge results for one tile.
    """
    (tile_id,) = params
    log = Logger.get()
    # file names
    infile = "{}_{:06d}.fits".format(outroot, tile_id)
    outfile = "{}_all_{:06d}.fits".format(outroot, tile_id)
    log.info("Reading tile data {}".format(infile))
    infd = fitsio.FITS(infile, "r")
    inhead = infd['FIBERASSIGN'].read_header()
    indata = infd['FIBERASSIGN'].read()
    # Mapping of target ID to row
    tgs = indata["TARGETID"]
    # Construct output recarray
    outdata = np.zeros(len(tgs), dtype=out_dtype)
    # Copy original data
    for field in assign_result_columns:
        outdata[field] = indata[field]
    # Loop over input target files and copy data
    targetfiles = list(merge_results_tile_tgbuffers.keys())
    for tf in targetfiles:
        tgview = np.frombuffer(merge_results_tile_tgbuffers[tf],
                               dtype=merge_results_tile_tgdtypes[tf])\
                               .reshape(merge_results_tile_tgshapes[tf])
        # Some columns may not exist in all target files (e.g. PRIORITY),
        # So we select the valid columns for this file and only copy those.
        log.debug("  {} Computing overlap for {}".format(infile, tf))
        inrows = np.where(np.isin(tgview["TARGETID"], tgs))[0]
        outrows = np.where(np.isin(tgs, tgview["TARGETID"]))[0]
        tfcols = [x for x in out_dtype.names
                  if x in tgview.dtype.names]
        log.debug("  {} Copying data from {}".format(infile, tf))
        if len(outrows) > 0:
            for irw, orw in zip(inrows, outrows):
                for c in tfcols:
                    outdata[c][orw] = tgview[c][irw]
    # Write this HDU
    if os.path.isfile(outfile):
        os.remove(outfile)
    fd = fitsio.FITS(outfile, "rw")
    log.info("Writing new data {}".format(outfile))
    fd.write(outdata, header=inhead, extname='FIBERASSIGN')
    # Copy the available targets HDU
    inhead = infd['POTENTIAL_ASSIGNMENTS'].read_header()
    indata = infd['POTENTIAL_ASSIGNMENTS'].read()
    fd.write(indata, header=inhead, extname='POTENTIAL_ASSIGNMENTS')
    del fd
    del infd
    return


def merge_results(targetfiles, result_dir=".", result_prefix="fiberassign",
                  columns=None):
    """Merge target files and assignment output.
    """
    # import multiprocessing as mp
    log = Logger.get()

    # Find all the per-tile files and get the tile IDs
    tiles = list()
    for root, dirs, files in os.walk(result_dir):
        for f in files:
            mat = re.match(r"{}_(\d+).fits".format(result_prefix), f)
            if mat is not None:
                # Matches the prefix
                tiles.append(int(mat.group(1)))
        break
    log.info("Found {} fiberassign tile files".format(len(tiles)))

    # Load the full set of target files into memory.  Also build a mapping of
    # target ID to row index.  We assume that the result columns have the same
    # dtype in any of the target files.  We take the first target file and
    # construct the output recarray dtype from the columns in that file.
    out_dtype = None
    tgdata = dict()
    tgdtype = dict()
    tgshape = dict()
    tghead = dict()
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

        tm.stop()
        tm.report("Read {} into shared memory".format(tf))
        if out_dtype is None:
            # This is the first target file
            if columns is None:
                columns = list(tgview.dtype.names)
            extra = [x for x in columns
                     if x not in assign_result_columns.keys()]
            dcols = [(x, y) for x, y in assign_result_columns.items()]
            for col in extra:
                subd = tgview.dtype[col].subdtype
                if subd is None:
                    dcols.extend([(col, tgview.dtype[col].str)])
                else:
                    dcols.extend([(col, subd[0], subd[1])])
            out_dtype = np.dtype(dcols)

    # For each tile, find the target IDs used.  Construct the output recarray
    # and copy data into place.

    outroot = os.path.join(result_dir, result_prefix)

    merge_tile = partial(merge_results_tile, outroot, out_dtype)

    tile_map_list = [(x,) for x in tiles]

    # for tid in tiles:
    #     merge_tile((tid,))

    with mp.Pool(processes=default_mp_proc,
                 initializer=merge_results_tile_initialize,
                 initargs=(tgdata, tgdtype, tgshape)) as pool:
        results = pool.map(merge_tile, tile_map_list)

    return
