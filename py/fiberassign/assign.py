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
from functools import partial

from collections import OrderedDict

import fitsio

from ._version import __version__

from .utils import Logger, Timer, default_mp_proc

from .tiles import load_tiles

from .targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                      TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                      TargetTree, TargetsAvailable, FibersAvailable,
                      load_target_file)

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


def write_assignment_fits_tile(outroot, tgs, tile_id, tile_ra, tile_dec,
                               tdata, avail):
    log = Logger()
    # The recarray dtype for the assignment and available targets
    assign_dtype = np.dtype([(x, y) for x, y in assign_result_columns.items()])
    avail_dtype = np.dtype([("FIBER", "i4"), ("TARGETID", "i8")])
    # Reverse lookup
    tgfiber = {y: x for x, y in tdata.items()}
    # Compute the total list of targets
    alltgid = set()
    navail = 0
    for f, av in avail.items():
        alltgid.update(av)
        navail += len(av)
    tgids = list(sorted(alltgid))

    # Output tile file
    tfile = "{}_{:06d}.fits".format(outroot, tile_id)
    if len(tgids) > 0:
        if os.path.isfile(tfile):
            raise RuntimeError("output file {} already exists"
                               .format(tfile))
        # This tile has some available targets
        log.info("Writing tile {}".format(tile_id))
        # Create the file
        fd = fitsio.FITS(tfile, "rw")
        # Construct the output recarray for the assignment and write
        fdata = np.zeros(len(tgids), dtype=assign_dtype)
        off = 0
        for tg in tgids:
            fid = -1
            if tg in tgfiber:
                fid = tgfiber[tg]
            props = tgs.get(tg)
            fdata[off] = (tg, fid, props.ra, props.dec, props.type,
                          props.priority, props.subpriority, props.obscond,
                          props.obs_remain)
            off += 1
        header = dict()
        header["TILE_RA"] = tile_ra
        header["TILE_DEC"] = tile_dec
        header["FBAVER"] = __version__
        fd.write(fdata, header=header)
        del fdata
        # Construct recarray for available targets and write
        fdata = np.zeros(navail, dtype=avail_dtype)
        off = 0
        for fid in sorted(avail.keys()):
            for tg in avail[fid]:
                fdata[off] = (fid, tg)
                off += 1
        fd.write(fdata, header=header)
        del fdata


def write_assignment_fits(tiles, asgn, outdir=".", out_prefix="fiberassign"):
    """Write out results in FITS format.

    """
    tm = Timer()
    tm.start()

    # Targets available for all tile / fibers
    tgsavail = asgn.targets_avail()

    # Target properties
    tgs = asgn.targets()

    # The tile IDs that were assigned
    tileids = asgn.tiles_assigned()
    tileorder = tiles.order
    tilera = tiles.ra
    tiledec = tiles.dec

    outroot = os.path.join(outdir, out_prefix)

    for tid in tileids:
        write_assignment_fits_tile(outroot, tgs, tid, tilera[tileorder[tid]],
                                   tiledec[tileorder[tid]],
                                   asgn.tile_fiber_target(tid),
                                   tgsavail.tile_data(tid))

    tm.stop()
    tm.report("Write output files")

    return


def write_assignment_ascii(tiles, asgn, outdir=".", out_prefix="fiberassign",
                           single=False):
    """Write out results.

    """
    log = Logger()
    # Go through the assignment, one tile at a time.  For each tile, get the
    # best assignment and potential targets.

    tileids = tiles.id

    outroot = os.path.join(outdir, out_prefix)

    if single:
        outfile = "{}.txt".format(outroot)
        with open(outfile, "w") as f:
            for t in tileids:
                tdata = asgn.tile_fiber_target(t)
                nfiber = len(tdata)
                if nfiber > 0:
                    log.debug("Writing tile {}".format(t))
                    for fid in sorted(tdata.keys()):
                        f.write("{} {} {}\n".format(t, fid, tdata[fid]))
    else:
        for t in tileids:
            tdata = asgn.tile_fiber_target(t)
            nfiber = len(tdata)
            tfile = "{}_{}.txt".format(outroot, t)
            if nfiber > 0:
                log.debug("Writing tile {}".format(t))
                with open(tfile, "w") as f:
                    for fid in sorted(tdata.keys()):
                        f.write("{} {}\n".format(fid, tdata[fid]))
    return


def read_assignment_fits_tile(outroot, params):
    """Read in results.
    """
    (tile_id,) = params
    log = Logger()
    # Output tile file
    tfile = "{}_{:06d}.fits".format(outroot, tile_id)
    log.debug("Reading tile data {}".format(tfile))
    if not os.path.isfile(tfile):
        raise RuntimeError("input file {} does not exist".format(tfile))
    # Open the file
    fd = fitsio.FITS(tfile, "r")
    header = fd[1].read_header()
    tgdata = fd[1].read()
    tgindx = {y: x for x, y in enumerate(tgdata["TARGETID"])}
    adata = fd[2].read()
    avail = dict()
    for row in adata:
        fid = row["FIBER"]
        tgid = row["TARGETID"]
        if fid in avail:
            avail[fid].append(tgid)
        else:
            avail[fid] = list([tgid])
    return header, tgindx, tgdata, avail


def merge_results_tile(outroot, out_dtype, tgdata, tgrow, params):
    """Merge results for one tile.
    """
    (tile_id,) = params
    log = Logger()
    # file names
    infile = "{}_{:06d}.fits".format(outroot, tile_id)
    outfile = "{}_all_{:06d}.fits".format(outroot, tile_id)
    log.info("Reading tile data {}".format(infile))
    infd = fitsio.FITS(infile, "r")
    inhead = infd[1].read_header()
    indata = infd[1].read()
    # Mapping of target ID to row
    tgs = indata["TARGETID"]
    trow = {x: y for y, x in enumerate(tgs)}
    # Construct output recarray
    outdata = np.zeros(len(tgs), dtype=out_dtype)
    # Copy original data
    for field in assign_result_columns:
        outdata[field] = indata[field]
    # Loop over input target files and copy data
    targetfiles = list(tgdata.keys())
    for tf in targetfiles:
        # Some columns may not exist in all target files (e.g. PRIORITY),
        # So we select the valid columns for this file and only copy those.
        tfcols = [x for x in out_dtype.names
                  if x in tgdata[tf].dtype.names]
        overlap = np.isin(tgs, tgdata[tf]["TARGETID"])
        inrows = np.array([tgrow[tf][x] for x in tgs[overlap]])
        outrows = np.array([trow[x] for x in tgs[overlap]])
        if len(outrows) > 0:
            for irw, orw in zip(inrows, outrows):
                for c in tfcols:
                    outdata[c][orw] = tgdata[tf][c][irw]
    # Write this HDU
    if os.path.isfile(outfile):
        os.remove(outfile)
    fd = fitsio.FITS(outfile, "rw")
    log.info("Writing new data {}".format(outfile))
    fd.write(outdata, header=inhead)
    # Copy the available targets HDU
    inhead = infd[2].read_header()
    indata = infd[2].read()
    fd.write(indata, header=inhead)
    del fd
    del infd
    return


def merge_results(targetfiles, resultdir=".", result_prefix="fiberassign",
                  columns=None):
    """Merge target files and assignment output.
    """
    # import multiprocessing as mp
    log = Logger()

    # Find all the per-tile files and get the tile IDs
    tiles = list()
    for root, dirs, files in os.walk(resultdir):
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
    tghead = dict()
    tgrow = dict()
    for tf in targetfiles:
        tm = Timer()
        tm.start()
        tgdata[tf], tghead[tf] = fitsio.read(tf, ext=1, header=True)
        tm.stop()
        tm.report("Read {} into memory".format(tf))
        if out_dtype is None:
            # This is the first target file
            if columns is None:
                columns = list(tgdata[tf].dtype.names)
            extra = [x for x in columns
                     if x not in assign_result_columns.keys()]
            dcols = [(x, y) for x, y in assign_result_columns.items()]
            dcols.extend([(x, tgdata[tf].dtype[x].str) for x in extra])
            out_dtype = np.dtype(dcols)
        tm.clear()
        tm.start()
        tgrow[tf] = {x: y for y, x in enumerate(tgdata[tf]["TARGETID"])}
        tm.stop()
        tm.report("Build row index for {}".format(tf))

    # For each tile, find the target IDs used.  Construct the output recarray
    # and copy data into place.

    outroot = os.path.join(resultdir, result_prefix)

    merge_tile = partial(merge_results_tile, outroot, out_dtype, tgdata, tgrow)

    tile_map_list = [(x,) for x in tiles]

    # for tid in tiles:
    #     merge_tile((tid,))

    with mp.Pool(processes=default_mp_proc) as pool:
        results = pool.map(merge_tile, tile_map_list)

    return
