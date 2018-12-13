# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.qa
=======================

Quality Assurance tools.

"""
from __future__ import absolute_import, division, print_function

import os
import re
import numpy as np

import multiprocessing as mp
from functools import partial

import json

from .utils import Logger, default_mp_proc

from .hardware import load_hardware

# from .tiles import load_tiles

from .targets import (Targets, append_target_table)

from .assign import read_assignment_fits_tile, read_assignment_fits_tile_old


def qa_parse_table(tgdata):
    """Extract target info from a table.
    """
    tgs = Targets()
    typecol = None
    if "FBATYPE" not in tgdata.dtype.names:
        typecol = "DESI_TARGET"
    append_target_table(tgs, tgdata, typecol=typecol)
    return tgs


def qa_tile(hw, tile_id, tgs, tile_assign, tile_avail):
    props = dict()
    # tile_idx = tiles.order[tile_id]
    # props["tile_ra"] = tiles.ra[tile_idx]
    # props["tile_dec"] = tiles.dec[tile_idx]
    # props["tile_obscond"] = tiles.obscond[tile_idx]
    fibers = np.array(hw.fiber_id)
    nassign = 0
    nscience = 0
    nstd = 0
    nsky = 0
    for fid in fibers:
        if fid not in tile_assign:
            continue
        tgid = tile_assign[fid]
        if tgid < 0:
            continue
        nassign += 1
        tg = tgs.get(tgid)
        if tg.is_science():
            nscience += 1
        if tg.is_standard():
            nstd += 1
        if tg.is_sky():
            nsky += 1
    nunassign = len(fibers) - nassign
    props["unassigned"] = nunassign
    props["assign_total"] = nassign
    props["assign_science"] = nscience
    props["assign_std"] = nstd
    props["assign_sky"] = nsky
    return props


def qa_tile_file(hw, inroot, old, params):
    (tile_id, ) = params
    log = Logger()

    log.info("Processing tile {}".format(tile_id))

    if old:
        header, tgdata, tavail = read_assignment_fits_tile_old(inroot,
                                                               (tile_id,))
    else:
        header, tgdata, tavail = read_assignment_fits_tile(inroot, (tile_id,))

    # Target properties
    tgs = qa_parse_table(tgdata)

    # Target assignment
    tassign = {x["FIBER"]: x["TARGETID"] for x in tgdata
               if (x["FIBER"] >= 0)}

    qa_data = qa_tile(hw, tile_id, tgs, tassign, tavail)

    return qa_data


def qa_tiles(resultdir=".", result_prefix="fiberassign", qa_out=None,
             tiles=None, old=False):
    log = Logger()
    # Find all the per-tile files and get the tile IDs
    alltiles = list()
    for root, dirs, files in os.walk(resultdir):
        for f in files:
            mat = re.match(r"{}_(\d+).fits".format(result_prefix), f)
            if mat is not None:
                # Matches the prefix
                alltiles.append(int(mat.group(1)))
        break
    log.info("Found {} fiberassign tile files".format(len(alltiles)))

    inroot = os.path.join(resultdir, result_prefix)
    if qa_out is None:
        qa_out = os.path.join(resultdir, "qa.json")

    hw = load_hardware()

    qa_tile = partial(qa_tile_file, hw, inroot, old)

    if tiles is None:
        tiles = alltiles
    tile_map_list = [(x,) for x in tiles]

    log.info("Selecting {} fiberassign tile files".format(len(tile_map_list)))
    with mp.Pool(processes=default_mp_proc) as pool:
        qa_result = pool.map(qa_tile, tile_map_list)

    qa_total = {x: y for x, y in zip(tiles, qa_result)}

    with open(qa_out, "w") as f:
        json.dump(qa_total, f, indent=4, sort_keys=True)

    return
