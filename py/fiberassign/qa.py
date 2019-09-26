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

from desitarget.targetmask import desi_mask

from .utils import Logger, default_mp_proc

from .targets import (Targets, load_target_table)

from .assign import (result_tiles, result_path, avail_table_to_dict, gfa_table_to_dict,
                     read_assignment_fits_tile)


def qa_parse_table(header, tgdata):
    """Extract target info from a table.
    """
    tgs = Targets()
    if "FA_SURV" in header:
        load_target_table(tgs, tgdata,
                          survey=str(header["FA_SURV"]).rstrip(),
                          typecol="FA_TYPE")
    else:
        load_target_table(tgs, tgdata)
    survey = tgs.survey()

    tgprops = dict()
    if survey == "main":
        lrgmask = int(desi_mask["LRG"].mask)
        elgmask = int(desi_mask["ELG"].mask)
        qsomask = int(desi_mask["QSO"].mask)
        badmask = int(desi_mask["BAD_SKY"].mask)
        for row in range(len(tgdata)):
            tgid = tgdata["TARGETID"][row]
            tg = tgs.get(tgid)
            dt = tg.bits
            tgprops[tgid] = dict()
            if dt & lrgmask:
                tgprops[tgid]["type"] = "LRG"
            elif dt & elgmask:
                tgprops[tgid]["type"] = "ELG"
            elif dt & qsomask:
                tgprops[tgid]["type"] = "QSO"
            elif dt & badmask:
                tgprops[tgid]["type"] = "BAD"
            else:
                tgprops[tgid]["type"] = "NA"
    else:
        # Could define similar things for other surveys here...
        for row in range(len(tgdata)):
            tgid = tgdata["TARGETID"][row]
            tgprops[tgid] = {"type": "NA"}
    return tgs, tgprops


def qa_tile_with_gfa(hw, tile_id, tgs, tgprops, tile_assign, tile_avail, tile_gfa):
    props = dict()
    
    locs = np.array(hw.device_locations("POS"))
    nassign = 0
    nscience = 0
    nstd = 0
    nsky = 0
    nsafe = 0
    unassigned = list()
    objtypes = dict()
       
    # SE to add GFA info to props
    petals = list(range(10))
    petals_wgfa = list(tile_gfa.keys())
    petals_wogfa = list(set(petals)-set(petals_wgfa))  
    gfas_per_tile = [nsafe]*10
    brightest_gfas = [nsafe]*10
    faintest_gfas = [nsafe]*10
    gfas_upto_18th = [nsafe]*10
    gfas_upto_19th = [nsafe]*10  
    gfas_upto_20th = [nsafe]*10  
   
    for cam in petals_wgfa:
        if (len(tile_gfa[cam])>0):
            gfas_per_tile[cam] = len(tile_gfa[cam])
            brightest_gfas[cam] = np.min(tile_gfa[cam])
            faintest_gfas[cam] = np.max(tile_gfa[cam])
            gfas_upto_18th[cam] = sum(np.asarray(tile_gfa[cam])<18)
            gfas_upto_19th[cam] = sum(np.asarray(tile_gfa[cam])<19)
            gfas_upto_20th[cam] = sum(np.asarray(tile_gfa[cam])<20)

        else:
            gfas_per_tile[cam] = [0]
            brightest_gfas[cam] = ['NaN']
            faintest_gfas[cam] = ['NaN']
            gfas_upto_18th[cam] = [0]
            gfas_upto_19th[cam] = [0]
            gfas_upto_20th[cam] = [0] 
               
    for lid in locs:
        if lid not in tile_assign:
            unassigned.append(int(lid))
            continue
        tgid = tile_assign[lid]
        if tgid < 0:
            unassigned.append(int(lid))
            continue
        nassign += 1
        tg = tgs.get(tgid)
        if tg.is_science():
            nscience += 1
            ot = tgprops[tgid]["type"]
            if ot in objtypes:
                objtypes[ot] += 1
            else:
                objtypes[ot] = 1
        if tg.is_standard():
            nstd += 1
        if tg.is_sky():
            nsky += 1
        if tg.is_safe():
            nsafe += 1
    props["assign_total"] = nassign
    props["assign_science"] = nscience
    props["assign_std"] = nstd
    props["assign_sky"] = nsky
    props["assign_safe"] = nsafe
    # SE: added this key for the number of GFA stars per camera: list of 10 integers per tile
    props["gfa_stars_percam"] = gfas_per_tile
    # SE: added this key for the list of 10 magnitudes of the brightest GFA stars available per camera  
    props["brightest_gfa_star_percam"] = [np.round(b*(nsafe+1)/(nsafe+1),2) for b in list(brightest_gfas)] 
    # SE: Number of gfa stars <18 per camera 
    props["gfa_stars_brighter_than_18th"] = [int(b*(nsafe+1)/(nsafe+1)) for b in list(gfas_upto_18th)]
    props["gfa_stars_brighter_than_19th"] = [int(b*(nsafe+1)/(nsafe+1)) for b in list(gfas_upto_19th)]
    props["gfa_stars_brighter_than_20th"] = [int(b*(nsafe+1)/(nsafe+1)) for b in list(gfas_upto_20th)]

    # SE: added this key for the list of 10 magnitudes of the faintest GFA stars available per camera  
    props["faintest_gfa_star_percam"] = [np.round(b*(nsafe+1)/(nsafe+1),2) for b in list(faintest_gfas)]
        
    for ot, cnt in objtypes.items():
        props["assign_obj_{}".format(ot)] = cnt
    props["unassigned"] = unassigned
    return props



def qa_tile(hw, tile_id, tgs, tgprops, tile_assign, tile_avail):
    props = dict()
    # tile_idx = tiles.order[tile_id]
    # props["tile_ra"] = tiles.ra[tile_idx]
    # props["tile_dec"] = tiles.dec[tile_idx]
    # props["tile_obscond"] = tiles.obscond[tile_idx]
    locs = np.array(hw.device_locations("POS"))
    nassign = 0
    nscience = 0
    nstd = 0
    nsky = 0
    nsafe = 0
    unassigned = list()
    objtypes = dict()
    
    for lid in locs:
        if lid not in tile_assign:
            unassigned.append(int(lid))
            continue
        tgid = tile_assign[lid]
        if tgid < 0:
            unassigned.append(int(lid))
            continue
        nassign += 1
        tg = tgs.get(tgid)
        if tg.is_science():
            nscience += 1
            ot = tgprops[tgid]["type"]
            if ot in objtypes:
                objtypes[ot] += 1
            else:
                objtypes[ot] = 1
        if tg.is_standard():
            nstd += 1
        if tg.is_sky():
            nsky += 1
        if tg.is_safe():
            nsafe += 1
    props["assign_total"] = nassign
    props["assign_science"] = nscience
    props["assign_std"] = nstd
    props["assign_sky"] = nsky
    props["assign_safe"] = nsafe
    for ot, cnt in objtypes.items():
        props["assign_obj_{}".format(ot)] = cnt
    props["unassigned"] = unassigned
    return props



def qa_tile_file(hw, params):
    (tile_id, tile_file) = params
    log = Logger.get()

    log.info("Processing tile {}".format(tile_id))

    header, fiber_data, targets_data, avail_data, gfa_data = \
        read_assignment_fits_tile((tile_id, tile_file))
       
    # Target properties
    tgs, tgprops = qa_parse_table(header, targets_data)

    # Only do QA on positioners.
    pos_rows = np.where(fiber_data["DEVICE_TYPE"].astype(str) == "POS")[0]

    # Target assignment
    tassign = {x["LOCATION"]: x["TARGETID"] for x in fiber_data[pos_rows]
               if (x["LOCATION"] >= 0)}

    tavail = avail_table_to_dict(avail_data)
    if gfa_data is not None:  
        tgfa = gfa_table_to_dict(gfa_data)

        qa_data = qa_tile_with_gfa(hw, tile_id, tgs, tgprops, tassign, tavail, tgfa)
    else:
        qa_data = qa_tile(hw, tile_id, tgs, tgprops, tassign, tavail)

    return qa_data

def qa_tiles(hw, tiles, result_dir=".", result_prefix="fiberassign_",
             result_split_dir=False, qa_out=None):
    """Run QA on a set of tiles.

    This will run QA on a set of output files.

    Args:
        hw (Hardware):  the hardware description.
        tiles (list):  List of tile IDs to process.
        result_dir (str):  Top-level directory of fiberassign results.
        result_prefix (str):  Prefix of each per-tile file name.
        result_split_dir (bool):  Results are in split tile directories.
        qa_out (str):  Override the output path.

    Returns:
        None.

    """
    log = Logger.get()

    foundtiles = result_tiles(dir=result_dir, prefix=result_prefix)
    log.info("Found {} fiberassign tile files".format(len(foundtiles)))

    if qa_out is None:
        qa_out = os.path.join(result_dir, "qa.json")

    qa_tile = partial(qa_tile_file, hw)

    avail_tiles = np.array(tiles.id)
    select_tiles = [x for x in foundtiles if x in avail_tiles]

    tile_map_list = [(x, result_path(x, dir=result_dir, prefix=result_prefix,
                                     split=result_split_dir))
                     for x in select_tiles]

    log.info("Selecting {} fiberassign tile files".format(len(tile_map_list)))
    with mp.Pool(processes=default_mp_proc) as pool:
        qa_result = pool.map(qa_tile, tile_map_list)

    qa_full = dict()
    for tid, props in zip(foundtiles, qa_result):
        tidx = tiles.order[tid]
        props["tile_ra"] = np.round(float(tiles.ra[tidx]),2)
        props["tile_dec"] = float(tiles.dec[tidx])
        props["tile_obscond"] = int(tiles.obscond[tidx])
        qa_full[tid] = props

    with open(qa_out, "w") as f:
        json.dump(qa_full, f, indent=4, sort_keys=True)

    return
