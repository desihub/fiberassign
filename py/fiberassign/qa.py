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

from collections import OrderedDict

import fitsio

from desitarget.targetmask import desi_mask

from .utils import Logger, default_mp_proc

from .targets import (
    Targets,
    load_target_table,
    default_main_sciencemask,
    default_main_stdmask,
    default_main_skymask,
    default_main_suppskymask,
    default_main_safemask,
    default_sv1_sciencemask,
    default_sv1_stdmask,
    default_sv1_skymask,
    default_sv1_suppskymask,
    default_sv1_safemask,
    default_cmx_sciencemask,
    default_cmx_stdmask,
    default_cmx_skymask,
    default_cmx_suppskymask,
    default_cmx_safemask,
)

from .assign import (
    result_tiles,
    result_path,
    avail_table_to_dict,
    gfa_table_to_dict,
    read_assignment_fits_tile,
)


def qa_parse_table(header, tgdata):
    """Extract target info from a table.
    """
    log = Logger.get()
    tgs = Targets()
    if "FA_SURV" in header:
        load_target_table(
            tgs, tgdata, survey=str(header["FA_SURV"]).rstrip(), typecol="FA_TYPE"
        )
    else:
        load_target_table(tgs, tgdata)
    survey = tgs.survey()

    lrgmask = int(desi_mask["LRG"].mask)
    elgmask = int(desi_mask["ELG"].mask)
    qsomask = int(desi_mask["QSO"].mask)
    stdfaintmask = int(desi_mask["STD_FAINT"].mask)
    stdwdmask = int(desi_mask["STD_WD"].mask)
    stdbrtmask = int(desi_mask["STD_BRIGHT"].mask)

    sciencemask = 0
    stdmask = 0
    skymask = 0
    suppskymask = 0
    safemask = 0

    tgprops = dict()
    if survey == "main":
        sciencemask = int(default_main_sciencemask())
        stdmask = int(default_main_stdmask())
        skymask = int(default_main_skymask())
        suppskymask = int(default_main_suppskymask())
        safemask = int(default_main_safemask())
    elif survey == "sv1":
        sciencemask = int(default_sv1_sciencemask())
        stdmask = int(default_sv1_stdmask())
        skymask = int(default_sv1_skymask())
        suppskymask = int(default_sv1_suppskymask())
        safemask = int(default_sv1_safemask())
    elif survey == "cmx":
        sciencemask = int(default_cmx_sciencemask())
        stdmask = int(default_cmx_stdmask())
        skymask = int(default_cmx_skymask())
        suppskymask = int(default_cmx_suppskymask())
        safemask = int(default_cmx_safemask())

    obscol = None
    if "NUMOBS_MORE" in tgdata.dtype.names:
        obscol = "NUMOBS_MORE"
    elif "NUMOBS_INIT" in tgdata.dtype.names:
        obscol = "NUMOBS_INIT"
    if obscol is None:
        msg = "Running QA on non-merged data- target QA will be pointless."
        log.warning(msg)

    for row in range(len(tgdata)):
        tgid = tgdata["TARGETID"][row]
        obsinit = -1
        if obscol is not None:
            obsinit = tgdata[obscol][row]
        tg = tgs.get(tgid)
        dt = tg.bits
        tgprops[tgid] = dict()
        if dt & stdmask:
            tgprops[tgid]["class"] = "standard"
            if dt & sciencemask:
                tgprops[tgid]["class"] = "science-standard"
            if dt & stdfaintmask:
                tgprops[tgid]["type"] = "STD_FAINT"
            elif dt & stdwdmask:
                tgprops[tgid]["type"] = "STD_WD"
            elif dt & stdbrtmask:
                tgprops[tgid]["type"] = "STD_BRIGHT"
            else:
                tgprops[tgid]["type"] = "NA"
        elif dt & sciencemask:
            tgprops[tgid]["class"] = "science"
            if dt & lrgmask:
                tgprops[tgid]["type"] = "LRG"
            elif dt & elgmask:
                tgprops[tgid]["type"] = "ELG"
            elif dt & qsomask:
                tgprops[tgid]["type"] = "QSO"
            else:
                tgprops[tgid]["type"] = "NA"
        elif dt & skymask:
            tgprops[tgid]["class"] = "sky"
            tgprops[tgid]["type"] = "SKY"
        elif dt & suppskymask:
            tgprops[tgid]["class"] = "suppsky"
            tgprops[tgid]["type"] = "SSKY"
        elif dt & safemask:
            tgprops[tgid]["class"] = "safe"
            tgprops[tgid]["type"] = "SAFE"
        else:
            tgprops[tgid]["class"] = "unknown"
            tgprops[tgid]["type"] = "NA"

        tgprops[tgid]["obsinit"] = obsinit
    return tgs, tgprops


def qa_tile_with_gfa(hw, tile_id, tgs, tgprops, tile_assign, tile_avail, tile_gfa):
    props = dict()

    locs = np.array(hw.device_locations("POS"))
    nassign = 0
    nscience = 0
    nstd = 0
    nsky = 0
    nsuppsky = 0
    nsafe = 0
    unassigned = list()
    objtypes = dict()

    # SE to add GFA info to props
    petals = list(range(10))
    petals_wgfa = list(tile_gfa.keys())
    petals_wogfa = list(set(petals) - set(petals_wgfa))
    gfas_per_tile = [nsafe] * 10
    brightest_gfas = [nsafe] * 10
    faintest_gfas = [nsafe] * 10
    gfas_upto_18th = [nsafe] * 10
    gfas_upto_19th = [nsafe] * 10
    gfas_upto_20th = [nsafe] * 10

    for cam in petals_wgfa:
        if len(tile_gfa[cam]) > 0:
            gfas_per_tile[cam] = len(tile_gfa[cam])
            brightest_gfas[cam] = np.min(tile_gfa[cam])
            faintest_gfas[cam] = np.max(tile_gfa[cam])
            gfas_upto_18th[cam] = sum(np.asarray(tile_gfa[cam]) < 18)
            gfas_upto_19th[cam] = sum(np.asarray(tile_gfa[cam]) < 19)
            gfas_upto_20th[cam] = sum(np.asarray(tile_gfa[cam]) < 20)

        else:
            gfas_per_tile[cam] = [0]
            brightest_gfas[cam] = ["NaN"]
            faintest_gfas[cam] = ["NaN"]
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
        if tg.is_suppsky():
            nsuppsky += 1
        if tg.is_safe():
            nsafe += 1
    props["assign_total"] = nassign
    props["assign_science"] = nscience
    props["assign_std"] = nstd
    props["assign_sky"] = nsky
    props["assign_suppsky"] = nsuppsky
    props["assign_safe"] = nsafe
    # SE: added this key for the number of GFA stars per camera: list of 10 integers per tile
    props["gfa_stars_percam"] = gfas_per_tile
    # SE: added this key for the list of 10 magnitudes of the brightest GFA stars available per camera
    props["brightest_gfa_star_percam"] = [
        np.round(b * (nsafe + 1) / (nsafe + 1), 2) for b in list(brightest_gfas)
    ]
    # SE: Number of gfa stars <18 per camera
    props["gfa_stars_brighter_than_18th"] = [
        int(b * (nsafe + 1) / (nsafe + 1)) for b in list(gfas_upto_18th)
    ]
    props["gfa_stars_brighter_than_19th"] = [
        int(b * (nsafe + 1) / (nsafe + 1)) for b in list(gfas_upto_19th)
    ]
    props["gfa_stars_brighter_than_20th"] = [
        int(b * (nsafe + 1) / (nsafe + 1)) for b in list(gfas_upto_20th)
    ]

    # SE: added this key for the list of 10 magnitudes of the faintest GFA stars available per camera
    props["faintest_gfa_star_percam"] = [
        np.round(b * (nsafe + 1) / (nsafe + 1), 2) for b in list(faintest_gfas)
    ]

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
    nsuppsky = 0
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
        if tg.is_suppsky():
            nsuppsky += 1
        if tg.is_safe():
            nsafe += 1
    props["assign_total"] = nassign
    props["assign_science"] = nscience
    props["assign_std"] = nstd
    props["assign_sky"] = nsky
    props["assign_suppsky"] = nsuppsky
    props["assign_safe"] = nsafe
    for ot, cnt in objtypes.items():
        props["assign_obj_{}".format(ot)] = cnt
    props["unassigned"] = unassigned
    return props


def qa_tile_file(hw, params):
    (tile_id, tile_file) = params
    log = Logger.get()

    log.info("Processing tile {}".format(tile_id))

    header, fiber_data, targets_data, avail_data, gfa_data = read_assignment_fits_tile(
        (tile_file)
    )

    # Target properties
    tgs, tgprops = qa_parse_table(header, targets_data)

    # Only do QA on positioners.
    pos_rows = np.where(fiber_data["DEVICE_TYPE"].astype(str) == "POS")[0]

    # Target assignment
    tassign = {
        x["LOCATION"]: x["TARGETID"]
        for x in fiber_data[pos_rows]
        if (x["LOCATION"] >= 0)
    }

    tavail = avail_table_to_dict(avail_data)
    if gfa_data is not None:
        tgfa = gfa_table_to_dict(gfa_data)
        qa_data = qa_tile_with_gfa(hw, tile_id, tgs, tgprops, tassign, tavail, tgfa)
    else:
        qa_data = qa_tile(hw, tile_id, tgs, tgprops, tassign, tavail)

    return qa_data


def qa_tiles(
    hw, tiles, result_dir=".", result_prefix="fba-", result_split_dir=False, qa_out=None
):
    """Run QA on a set of tiles.

    This will run QA on a set of output files.

    Args:
        hw (Hardware):  the hardware description.
        tiles (Tiles):  a Tiles object.
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

    tile_map_list = [
        (
            x,
            result_path(
                x, dir=result_dir, prefix=result_prefix, split=result_split_dir
            ),
        )
        for x in select_tiles
    ]

    log.info("Selecting {} fiberassign tile files".format(len(tile_map_list)))
    with mp.Pool(processes=default_mp_proc) as pool:
        qa_result = pool.map(qa_tile, tile_map_list)

    qa_full = dict()
    for tid, props in zip(foundtiles, qa_result):
        tidx = tiles.order[tid]
        props["tile_ra"] = np.round(float(tiles.ra[tidx]), 2)
        props["tile_dec"] = float(tiles.dec[tidx])
        props["tile_obscond"] = int(tiles.obscond[tidx])
        qa_full[tid] = props

    with open(qa_out, "w") as f:
        json.dump(qa_full, f, indent=4, sort_keys=True)

    return


def qa_targets(
    hw, tiles, result_dir=".", result_prefix="fiberassign-", result_split_dir=False, qa_out=None
):
    """Run target QA on a set of tiles.

    This will run go through a list of fiberassign tiles and track the
    assignment of targets through the survey.

    Args:
        hw (Hardware):  the hardware description.
        tiles (Tiles):  a Tiles object.
        result_dir (str):  Top-level directory of fiberassign results.
        result_prefix (str):  Prefix of each per-tile file name.
        result_split_dir (bool):  Results are in split tile directories.
        out (str):  Override the output path.

    Returns:
        None.

    """
    log = Logger.get()

    if qa_out is None:
        qa_out = os.path.join(result_dir, "qa_target_tilecount.json")

    # We cannot have missing tiles or else it will skew our picture
    # of how targets are assigned.

    tile_ids = np.array(tiles.id)
    log.info("Using survey with {} fiberassign tiles".format(len(tile_ids)))

    tile_files = [
        (
            x,
            result_path(
                x, dir=result_dir, prefix=result_prefix, split=result_split_dir
            ),
        )
        for x in np.array(tiles.id)
    ]

    # Pre-check that all files exists

    for tid, tf in tile_files:
        if not os.path.isfile(tf):
            msg = "fiberassign output file '{}' does not exist.".format(tf)
            raise RuntimeError(msg)

    # We want to track assignment history for each target and for classes of
    # targets.

    qadata = dict()
    qatargets = dict()
    qaavail = dict()

    # We must process the tiles in sequence, so do not use multiprocessing.

    for tid, tf in tile_files:
        log.info("Processing tile {}".format(tid))
        itid = int(tid)

        header, fiber_data, targets_data, avail_data, gfa_data = \
            read_assignment_fits_tile((tf))

        # Target properties
        tgs, tgprops = qa_parse_table(header, targets_data)

        # Only do QA on positioners.
        pos_rows = np.where(fiber_data["DEVICE_TYPE"].astype(str) == "POS")[0]

        # Assigned Targets
        tassign = np.unique([
            x["TARGETID"] for x in fiber_data[pos_rows]
            if (x["LOCATION"] >= 0)
        ])

        # Unique set of targets available for this tile
        tavail = np.unique(avail_data["TARGETID"])

        for tgid in tavail:
            if tgid < 0:
                continue
            if tgid not in tgprops:
                msg = "Available target not in target props.  "\
                    "Rerun fba_run with the --write_all_targets option."
                raise RuntimeError()
            cls = str(tgprops[tgid]["class"])
            typ = str(tgprops[tgid]["type"])
            obi = int(tgprops[tgid]["obsinit"])

            if cls not in qadata:
                # first occurrence of this class (science, sky, etc)
                qadata[cls] = dict()
                qatargets[cls] = dict()
                qaavail[cls] = dict()
            if typ not in qadata[cls]:
                # first occurrence of this type (elg, qso, lrg, etc)
                qadata[cls][typ] = dict()
                qatargets[cls][typ] = dict()
                qaavail[cls][typ] = dict()
            if obi not in qadata[cls][typ]:
                # first occurrence of this NUMOBS_INIT
                qadata[cls][typ][obi] = dict()
                qatargets[cls][typ][obi] = dict()
                qaavail[cls][typ][obi] = dict()
            if tgid not in qaavail[cls][typ][obi]:
                # first occurrence of this target
                qaavail[cls][typ][obi][tgid] = int(0)
            qaavail[cls][typ][obi][tgid] += int(1)
            if tgid in tassign:
                if tgid not in qatargets[cls][typ][obi]:
                    # first occurrence of this target
                    qatargets[cls][typ][obi][tgid] = int(0)
                qatargets[cls][typ][obi][tgid] += int(1)
                if tid not in qadata[cls][typ][obi]:
                    # first tile for this NUMOBS_INIT
                    qadata[cls][typ][obi][itid] = int(0)
                qadata[cls][typ][obi][itid] += int(1)

    with open(qa_out, "w") as f:
        json.dump(qadata, f, indent=4, sort_keys=True)

    tcols = OrderedDict([
        ("TARGETID", "i8"),
        ("NUMOBS_AVAIL", "i4"),
        ("NUMOBS_DONE", "i4"),
    ])
    tdtype = np.dtype([(x, y) for x, y in tcols.items()])
    for tcls, clsprops in qaavail.items():
        aclsprops = None
        if tcls in qatargets:
            aclsprops = qatargets[tcls]
        for ttyp, typprops in clsprops.items():
            atypprops = None
            if aclsprops is not None:
                if ttyp in aclsprops:
                    atypprops = aclsprops[ttyp]
            for tobi, obitgs in typprops.items():
                aobitgs = None
                if atypprops is not None:
                    if tobi in atypprops:
                        aobitgs = atypprops[tobi]
                fobi = tobi
                if fobi < 0:
                    fobi = "NA"
                tout = os.path.join(
                    result_dir, "qa_target_count_{}-{}_init-{}.fits".format(
                        tcls, ttyp, fobi
                    )
                )
                if os.path.isfile(tout):
                    os.remove(tout)
                tdata = np.zeros(len(obitgs), dtype=tdtype)
                for row, (tgid, hits) in enumerate(obitgs.items()):
                    tdata[row]["TARGETID"] = tgid
                    tdata[row]["NUMOBS_AVAIL"] = hits
                    tdata[row]["NUMOBS_DONE"] = 0
                    if aobitgs is not None and tgid in aobitgs:
                        tdata[row]["NUMOBS_DONE"] = aobitgs[tgid]

                fd = fitsio.FITS(tout, "rw")
                fd.write(None, header=None, extname="PRIMARY")
                fd.write(tdata, header=None, extname="COUNTS")
                fd.close()

    return
