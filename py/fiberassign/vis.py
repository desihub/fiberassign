# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.vis
=======================

Visualization tools.

"""
from __future__ import absolute_import, division, print_function

import os

import numpy as np

import matplotlib
# try:
#     matplotlib.use("module://mplcairo.base")
# except ModuleNotFoundError:
#     try:
#         matplotlib.use("cairo")
#     except ModuleNotFoundError:
#         matplotlib.use("svg")
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import multiprocessing as mp
from functools import partial

import fitsio

from .utils import Logger, default_mp_proc

from .hardware import load_hardware

from .tiles import load_tiles

from .targets import (Targets, append_target_table, TARGET_TYPE_SCIENCE,
                      TARGET_TYPE_SKY, TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE)

from .assign import (read_assignment_fits_tile, result_tiles, result_path,
                     avail_table_to_dict)


def plot_target_type_color(tgtype):
    color = "gray"
    tp = int(tgtype)
    if (tp & TARGET_TYPE_SAFE) != 0:
        color = "black"
    elif (tp & TARGET_TYPE_SKY) != 0:
        color = "blue"
    elif (tp & TARGET_TYPE_STANDARD) != 0:
        color = "gold"
        if (tp & TARGET_TYPE_SCIENCE) != 0:
            color = "orange"
    elif (tp & TARGET_TYPE_SCIENCE) != 0:
        color = "red"
    return color


def plot_positioner(ax, patrol_rad, fiber, center, cb, fh, color="k",
                    linewidth=0.2):
    """Plot one fiber positioner.
    """
    patrol = plt.Circle((center[0], center[1]), radius=patrol_rad, fc=color,
                        ec="none", alpha=0.07)
    ax.add_artist(patrol)
    # Plot the arm from the center to the body
    cbcent = cb.circles[0].center
    armwidth = 0.25
    armlen = np.sqrt((cbcent[1] - center[1])**2 + (cbcent[0] - center[0])**2)
    armang = np.arctan2(cbcent[1] - center[1], cbcent[0] - center[0])
    sinarm = np.sin(armang)
    cosarm = np.cos(armang)
    arm_xoff = center[0] + (0.5*armwidth) * sinarm
    arm_yoff = center[1] - (0.5*armwidth) * cosarm
    armang_deg = armang * 180.0 / np.pi
    arm = plt.Rectangle((arm_xoff, arm_yoff), armlen, armwidth,
                        angle=armang_deg, color=color, linewidth=2*linewidth,
                        fill=False)
    ax.add_artist(arm)
    for piece in [cb, fh]:
        for circle in piece.circles:
            xcent, ycent = circle.center
            rad = circle.radius
            circ = plt.Circle((xcent, ycent), radius=rad, fc="none", ec=color,
                              linewidth=linewidth)
            ax.add_artist(circ)
        for segs in piece.segments:
            xpts = np.array([p[0] for p in segs.points])
            ypts = np.array([p[1] for p in segs.points])
            ax.plot(xpts, ypts, linewidth=linewidth, color=color)
    fontpix = armwidth * 2
    fontpt = int(0.25 * fontpix)
    xtxt = center[0] - 2 * armwidth * cosarm
    ytxt = center[1] - 2 * armwidth * sinarm
    ax.text(xtxt, ytxt, "{}".format(fiber),
            color='k', fontsize=fontpt,
            horizontalalignment='center',
            verticalalignment='center',
            bbox=None)
    #        bbox=dict(fc='w', ec='none', pad=1, alpha=1.0))
    return


def plot_positioner_simple(ax, patrol_rad, fiber, center, cb, fh, color="k",
                           linewidth=0.2):
    """Plot one fiber positioner.

    This uses a simpler representation of the positioner geometry, in order to
    speed up the plotting.

    """
    patrol = plt.Circle((center[0], center[1]), radius=patrol_rad, fc=color,
                        ec="none", alpha=0.07)
    ax.add_artist(patrol)
    # Plot the arm from the center to the body
    cbcent = cb.circles[0].center

    ax.plot([center[0], cbcent[0]], [center[1], cbcent[1]], color=color,
            linewidth=4*linewidth)

    tgloc = fh.circles[0].center

    ax.plot([cbcent[0], tgloc[0]], [cbcent[1], tgloc[1]], color=color,
            linewidth=linewidth)

    fontpt = 0.125
    xtxt = center[0]
    ytxt = center[1] + 0.5
    ax.text(xtxt, ytxt, "{}".format(fiber),
            color='k', fontsize=fontpt,
            horizontalalignment='center',
            verticalalignment='center',
            bbox=None)
    return


def plot_parse_table(tgdata):
    """Create a Targets object from a recarray.
    """
    tgs = Targets()
    typecol = None
    if "FBATYPE" not in tgdata.dtype.names:
        typecol = "DESI_TARGET"
    append_target_table(tgs, tgdata, typecol=typecol)
    return tgs


def plot_tile_targets_props(hw, tile_ra, tile_dec, tgs):
    avail_tgid = tgs.ids()
    # print("  DBG:  avail_tgid len = ", len(avail_tgid), flush=True)
    ra = np.empty(len(avail_tgid), dtype=np.float64)
    dec = np.empty(len(avail_tgid), dtype=np.float64)
    color = list()
    for idx, tgid in enumerate(avail_tgid):
        tg = tgs.get(tgid)
        # print("  DBG:  ",idx," tgid ",tgid," = ",tg, flush=True)
        ra[idx] = tg.ra
        dec[idx] = tg.dec
        color.append(plot_target_type_color(tg.type))
    # We disable threading here, since it does not interact well with
    # multiprocessing.
    tgxy = hw.radec2xy_multi(tile_ra, tile_dec, ra, dec, 1)
    props = {tgid: {"xy": xy, "color": cl} for tgid, xy, cl
             in zip(avail_tgid, tgxy, color)}
    return props


def plot_available(ax, targetprops, selected, linewidth=0.1):
    mwidth = 5.0 * linewidth
    xdata = np.empty(len(selected), dtype=np.float64)
    ydata = np.empty(len(selected), dtype=np.float64)
    color = list()
    for idx, tgid in enumerate(selected):
        xdata[idx] = targetprops[tgid]["xy"][0]
        ydata[idx] = targetprops[tgid]["xy"][1]
        color.append(targetprops[tgid]["color"])
    ax.scatter(xdata, ydata, color=color, marker=".",
               linewidth=linewidth, s=mwidth)
    return


def plot_assignment(ax, hw, targetprops, tile_assigned, linewidth=0.1,
                    real_shapes=False):
    center_mm = hw.fiber_pos_xy_mm
    patrol_rad = hw.patrol_mm
    assigned = np.array(sorted(tile_assigned.keys()), dtype=np.int32)
    for fid in assigned:
        center = center_mm[fid]
        tgid = tile_assigned[fid]
        if tgid >= 0:
            # This fiber is assigned.  Plot the positioner located at the
            # assigned target.
            cb, fh = hw.fiber_position(fid, targetprops[tgid]["xy"])
            if real_shapes:
                plot_positioner(ax, patrol_rad, fid, center, cb, fh,
                                color=targetprops[tgid]["color"],
                                linewidth=linewidth)
            else:
                plot_positioner_simple(ax, patrol_rad, fid, center, cb, fh,
                                       color=targetprops[tgid]["color"],
                                       linewidth=linewidth)
        else:
            # This fiber is unassigned.  Plot the positioner in its home
            # position.
            color = "gray"
            cb, fh = hw.fiber_position(fid, center)
            if real_shapes:
                plot_positioner(ax, patrol_rad, fid, center, cb, fh,
                                color=color, linewidth=linewidth)
            else:
                plot_positioner_simple(ax, patrol_rad, fid, center, cb, fh,
                                       color=color, linewidth=linewidth)

    return


plot_assignment_tile_file_hw = None


def plot_assignment_tile_file_initialize(hw):
    global plot_assignment_tile_file_hw
    plot_assignment_tile_file_hw = hw
    return


def plot_assignment_tile_file(fibers, real_shapes, params):
    (tile_id, tile_ra, tile_dec, infile, outfile) = params
    log = Logger.get()

    if os.path.isfile(outfile):
        log.info("Skipping existing plot {}".format(outfile))
        return
    else:
        log.info("Creating {}".format(outfile))

    header, fiber_data, targets_data, avail_data = read_assignment_fits_tile(
        (tile_id, infile))

    tavail = avail_table_to_dict(avail_data)

    log.debug("  tile {} at RA/DEC {} / {}".format(tile_id, tile_ra,
                                                   tile_dec))

    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_aspect("equal")

    # Target properties (x, y, color) for plotting
    tgs = plot_parse_table(targets_data)

    targetprops = plot_tile_targets_props(plot_assignment_tile_file_hw,
                                          tile_ra, tile_dec, tgs)

    log.debug("  tile {} has {} targets with properties"
              .format(tile_id, len(targets_data)))

    # When plotting available targets, we only consider those which have
    # RA/DEC information.  Depending on how the assignment was written out,
    # this might include only assigned targets or all targets available to
    # the tile.

    # Available targets for our selected fibers.
    avtg_fibers = [f for f in fibers if f in tavail]
    avtg_ids = np.unique([x for f in avtg_fibers for x in tavail[f]])

    # Downselect to include only targets with properties in the file.
    avtg = avtg_ids[np.isin(avtg_ids, targets_data["TARGETID"],
                            assume_unique=True)]

    plot_available(ax, targetprops, avtg, linewidth=0.1)

    # Assigned targets for our selected fibers
    tassign = {x["FIBER"]: x["TARGETID"] for x in fiber_data
               if (x["FIBER"] in fibers)}

    log.debug("  tile {} plotting {} assigned fibers"
              .format(tile_id, len(tassign)))

    fassign = {f: tassign[f] if f in tassign else -1 for f in fibers}

    plot_assignment(ax, plot_assignment_tile_file_hw, targetprops, fassign,
                    linewidth=0.1, real_shapes=real_shapes)

    ax.set_xlabel("Millimeters", fontsize="large")
    ax.set_ylabel("Millimeters", fontsize="large")
    plt.savefig(outfile, dpi=300, format="pdf")
    return


def plot_tiles(hw, tiles, result_dir=".", result_prefix="fiberassign",
               result_split_dir=False, plot_dir=".", plot_prefix="fiberassign",
               plot_split_dir=False, petals=None, real_shapes=False):
    """Plot assignment output.

    Args:
        hw (Hardware):  the hardware description.
        tiles (list):  List of tile IDs to process.
        result_dir (str):  Top-level directory of fiberassign results.
        result_prefix (str):  Prefix of each per-tile file name.
        result_split_dir (bool):  Results are in split tile directories.
        plot_dir (str):  Top-level directory for plots.
        plot_prefix (str):  Prefix of each per-tile output file name.
        plot_split_dir (bool):  Write outputs in split tile directories.
        petals (list):  List of petals to plot.
        real_shapes (bool):  If True, plot the full positioner shapes.

    Returns:
        None.

    """
    log = Logger.get()

    foundtiles = result_tiles(dir=result_dir, prefix=result_prefix)

    log.info("Found {} fiberassign tile files".format(len(foundtiles)))

    fibers = None
    if petals is None:
        fibers = [x for x in hw.fiber_id]
    else:
        fibers = list()
        for p in petals:
            fibers.extend([x for x in hw.petal_fibers[p]])
    fibers = np.array(fibers)

    plot_tile = partial(plot_assignment_tile_file, fibers, real_shapes)

    avail_tiles = np.array(tiles.id)
    select_tiles = [x for x in foundtiles if x in avail_tiles]

    tile_map_list = [(x, tiles.ra[tiles.order[x]], tiles.dec[tiles.order[x]],
                      result_path(x, dir=result_dir, prefix=result_prefix,
                                  split=result_split_dir),
                      result_path(x, dir=plot_dir, prefix=result_prefix,
                                  ext="pdf", create=True,
                                  split=plot_split_dir))
                     for x in select_tiles]

    log.info("Selecting {} fiberassign tile files".format(len(tile_map_list)))

    with mp.Pool(processes=default_mp_proc,
                 initializer=plot_assignment_tile_file_initialize,
                 initargs=(hw,)) as pool:
        pool.map(plot_tile, tile_map_list)

    return


def plot_qa_tile_color(desired, value, incr):
    des_color = "green"
    low_one_color = "gold"
    low_two_color = "red"
    low_color = "black"
    high_color = "cyan"
    if value == desired:
        return des_color
    if value > desired:
        return high_color
    if value < (desired - 2 * incr):
        return low_color
    if value < (desired - incr):
        return low_two_color
    return low_one_color


def plot_qa(data, outroot, outformat="pdf", labels=False):
    """Make plots of QA data.
    """
    hw = load_hardware()
    tile_radius = hw.focalplane_radius_deg

    fontpt = 1
    linewidth = 0.1

    fig = plt.figure(figsize=(12, 10))

    plot_param = [
        ("Total Fibers Assigned Per Tile", "assign_total", 5000, 5),
        ("Standards Assigned Per Tile", "assign_std", 100, 2),
        ("Sky Assigned Per Tile", "assign_sky", 400, 2),
    ]

    pindx = 1
    for title, key, desired, incr in plot_param:
        ax = fig.add_subplot(3, 1, pindx)
        ax.set_aspect("equal")
        xmin = 360.0
        xmax = 0.0
        ymin = 90.0
        ymax = -90.0
        for tid, props in data.items():
            xcent = props["tile_ra"]
            ycent = props["tile_dec"]
            if xcent > xmax:
                xmax = xcent
            if xcent < xmin:
                xmin = xcent
            if ycent > ymax:
                ymax = ycent
            if ycent < ymin:
                ymin = ycent
            color = plot_qa_tile_color(desired, props[key], incr)
            circ = plt.Circle((xcent, ycent), radius=tile_radius, fc="none",
                              ec=color, linewidth=linewidth)
            ax.add_artist(circ)
            if labels:
                ax.text(xcent, ycent, "{}".format(tid),
                        color=color, fontsize=fontpt,
                        horizontalalignment='center',
                        verticalalignment='center',
                        bbox=None)

        margin = 1.1 * tile_radius

        xmin -= margin
        xmax += margin
        ymin -= margin
        ymax += margin
        if xmin < 0.0:
            xmin = 0.0
        if xmax > 360.0:
            xmax = 360.0
        if ymin < -90.0:
            ymin = -90.0
        if ymax > 90.0:
            ymax = 90.0

        ax.set_xlim(left=xmin, right=xmax)
        ax.set_ylim(bottom=ymin, top=ymax)
        ax.set_xlabel("RA (degrees)", fontsize="large")
        ax.set_ylabel("DEC (degrees)", fontsize="large")
        ax.set_title(title)

        c_high = plot_qa_tile_color(desired, desired+1, incr)
        c_exact = plot_qa_tile_color(desired, desired, incr)
        c_low_one = plot_qa_tile_color(desired, desired-incr, incr)
        c_low_two = plot_qa_tile_color(desired, desired-2*incr, incr)
        c_low = plot_qa_tile_color(desired, 0, incr)

        c_low_two_val = desired - incr
        c_low_val = desired - 2 * incr

        legend_elements = [
            Patch(facecolor=c_high, edgecolor="none",
                  label="> {} assigned".format(desired)),
            Patch(facecolor=c_exact, edgecolor="none",
                  label="Exactly {} assigned".format(desired)),
            Patch(facecolor=c_low_one, edgecolor="none",
                  label="< {} assigned".format(desired)),
            Patch(facecolor=c_low_two, edgecolor="none",
                  label="< {} assigned".format(c_low_two_val)),
            Patch(facecolor=c_low, edgecolor="none",
                  label="< {} assigned".format(c_low_val)),
        ]
        ax.legend(handles=legend_elements, loc="best",
                  fontsize="x-small")
        pindx += 1

    outfile = "{}.{}".format(outroot, outformat)
    plt.savefig(outfile, dpi=300, format="pdf")

    return
