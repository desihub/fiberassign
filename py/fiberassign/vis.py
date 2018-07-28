# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.vis
=======================

Visualization tools.

"""
from __future__ import absolute_import, division, print_function

import os
import re
import numpy as np

import matplotlib
# try:
#     matplotlib.use("module://mplcairo.base")
# except ModuleNotFoundError:
#     try:
#         matplotlib.use("cairo")
#     except ModuleNotFoundError:
#         matplotlib.use("svg")
matplotlib.use("svg")
import matplotlib.pyplot as plt

import multiprocessing as mp
from functools import partial

from .utils import Logger, default_mp_proc

from .hardware import load_hardware

from .targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                      TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE)

from .assign import read_assignment_fits_tile


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


def plot_positioner(ax, fiber, center, cb, fh, color="k", linewidth=0.2):
    """Plot one fiber positioner.
    """
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
    xtxt = center[0] - armwidth * cosarm
    ytxt = center[1] - armwidth * sinarm
    ax.text(xtxt, ytxt, "{}".format(fiber),
            color='k', fontsize=fontpt,
            horizontalalignment='center',
            verticalalignment='center',
            bbox=None)
    #        bbox=dict(fc='w', ec='none', pad=1, alpha=1.0))
    return


def plot_tile_targets_props(hw, tile_ra, tile_dec, tgs, tile_avail, fibers):
    alltgid = np.unique([x for f in fibers for x in tile_avail[f]])
    ra = np.empty(len(alltgid), dtype=np.float64)
    dec = np.empty(len(alltgid), dtype=np.float64)
    color = list()
    for idx, tgid in enumerate(alltgid):
        tg = tgs.get(tgid)
        ra[idx] = tg.ra
        dec[idx] = tg.dec
        color.append(plot_target_type_color(tg.type))
    tgxy = hw.radec2xy_multi(tile_ra, tile_dec, ra, dec)
    props = dict()
    for tgid, xy, cl in zip(alltgid, tgxy, color):
        props[tgid] = dict()
        props[tgid]["xy"] = xy
        props[tgid]["color"] = cl
    return props


def plot_tile_target_table_props(hw, tile_ra, tile_dec, tgidx, tgs,
                                 tile_avail, fibers):
    tempid = list()
    for f in fibers:
        if f in tile_avail:
            tempid.extend(tile_avail[f])
    alltgid = np.unique(tempid)
    ra = np.empty(len(alltgid), dtype=np.float64)
    dec = np.empty(len(alltgid), dtype=np.float64)
    color = list()
    for idx, tgid in enumerate(alltgid):
        tgrow = tgidx[tgid]
        tgentry = tgs[tgrow]
        ra[idx] = tgentry["RA"]
        dec[idx] = tgentry["DEC"]
        color.append(plot_target_type_color(tgentry["FBATYPE"]))
    tgxy = hw.radec2xy_multi(tile_ra, tile_dec, ra, dec)
    props = dict()
    for tgid, xy, cl in zip(alltgid, tgxy, color):
        props[tgid] = dict()
        props[tgid]["xy"] = xy
        props[tgid]["color"] = cl
    return props


def plot_available(ax, targetprops, linewidth=0.1):
    mwidth = 5.0 * linewidth
    xdata = np.empty(len(targetprops), dtype=np.float64)
    ydata = np.empty(len(targetprops), dtype=np.float64)
    color = list()
    for idx, tgid in enumerate(sorted(targetprops.keys())):
        xdata[idx] = targetprops[tgid]["xy"][0]
        ydata[idx] = targetprops[tgid]["xy"][1]
        color.append(targetprops[tgid]["color"])
    ax.scatter(xdata, ydata, color=color, marker="x",
               linewidth=linewidth, s=mwidth)
    return


def plot_assignment(ax, hw, targetprops, fibers, tile_assigned, linewidth=0.1):
    center_mm = hw.center_mm
    for fid in fibers:
        center = center_mm[fid]
        if fid in tile_assigned:
            # This fiber is assigned.  Plot the positioner located at the
            # assigned target.
            tgid = tile_assigned[fid]
            cb, fh = hw.fiber_position(fid, targetprops[tgid]["xy"])
            plot_positioner(ax, fid, center, cb, fh,
                            color=targetprops[tgid]["color"],
                            linewidth=linewidth)
        else:
            # This fiber is unassigned.  Plot the positioner in its home
            # position.
            color = "gray"
            cb, fh = hw.fiber_position(fid, center)
            plot_positioner(ax, fid, center, cb, fh, color=color,
                            linewidth=linewidth)
    return


def plot_assignment_tile_file(hw, outroot, fibers, params):
    (tile_id,) = params
    log = Logger()

    outformat = "svg"
    outfile = "{}_{:06d}.{}".format(outroot, tile_id, outformat)
    if os.path.isfile(outfile):
        log.info("Skipping existing plot {}".format(outfile))
        return
    else:
        log.info("Creating {}".format(outfile))

    header, tgindx, tgdata, tavail = read_assignment_fits_tile(outroot,
                                                               (tile_id,))
    tile_ra = header["TILE_RA"]
    tile_dec = header["TILE_DEC"]

    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_aspect("equal")

    tassign = {x["FIBER"]: x["TARGETID"] for x in tgdata if x["FIBER"] >= 0}

    # Target properties for plotting
    targetprops = plot_tile_target_table_props(hw, tile_ra, tile_dec, tgindx,
                                               tgdata, tavail, fibers)

    plot_available(ax, targetprops, linewidth=0.1)

    plot_assignment(ax, hw, targetprops, fibers, tassign, linewidth=0.1)

    ax.set_xlabel("Millimeters", fontsize="large")
    ax.set_ylabel("Millimeters", fontsize="large")
    plt.savefig(outfile, dpi=300, format=outformat)
    return


def plot_tiles(resultdir=".", result_prefix="fiberassign", tiles=None,
               petals=None):
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

    outroot = os.path.join(resultdir, result_prefix)

    hw = load_hardware()
    fiber_id = hw.fiber_id
    petal_fibers = hw.petal_fibers
    fibers = None
    if petals is None:
        fibers = [x for x in fiber_id]
    else:
        fibers = list()
        for p in petals:
            fibers.extend([x for x in petal_fibers[p]])

    plot_tile = partial(plot_assignment_tile_file, hw, outroot, fibers)

    tile_map_list = None
    if tiles is None:
        tile_map_list = [(x,) for x in alltiles]
    else:
        tile_map_list = [(x,) for x in tiles]
    log.info("Selecting {} fiberassign tile files".format(len(tile_map_list)))

    with mp.Pool(processes=default_mp_proc) as pool:
        pool.map(plot_tile, tile_map_list)

    return
