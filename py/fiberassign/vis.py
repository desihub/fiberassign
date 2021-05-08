# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.vis
=======================

Visualization tools.

"""
from __future__ import absolute_import, division, print_function

import os
import warnings

import numpy as np

import multiprocessing as mp
from functools import partial

import fitsio

from ._internal import Shape

from .utils import Logger, default_mp_proc

from .hardware import load_hardware, FIBER_STATE_STUCK, FIBER_STATE_BROKEN

from .tiles import load_tiles

from .targets import (Targets, load_target_table,
                      TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                      TARGET_TYPE_SUPPSKY,
                      TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE)

from .assign import (read_assignment_fits_tile, result_tiles, result_path,
                     avail_table_to_dict, get_parked_thetaphi)

plt = None

def set_matplotlib_pdf_backend():
    """Set the matplotlib backend to PDF.

    This is necessary to render high resolution figures.
    """
    global plt
    if plt is not None:
        return
    try:
        import matplotlib
        matplotlib.use("pdf")
        import matplotlib.pyplot as plt
    except ValueError:
        warnings.warn(
            """Couldn't set the PDF matplotlib backend,
positioner plots may be low resolution.
Proceeding with the default matplotlib backend."""
        )
        import matplotlib.pyplot as plt


def plot_target_type_color(tgtype):
    color = "gray"
    tp = int(tgtype)
    if (tp & TARGET_TYPE_SAFE) != 0:
        color = "black"
    elif (tp & TARGET_TYPE_SKY) != 0:
        color = "blue"
    elif (tp & TARGET_TYPE_SUPPSKY) != 0:
        color = "blue"
    elif (tp & TARGET_TYPE_STANDARD) != 0:
        color = "gold"
        if (tp & TARGET_TYPE_SCIENCE) != 0:
            color = "green"
    elif (tp & TARGET_TYPE_SCIENCE) != 0:
        color = "red"
    return color


def plot_positioner(ax, patrol_rad, loc, center, shptheta, shpphi, color="k",
                    linewidth=0.2):
    """Plot one fiber positioner.
    """
    set_matplotlib_pdf_backend()
    patrol = plt.Circle((center[0], center[1]), radius=patrol_rad, fc=color,
                        ec="none", alpha=0.1)
    ax.add_artist(patrol)
    # Plot the arm from the center to the body
    thetacent = shptheta.axis
    armwidth = 0.25
    armlen = np.sqrt((thetacent[1] - center[1])**2
                     + (thetacent[0] - center[0])**2)
    armang = np.arctan2(thetacent[1] - center[1], thetacent[0] - center[0])
    sinarm = np.sin(armang)
    cosarm = np.cos(armang)
    arm_xoff = center[0] + (0.5*armwidth) * sinarm
    arm_yoff = center[1] - (0.5*armwidth) * cosarm
    armang_deg = armang * 180.0 / np.pi
    arm = plt.Rectangle((arm_xoff, arm_yoff), armlen, armwidth,
                        angle=armang_deg, color=color, linewidth=2*linewidth,
                        fill=False)
    ax.add_artist(arm)
    for piece in [shptheta, shpphi]:
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
    fontpt = 2.0
    xtxt = center[0] - 2 * armwidth * cosarm
    ytxt = center[1] - 2 * armwidth * sinarm
    ax.text(xtxt, ytxt, "{}".format(loc),
            color='k', fontsize=fontpt,
            horizontalalignment='center',
            verticalalignment='center',
            bbox=None)
    #        bbox=dict(fc='w', ec='none', pad=1, alpha=1.0))
    return


def plot_positioner_simple(ax, patrol_rad, loc, center, theta_ang, theta_arm,
                           phi_ang, phi_arm, color="k", linewidth=0.2):
    """Plot one fiber positioner.

    This uses a simpler representation of the positioner geometry, in order to
    speed up the plotting.

    """
    set_matplotlib_pdf_backend()
    patrol = plt.Circle((center[0], center[1]), radius=patrol_rad, fc=color,
                        ec="none", alpha=0.1)
    ax.add_artist(patrol)

    # Plot the arm from the center to the phi body
    theta_x = theta_arm * np.cos(theta_ang) + center[0]
    theta_y = theta_arm * np.sin(theta_ang) + center[1]

    ax.plot([center[0], theta_x], [center[1], theta_y], color=color,
            linewidth=5*linewidth)

    # Plot the phi arm.
    phi_x = phi_arm * np.cos(phi_ang + theta_ang) + theta_x
    phi_y = phi_arm * np.sin(phi_ang + theta_ang) + theta_y

    ax.plot([theta_x, phi_x], [theta_y, phi_y], color=color,
            linewidth=linewidth)

    fontpt = 2.0
    xtxt = center[0]
    ytxt = center[1] + 0.5
    ax.text(xtxt, ytxt, "{}".format(loc),
            color='k', fontsize=fontpt,
            horizontalalignment='center',
            verticalalignment='center',
            bbox=None)
    return


def plot_tile_targets_props(hw, tile_ra, tile_dec, tile_theta, tgs,
                            avail_tgid=None):
    if avail_tgid is None:
        avail_tgid = tgs.ids()
    ra = np.full(len(avail_tgid), 9999.9, dtype=np.float64)
    dec = np.full(len(avail_tgid), 9999.9, dtype=np.float64)
    color = list()
    for idx, tgid in enumerate(avail_tgid):
        tg = tgs.get(tgid)
        ra[idx] = tg.ra
        dec[idx] = tg.dec
        color.append(plot_target_type_color(tg.type))
    # We disable threading here, since it does not interact well with
    # multiprocessing.

    tgxy = hw.radec2xy_multi(tile_ra, tile_dec, tile_theta, ra, dec, False, 1)
    props = {tgid: {"xy": xy, "color": cl} for tgid, xy, cl
             in zip(avail_tgid, tgxy, color)}

    return props


def plot_available(ax, targetprops, selected, linewidth=0.1):
    mwidth = 5.0 * linewidth
    xdata = np.full(len(selected), 9999.9, dtype=np.float64)
    ydata = np.full(len(selected), 9999.9, dtype=np.float64)
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
    log = Logger.get()
    center_mm = hw.loc_pos_curved_mm
    theta_arm = hw.loc_theta_arm
    phi_arm = hw.loc_phi_arm
    theta_offset = hw.loc_theta_offset
    theta_min = hw.loc_theta_min
    theta_max = hw.loc_theta_max
    theta_pos = hw.loc_theta_pos
    phi_offset = hw.loc_phi_offset
    phi_min = hw.loc_phi_min
    phi_max = hw.loc_phi_max
    phi_pos = hw.loc_phi_pos
    state = hw.state
    loc_petal = dict(hw.loc_petal)
    device_type = dict(hw.loc_device_type)
    assigned = np.array(sorted(tile_assigned.keys()), dtype=np.int32)

    # Plot GFA / Petal edges.  Only plot one shape per petal, although
    # the code formally allows unique petal / GFA boundaries per device.

    if len(assigned) > 0:
        edge_gfa = dict()
        edge_petal = dict()
        for loc in assigned:
            pt = loc_petal[loc]
            if pt not in edge_gfa:
                edge_gfa[pt] = hw.loc_gfa_excl[loc]
                edge_petal[pt] = hw.loc_petal_excl[loc]
        for pt, shp in edge_gfa.items():
            for segs in shp.segments:
                xpts = np.array([p[0] for p in segs.points])
                ypts = np.array([p[1] for p in segs.points])
                ax.plot(xpts, ypts, linewidth=0.2*linewidth, color="gray")
        for pt, shp in edge_petal.items():
            for segs in shp.segments:
                xpts = np.array([p[0] for p in segs.points])
                ypts = np.array([p[1] for p in segs.points])
                ax.plot(xpts, ypts, linewidth=0.2*linewidth, color="gray")

    for lid in assigned:
        color = "gray"
        if (device_type[lid] != "POS") and (device_type[lid] != "ETC"):
            continue
        shptheta = Shape()
        shpphi = Shape()
        theta = None
        phi = None
        center = center_mm[lid]
        tgid = tile_assigned[lid][0]
        is_stuck_sky = tile_assigned[lid][1]
        patrol_rad = theta_arm[lid] + phi_arm[lid]
        failed = False
        is_assigned = (tgid >= 0)
        if is_assigned:
            # This fiber is assigned.  Plot the positioner located at the
            # assigned target.
            failed = hw.loc_position_xy(lid, targetprops[tgid]["xy"],
                                        shptheta, shpphi)
            if failed:
                msg = "Positioner at location {} cannot move to target {} at (x, y) = ({}, {}).  This should have been dected during assignment!".format(lid, tgid, targetprops[tgid]["xy"][0], targetprops[tgid]["xy"][1])
                log.warning(msg)
                raise RuntimeError(msg)
                is_assigned = False
                failed = False
            else:
                color = targetprops[tgid]["color"]
                theta, phi = hw.xy_to_thetaphi(
                    center, targetprops[tgid]["xy"],
                    theta_arm[lid], phi_arm[lid],
                    theta_offset[lid], phi_offset[lid],
                    theta_min[lid], phi_min[lid],
                    theta_max[lid], phi_max[lid],
                )
        if not is_assigned:
            # This fiber is unassigned.
            if (state[lid] & FIBER_STATE_STUCK) or (state[lid] & FIBER_STATE_BROKEN):
                # The positioner is stuck or fiber broken.  Plot it at its current
                # location.
                theta = theta_pos[lid] + theta_offset[lid]
                phi   = phi_pos  [lid] + phi_offset  [lid]
                msg = "Device location {}, state {} is stuck / broken, plotting fixed theta = {}, phi = {}".format(
                    lid, state[lid], theta, phi
                )
                log.debug(msg)
                failed = hw.loc_position_thetaphi(
                    lid, theta, phi, shptheta, shpphi, True
                )
                if is_stuck_sky:
                    color = "cyan"
            else:
                # Plot the positioner in its home (parked) position
                theta, phi = get_parked_thetaphi(theta_offset[lid],
                                                theta_min[lid], theta_max[lid],
                                                phi_offset[lid],
                                                phi_min[lid], phi_max[lid])
                msg = "Device location {}, state {} is unassigned, plotting parked theta = {}, phi = {}".format(
                    lid, state[lid], theta, phi
                )
                log.debug(msg)
                failed = hw.loc_position_thetaphi(lid, theta, phi, shptheta, shpphi, True)
            if failed:
                msg = "Positioner at location {} cannot move to its stuck or home position.  This should never happen!".format(lid)
                log.warning(msg)
        if not failed:
            if real_shapes:
                plot_positioner(
                    ax, patrol_rad, lid, center, shptheta, shpphi,
                    color=color, linewidth=linewidth
                )
            else:
                plot_positioner_simple(
                    ax, patrol_rad, lid, center, theta, theta_arm[lid], phi,
                    phi_arm[lid], color=color, linewidth=linewidth
                )
    return


def plot_assignment_tile_file(petals, real_shapes, params):
    (infile, outfile) = params
    set_matplotlib_pdf_backend()
    log = Logger.get()

    if os.path.isfile(outfile):
        log.info("Skipping existing plot {}".format(outfile))
        return
    else:
        log.info("Creating {}".format(outfile))

    header, fiber_data, targets_data, avail_data, gfa_data = \
        read_assignment_fits_tile((infile))

    tile_id = int(header["TILEID"])
    tile_ra = float(header["TILERA"])
    tile_dec = float(header["TILEDEC"])
    tile_theta = float(header["FIELDROT"])

    run_date = header["FA_RUN"]

    hw = load_hardware(rundate=run_date)

    locs = None
    if petals is None:
        locs = [x for x in hw.locations]
    else:
        locs = list()
        for p in petals:
            locs.extend([x for x in hw.petal_locations[p]])
    locs = np.array(locs)

    tavail = avail_table_to_dict(avail_data)

    log.debug("  tile {} at RA/DEC {} / {} with rotation {}".format(
        tile_id, tile_ra, tile_dec, tile_theta)
    )

    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_aspect("equal")

    # Target properties (x, y, color) for plotting
    tgs = Targets()
    if "FA_SURV" in header:
        load_target_table(tgs, targets_data,
                          survey=str(header["FA_SURV"]).rstrip(),
                          typecol="FA_TYPE")
    else:
        load_target_table(tgs, targets_data)

    targetprops = plot_tile_targets_props(hw, tile_ra, tile_dec, tile_theta, tgs)

    log.debug("  tile {} has {} targets with properties"
              .format(tile_id, len(targets_data)))

    # When plotting available targets, we only consider those which have
    # RA/DEC information.  Depending on how the assignment was written out,
    # this might include only assigned targets or all targets available to
    # the tile.

    # Available targets for our selected fibers.
    avtg_locs = [f for f in locs if f in tavail]
    avtg_ids = np.unique([x for f in avtg_locs for x in tavail[f]])

    # Downselect to include only targets with properties in the file.
    avtg = avtg_ids[np.isin(avtg_ids, targets_data["TARGETID"],
                            assume_unique=True)]

    plot_available(ax, targetprops, avtg, linewidth=0.1)

    # Assigned targets for our selected fibers.  We handle the special case of fibers
    # being used as sky but not formally assigned to a target.
    tassign = {
        x["LOCATION"]: (x["TARGETID"], (x["FA_TYPE"] & TARGET_TYPE_SKY))
        for x in fiber_data if (x["LOCATION"] in locs)
    }

    log.debug("  tile {} plotting {} assigned fibers"
              .format(tile_id, len(tassign)))

    fassign = {f: tassign[f] if f in tassign else (-1, False) for f in locs}

    plot_assignment(
        ax,
        hw,
        targetprops,
        fassign,
        linewidth=0.1,
        real_shapes=real_shapes
    )

    ax.set_xlabel("Curved Focal Surface Millimeters", fontsize="large")
    ax.set_ylabel("Curved Focal Surface Millimeters", fontsize="large")
    plt.savefig(outfile, dpi=300, format="pdf")
    plt.close()
    return


def plot_tiles(files, petals=None, real_shapes=False, serial=False):
    """Plot assignment output.

    Args:
        files (list):  The list of fiberassign files.
        petals (list):  List of petals to plot.
        real_shapes (bool):  If True, plot the full positioner shapes.
        serial (bool):  If True, disable use of multiprocessing.

    Returns:
        None.

    """
    log = Logger.get()

    log.info("Plotting {} fiberassign tile files".format(len(files)))

    plot_tile = partial(plot_assignment_tile_file, petals, real_shapes)

    file_map_list = list()
    for f in files:
        d, base = os.path.split(f)
        root = f.split(".")[0]
        file_map_list.append((f, os.path.join(d, "{}.pdf".format(root))))

    if serial:
        for params in file_map_list:
            plot_tile(params)
    else:
        with mp.Pool(processes=default_mp_proc) as pool:
            pool.map(plot_tile, file_map_list)

    return


def plot_assignment_tile(hw, tgs, tile_id, tile_ra, tile_dec, tile_theta,
                         tile_assign, tile_avail=None, petals=None,
                         real_shapes=False, outfile=None, figsize=8):
    set_matplotlib_pdf_backend()
    # Get selected fibers
    locs = None
    if petals is None:
        locs = [x for x in hw.locations]
    else:
        locs = list()
        for p in petals:
            locs.extend([x for x in hw.petal_locations[p]])
    locs = np.array(locs)

    # Available targets for our selected fibers.
    avtg_locs = None
    avtg_ids = None
    if tile_avail is None:
        # Just plot assigned targets
        avtg_locs = [f for f in locs if f in tile_assign]
        avtg_ids = [tile_assign[f] for f in avtg_locs]
    else:
        # Plot all available targets
        avtg_locs = [f for f in locs if f in tile_avail]
        avtg_ids = np.unique([x for f in avtg_locs for x in tile_avail[f]])

    # Target properties
    targetprops = plot_tile_targets_props(hw, tile_ra, tile_dec, tile_theta,
                                          tgs, avail_tgid=avtg_ids)

    fig = plt.figure(figsize=(figsize, figsize))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_aspect("equal")

    plot_available(ax, targetprops, avtg_ids, linewidth=0.1)

    # Assigned targets for our selected fibers
    tassign = {x: tile_assign[x] for x in locs if x in tile_assign}

    fassign = {f: tassign[f] if f in tassign else -1 for f in locs}

    plot_assignment(ax, hw, targetprops, fassign,
                    linewidth=0.1, real_shapes=real_shapes)

    ax.set_xlabel("Curved Focal Surface Millimeters", fontsize="large")
    ax.set_ylabel("Curved Focal Surface Millimeters", fontsize="large")
    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile, dpi=300, format="pdf")
        plt.close()
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
    set_matplotlib_pdf_backend()
    # Imported here, to ensure that the backend has been set.
    from matplotlib.patches import Patch

    hw = load_hardware()
    tile_radius = hw.focalplane_radius_deg

    fontpt = 1
    linewidth = 0.1

    fig = plt.figure(figsize=(12, 10))

    plot_param = [
        ("Total Fibers Assigned Per Tile", ["assign_total"], 5000, 5),
        ("Standards Assigned Per Tile", ["assign_std"], 100, 2),
        ("Sky Assigned Per Tile", ["assign_sky", "assign_suppsky"], 400, 2),
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
            keytot = np.sum([props[x] for x in key])
            color = plot_qa_tile_color(desired, keytot, incr)
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
