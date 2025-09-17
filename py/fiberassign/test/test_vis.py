"""
Test fiberassign tile operations.
"""
import os

import unittest

from datetime import datetime

import numpy as np

from fiberassign.hardware import load_hardware

from fiberassign.vis import (
    plot_positioner,
    plot_positioner_simple,
    Shape,
    set_matplotlib_pdf_backend
)

from .simulate import sim_data_subdir_create, sim_focalplane, sim_assign_date

set_matplotlib_pdf_backend()
import matplotlib.pyplot as plt


class TestVis(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def _load_and_plotpos(self, time, dir, suffix, simple=False):
        hw = load_hardware(rundate=time)
        locs = hw.locations
        center_mm = hw.loc_pos_curved_mm
        theta_arm = hw.loc_theta_arm
        phi_arm = hw.loc_phi_arm
        theta_offset = hw.loc_theta_offset
        theta_min = hw.loc_theta_min
        theta_max = hw.loc_theta_max
        phi_offset = hw.loc_phi_offset
        phi_min = hw.loc_phi_min
        phi_max = hw.loc_phi_max

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_aspect("equal")

        # Compute the font size to use for detector labels
        figdpi = 75
        fontpix = 0.2 * figdpi
        fontpt = int(0.75 * fontpix)

        # angle increments
        nincr = 8
        configincr = 2.0 * np.pi / nincr

        # Plot the first fiber in a variety of positions
        lid = locs[0]
        patrol_mm = theta_arm[lid] + phi_arm[lid]
        center = center_mm[lid]

        # Plot data range in mm
        width = 1.2 * (2.0 * patrol_mm)
        height = 1.2 * (2.0 * patrol_mm)

        linewidth = 0.1

        # phipos = [np.pi/2.0]
        # phicol = ["r"]
        phipos = [0.0, np.pi/6.0]
        phicol = ["r", "b"]

        for configindx, (angphi, col) in enumerate(zip(phipos, phicol)):
            shptheta = Shape()
            shpphi = Shape()

            theta = theta_offset[lid] + theta_min[lid]
            phi = phi_offset[lid] + phi_min[lid] + angphi

            failed = hw.loc_position_thetaphi(
                lid, theta, phi, shptheta, shpphi
            )
            if failed:
                print(
                    "Failed to move positioner {} to theta = {}, phi = {}"
                    .format(lid, theta, phi)
                )
            else:
                if simple:
                    plot_positioner_simple(
                        ax, patrol_mm, lid, center, theta, theta_arm[lid], phi,
                        phi_arm[lid], color="black", linewidth=linewidth
                    )
                else:
                    plot_positioner(
                        ax, patrol_mm, lid, center, shptheta, shpphi,
                        color="black", linewidth=linewidth
                    )
            for inc in range(1, nincr):
                ang = inc * configincr
                effrad = 0.5 * patrol_mm * np.sin(angphi)
                theta = theta_offset[lid] + theta_min[lid] + ang
                xoff = effrad * np.cos(theta) + center[0]
                yoff = effrad * np.sin(theta) + center[1]
                failed = hw.loc_position_thetaphi(
                    lid, theta, phi, shptheta, shpphi
                )
                if failed:
                    print(
                        "Failed to move positioner {} to theta = {}, phi = {}"
                        .format(lid, theta, phi)
                    )
                else:
                    if simple:
                        plot_positioner_simple(
                            ax, patrol_mm, lid, center, theta, theta_arm[lid],
                            phi, phi_arm[lid], color=col, linewidth=linewidth
                        )
                    else:
                        plot_positioner(
                            ax, patrol_mm, lid, center, shptheta, shpphi,
                            color=col, linewidth=linewidth
                        )
                    xend = xoff
                    yend = yoff
                    ax.text(xend, yend, "{}".format(inc),
                            color='k', fontsize=fontpt,
                            horizontalalignment='center',
                            verticalalignment='center',
                            bbox=dict(fc='w', ec='none', pad=1, alpha=1.0))
        pxcent = center[0]
        pycent = center[1]
        half_width = 0.5 * width
        half_height = 0.5 * height
        ax.set_xlabel("Millimeters", fontsize="large")
        ax.set_ylabel("Millimeters", fontsize="large")
        ax.set_xlim([pxcent-half_width, pxcent+half_width])
        ax.set_ylim([pycent-half_height, pycent+half_height])
        outfile = os.path.join(dir, "test_plotpos_{}.pdf".format(suffix))
        plt.savefig(outfile, dpi=300, format="pdf")
        plt.close()

    def test_plotpos(self):
        test_dir = sim_data_subdir_create("vis_test_plotpos")
        time = sim_assign_date
        suffix = "{}_simple".format(time)
        self._load_and_plotpos(time, test_dir, suffix, simple=True)
        suffix = "{}".format(time)
        self._load_and_plotpos(time, test_dir, suffix, simple=False)
        time = "2012-12-12T00:00:00"
        suffix = "{}_simple".format(time)
        self._load_and_plotpos(time, test_dir, suffix, simple=True)
        suffix = "{}".format(time)
        self._load_and_plotpos(time, test_dir, suffix, simple=False)
        return


    def _load_and_plotfp(self, hw, dir, suffix, simple=False):
        locs = hw.locations
        center_mm = hw.loc_pos_curved_mm
        theta_arm = hw.loc_theta_arm
        phi_arm = hw.loc_phi_arm
        theta_offset = hw.loc_theta_offset
        theta_min = hw.loc_theta_min
        theta_max = hw.loc_theta_max
        phi_offset = hw.loc_phi_offset
        phi_min = hw.loc_phi_min
        phi_max = hw.loc_phi_max
        loc_petal = hw.loc_petal

        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_aspect("equal")

        linewidth = 0.1
        color = "black"

        plot_petals = [7]

        for lid in locs:
            if loc_petal[lid] not in plot_petals:
                continue
            shptheta = Shape()
            shpphi = Shape()
            center = center_mm[lid]
            # Many positioners can actually not reach to the center.  For
            # testing visualization, move theta to its minimum value and phi
            # to 180 degrees.
            theta = theta_offset[lid] + theta_min[lid]
            phi = phi_offset[lid] + np.pi
            # print("loc {} theta_0 {} + min {} = {}".format(lid, theta_offset[lid], theta_min[lid], theta))
            patrol_rad = theta_arm[lid] + phi_arm[lid]
            failed = hw.loc_position_thetaphi(
                lid, theta, phi, shptheta, shpphi
            )
            if failed:
                print(
                    "Failed to move positioner {} to theta = {}, phi = {}"
                    .format(lid, theta, phi)
                )
            else:
                if simple:
                    plot_positioner_simple(
                        ax, patrol_rad, lid, center, theta, theta_arm[lid],
                        phi, phi_arm[lid], color=color, linewidth=linewidth
                    )
                else:
                    plot_positioner(
                        ax, patrol_rad, lid, center, shptheta, shpphi,
                        color=color, linewidth=linewidth
                    )

        outfile = os.path.join(dir, "test_plotfp_{}.pdf".format(suffix))
        ax.set_xlabel("Millimeters", fontsize="large")
        ax.set_ylabel("Millimeters", fontsize="large")
        plt.savefig(outfile, dpi=300, format="pdf")
        plt.close()

    def test_plotfp(self):
        test_dir = sim_data_subdir_create("vis_test_plotfp")
        time = sim_assign_date
        hw = load_hardware(rundate=time)
        suffix = "{}_simple".format(time)
        self._load_and_plotfp(hw, test_dir, suffix, simple=True)
        suffix = "{}".format(time)
        self._load_and_plotfp(hw, test_dir, suffix, simple=False)
        # time = "2012-12-12T00:00:00"
        # hw = load_hardware(rundate=time)
        # suffix = "{}_simple".format(time)
        # self._load_and_plotfp(hw, test_dir, suffix, simple=True)
        # suffix = "{}".format(time)
        # self._load_and_plotfp(hw, test_dir, suffix, simple=False)
        return

    def test_plot_fakefp(self):
        test_dir = sim_data_subdir_create("vis_test_fakefp")
        time = sim_assign_date

        # Simulate a fake focalplane
        fp, exclude, state = sim_focalplane(rundate=sim_assign_date, fakepos=True)

        # Load the focalplane
        hw = load_hardware(focalplane=(fp, exclude, state))

        suffix = "{}".format(time)
        self._load_and_plotfp(hw, test_dir, suffix, simple=False)

        time = "2012-12-12T00:00:00"
        suffix = "{}".format(time)
        self._load_and_plotfp(hw, test_dir, suffix, simple=False)
        return
