"""
Test fiberassign tile operations.
"""
import os

import unittest

import numpy as np

from fiberassign.hardware import load_hardware

from fiberassign.vis import (plot_positioner, )

from .simulate import test_subdir_create

import matplotlib.pyplot as plt


class TestVis(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_plotpos(self):
        test_dir = test_subdir_create("vis_test_plotpos")
        hw = load_hardware()
        patrol_mm = hw.patrol_mm
        fiber_id = hw.fiber_id
        center_mm = hw.fiber_pos_xy_mm

        # Plot data range in mm
        width = 1.2 * (2.0 * patrol_mm)
        height = 1.2 * (2.0 * patrol_mm)
        # Plot size in inches
        xfigsize = 8
        yfigsize = 8
        figdpi = 75

        # Compute the font size to use for detector labels
        fontpix = 0.2 * figdpi
        fontpt = int(0.75 * fontpix)

        # Plot the first fiber in a variety of positions
        fid = fiber_id[0]
        nincr = 8
        configincr = 2.0 * np.pi / nincr

        center = center_mm[fid]

        for configindx, (configrad, col) in \
                enumerate(zip([0.5*patrol_mm, patrol_mm], ["r", "b"])):
            fig = plt.figure(figsize=(xfigsize, yfigsize), dpi=figdpi)
            ax = fig.add_subplot(1, 1, 1)
            ax.set_aspect("equal")
            cb, fh = hw.fiber_position(fid, center_mm[fid])
            plot_positioner(ax, patrol_mm, fid, center, cb, fh, color="k")
            for inc in range(nincr):
                ang = inc * configincr
                xoff = configrad * np.cos(ang) + center_mm[fid][0]
                yoff = configrad * np.sin(ang) + center_mm[fid][1]
                cb, fh = hw.fiber_position(fid, (xoff, yoff))
                plot_positioner(ax, patrol_mm, fid, center, cb, fh, color=col)
                xend = xoff
                yend = yoff
                ax.plot([center_mm[fid][0], xend], [center_mm[fid][1], yend],
                        color="k", linewidth="0.5")
                ax.text(xend, yend, "{}".format(inc),
                        color='k', fontsize=fontpt,
                        horizontalalignment='center',
                        verticalalignment='center',
                        bbox=dict(fc='w', ec='none', pad=1, alpha=1.0))
            pxcent = center_mm[fid][0]
            pycent = center_mm[fid][1]
            half_width = 0.5 * width
            half_height = 0.5 * height
            ax.set_xlabel("Millimeters", fontsize="large")
            ax.set_ylabel("Millimeters", fontsize="large")
            ax.set_xlim([pxcent-half_width, pxcent+half_width])
            ax.set_ylim([pycent-half_height, pycent+half_height])

            outfile = os.path.join(test_dir,
                                   "test_plotpos_{}.png".format(configindx))
            plt.savefig(outfile)
            plt.close()
        return
