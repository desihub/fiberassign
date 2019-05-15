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
        locs = hw.locations
        center_mm = hw.loc_pos_xy_mm

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
        lid = locs[0]
        nincr = 8
        configincr = 2.0 * np.pi / nincr

        center = center_mm[lid]

        for configindx, (configrad, col) in \
                enumerate(zip([0.5*patrol_mm, patrol_mm], ["r", "b"])):
            fig = plt.figure(figsize=(xfigsize, yfigsize), dpi=figdpi)
            ax = fig.add_subplot(1, 1, 1)
            ax.set_aspect("equal")
            cb, fh = hw.loc_position(lid, center_mm[lid])
            plot_positioner(ax, patrol_mm, lid, center, cb, fh, color="k")
            for inc in range(nincr):
                ang = inc * configincr
                xoff = configrad * np.cos(ang) + center_mm[lid][0]
                yoff = configrad * np.sin(ang) + center_mm[lid][1]
                cb, fh = hw.loc_position(lid, (xoff, yoff))
                plot_positioner(ax, patrol_mm, lid, center, cb, fh, color=col)
                xend = xoff
                yend = yoff
                ax.plot([center_mm[lid][0], xend], [center_mm[lid][1], yend],
                        color="k", linewidth="0.5")
                ax.text(xend, yend, "{}".format(inc),
                        color='k', fontsize=fontpt,
                        horizontalalignment='center',
                        verticalalignment='center',
                        bbox=dict(fc='w', ec='none', pad=1, alpha=1.0))
            pxcent = center_mm[lid][0]
            pycent = center_mm[lid][1]
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

def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
