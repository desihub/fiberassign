"""
Test fiberassign tile operations.
"""

import os

import unittest

import numpy as np

from pkg_resources import resource_filename

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles

from fiberassign.targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                                 TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                                 Targets, TargetsAvailable, TargetTree,
                                 FibersAvailable, load_target_file)

from fiberassign.assign import Assignment

from fiberassign.vis import (plot_positioner, plot_available, plot_assignment,
                             plot_tile_targets_props)


import matplotlib.pyplot as plt


class TestVis(unittest.TestCase):

    def setUp(self):
        # Find test data
        self.has_data = True
        self.input_mtl = resource_filename("fiberassign", "test/data/mtl.fits")
        self.input_sky = resource_filename("fiberassign", "test/data/sky.fits")
        self.input_std = resource_filename("fiberassign",
                                           "test/data/standards-dark.fits")
        if not os.path.isfile(self.input_mtl):
            self.has_data = False
        if not os.path.isfile(self.input_sky):
            self.has_data = False
        if not os.path.isfile(self.input_std):
            self.has_data = False
        return

    def tearDown(self):
        pass

    def test_plotpos(self):
        hw = load_hardware()
        patrol_mm = hw.patrol_mm
        fiber_id = hw.fiber_id
        center_mm = hw.center_mm

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
            plot_positioner(ax, fid, center, cb, fh, color="k")
            for inc in range(nincr):
                ang = inc * configincr
                xoff = configrad * np.cos(ang) + center_mm[fid][0]
                yoff = configrad * np.sin(ang) + center_mm[fid][1]
                cb, fh = hw.fiber_position(fid, (xoff, yoff))
                plot_positioner(ax, fid, center, cb, fh, color=col)
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
            outfile = "test_plotpos_{}.png".format(configindx)
            plt.savefig(outfile)
            plt.close()
        return

    def test_plotassigned(self):
        if not self.has_data:
            print("Test data not found- skipping")
            return

        # Read hardware properties and tiles
        hw = load_hardware()
        fiber_id = hw.fiber_id
        fiber_petal = hw.fiber_petal

        tiles = load_tiles(hw)

        # Read targets
        tgs = Targets()

        load_target_file(tgs, self.input_sky)
        load_target_file(tgs, self.input_std)
        load_target_file(tgs, self.input_mtl)

        # Create a hierarchical triangle mesh lookup of the targets positions
        tree = TargetTree(tgs)

        # Compute the targets available to each fiber for each tile.
        tgsavail = TargetsAvailable(tgs, tiles, tree)

        # Free the tree
        del tree

        # Compute the fibers on all tiles available for each target and sky
        favail = FibersAvailable(tgsavail)

        # Create assignment object
        asgn = Assignment(tgs, tgsavail, favail)

        # First-pass assignment of science targets
        asgn.assign_unused(TARGET_TYPE_SCIENCE)

        # Redistribute science targets across available petals
        asgn.redistribute_science()

        # Assign standards to unused fibers, up to 10 per petal
        asgn.assign_unused(TARGET_TYPE_STANDARD, 10)

        # Assign sky to unused fibers, up to 40 per petal
        asgn.assign_unused(TARGET_TYPE_SKY, 40)

        # Force assignment of sufficient standards: 10 per petal
        asgn.assign_force(TARGET_TYPE_STANDARD, 10)

        # Force assignment of sufficient standards: 40 per petal
        asgn.assign_force(TARGET_TYPE_SKY, 40)

        # Assign sky to unused fibers (no maximum)
        asgn.assign_unused(TARGET_TYPE_SKY)

        # Assign safe location to unused fibers (no maximum).  There should
        # always be at least one safe location (i.e. "BAD_SKY") for each fiber.
        # So after this is run every fiber should be assigned to something.
        asgn.assign_unused(TARGET_TYPE_SAFE)

        # Plot data range in mm
        tile_radius_deg = hw.focalplane_radius_deg
        platescale = 200.0
        tile_radius_mm = platescale * tile_radius_deg

        width = 1.2 * (2.0 * tile_radius_mm)
        height = 1.2 * (2.0 * tile_radius_mm)
        # Plot size in inches
        xfigsize = 12
        yfigsize = 12
        figdpi = 300

        # Compute the font size to use for detector labels
        fontpix = 0.2 * figdpi
        fontpt = int(0.75 * fontpix)

        fig = plt.figure(figsize=(xfigsize, yfigsize), dpi=figdpi)
        ax = fig.add_subplot(1, 1, 1)
        ax.set_aspect("equal")

        tile_order = tiles.order
        tile_id = 1148
        tile_ra = tiles.ra[tile_order[tile_id]]
        tile_dec = tiles.dec[tile_order[tile_id]]

        petal = 3
        fibers = [x for x in fiber_id if fiber_petal[x] == petal]
        # fibers = [x for x in fiber_id]

        tavail = tgsavail.tile_data(tile_id)
        tassign = asgn.tile_fiber_target(tile_id)

        # Target properties for plotting
        targetprops = plot_tile_targets_props(hw, tile_ra, tile_dec, tgs,
                                              tavail, fibers)

        plot_available(ax, targetprops, linewidth=0.1)

        plot_assignment(ax, hw, targetprops, fibers, tassign, linewidth=0.1)

        # pxcent = 0.0
        # pycent = 0.0
        # half_width = 0.5 * width
        # half_height = 0.5 * height
        ax.set_xlabel("Millimeters", fontsize="large")
        ax.set_ylabel("Millimeters", fontsize="large")
        # ax.set_xlim([pxcent-half_width, pxcent+half_width])
        # ax.set_ylim([pycent-half_height, pycent+half_height])

        outfile = "test_plotassigned.svg"
        plt.savefig(outfile, dpi=figdpi, format="svg")

        # outfile = "test_plotassigned.png"
        # plt.savefig(outfile, dpi=figdpi, format="png")

        plt.close()

        return
