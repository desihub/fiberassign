"""
Test fiberassign target operations.
"""

import os

import unittest

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles

from fiberassign.targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                                 TARGET_TYPE_STANDARD, load_target_file,
                                 Targets, TargetTree, TargetsAvailable,
                                 FibersAvailable)

from .simulate import (test_subdir_create, sim_tiles, sim_targets)


class TestTargets(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_available(self):
        test_dir = test_subdir_create("targets_test_available")
        input_mtl = os.path.join(test_dir, "mtl.fits")
        input_std = os.path.join(test_dir, "standards.fits")
        input_sky = os.path.join(test_dir, "sky.fits")
        nscience = sim_targets(input_mtl, TARGET_TYPE_SCIENCE, 0)
        nstd = sim_targets(input_std, TARGET_TYPE_STANDARD, nscience)
        nsky = sim_targets(input_sky, TARGET_TYPE_SKY, (nscience + nstd))

        tgs = Targets()
        load_target_file(tgs, input_mtl)
        load_target_file(tgs, input_std)
        load_target_file(tgs, input_sky)
        print(tgs)

        # Create a hierarchical triangle mesh lookup of the targets positions
        tree = TargetTree(tgs, 0.01)

        # Compute the targets available to each fiber for each tile.
        hw = load_hardware()
        tfile = os.path.join(test_dir, "footprint.fits")
        sim_tiles(tfile)
        tiles = load_tiles(hw, tiles_file=tfile)
        tgsavail = TargetsAvailable(tgs, tiles, tree)

        # Free the tree
        del tree

        # Compute the fibers on all tiles available for each target
        favail = FibersAvailable(tgsavail)

        return
