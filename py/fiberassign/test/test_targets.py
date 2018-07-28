"""
Test fiberassign target operations.
"""

import os
import unittest

from pkg_resources import resource_filename

from fiberassign._internal import (Targets, TargetTree, TargetsAvailable,
                                   FibersAvailable)

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles

from fiberassign.targets import load_target_file


class TestTargets(unittest.TestCase):

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

    # def test_read(self):
    #     testfile = "/global/cscratch1/sd/kisner/desi/fiberassign/input/mtl_large.fits"
    #     tgs = Targets()
    #     load_target_file(tgs, testfile)
    #     return


    # def test_available(self):
    #     input_mtl = resource_filename("fiberassign", "test/data/mtl.fits")
    #     input_sky = resource_filename("fiberassign", "test/data/sky.fits")
    #     input_std = resource_filename("fiberassign", "test/data/standards-dark.fits")
    #     #tgs = read_targets([input_mtl])
    #     tgs = read_targets([input_mtl, input_sky, input_std])
    #     print(tgs)
    #
    #     # Create a hierarchical triangle mesh lookup of the targets positions
    #     tree = TargetsTree(tgs, 0.01)
    #
    #     # Compute the targets available to each fiber for each tile.
    #     hw = load_hardware()
    #     tiles = read_tiles(hw)
    #     tgsavail = TargetsAvailable(tiles, tree)
    #     print("ASS: done with targets available", flush=True)
    #
    #     # Free the tree
    #     del tree
    #     print("ASS: free tree", flush=True)
    #
    #     # Compute the fibers on all tiles available for each target
    #     favail = FibersAvailable(tgsavail)
    #     print("ASS: done with fibers available", flush=True)
    #
    #     return
