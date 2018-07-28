"""
Test fiberassign target operations.
"""

import os
import unittest

from pkg_resources import resource_filename

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles

from fiberassign.targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                                 TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                                 Targets, TargetsAvailable, TargetTree,
                                 FibersAvailable, load_target_file)

from fiberassign.assign import (Assignment, write_assignment_fits,
                                write_assignment_ascii)


class TestAssign(unittest.TestCase):

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


    # def test_io(self):
    #     if not self.has_data:
    #         print("Test data not found- skipping")
    #         return
    #
    #     # Read hardware properties and tiles
    #     hw = load_hardware()
    #     tiles = load_tiles(hw)
    #     print("done with hardware + tiles", flush=True)
    #
    #     # Read targets
    #     tgs = Targets()
    #
    #     # FIXME:  stop using the force option once we have new test files with
    #     # correct bits.
    #
    #     load_target_file(tgs, self.input_sky, typeforce=TARGET_TYPE_SKY)
    #     load_target_file(tgs, self.input_std, typeforce=TARGET_TYPE_STANDARD)
    #     load_target_file(tgs, self.input_mtl, typeforce=TARGET_TYPE_SCIENCE)
    #
    #     # Create a hierarchical triangle mesh lookup of the targets positions
    #     tree = TargetTree(tgs)
    #
    #     # Compute the targets available to each fiber for each tile.
    #     tgsavail = TargetsAvailable(tgs, tiles, tree)
    #
    #     # Free the tree
    #     del tree
    #
    #     # Compute the fibers on all tiles available for each target and sky
    #     favail = FibersAvailable(tgsavail)
    #
    #     # First pass assignment
    #     asgn = Assignment(tgs, tgsavail, favail)
    #     asgn.assign_unused(TARGET_TYPE_SCIENCE)
    #
    #     # Write out
    #
    #     #mtls = [self.input_mtl, self.input_std, self.input_sky]
    #     write_assignment_ascii(tiles, asgn, outroot="out-test_io", single=True)
    #
    #     return


    # def test_full(self):
    #     # if not self.has_data:
    #     #     print("Test data not found- skipping")
    #     #     return
    #
    #     # Read hardware properties and tiles
    #     hw = load_hardware()
    #     tiles = load_tiles(hw)
    #
    #     tgs = Targets()
    #
    #     # FIXME:  stop using the force option once we have new test files with
    #     # correct bits.
    #
    #     # load_target_file(tgs, self.input_sky, typeforce=TARGET_TYPE_SKY)
    #     # load_target_file(tgs, self.input_std, typeforce=TARGET_TYPE_STANDARD)
    #     # load_target_file(tgs, self.input_mtl, typeforce=TARGET_TYPE_SCIENCE)
    #
    #     load_target_file(tgs, "/global/cscratch1/sd/kisner/desi/fiberassign/input/mtl_large.fits")
    #     load_target_file(tgs, "/global/cscratch1/sd/kisner/desi/fiberassign/input/std_large.fits")
    #     load_target_file(tgs, "/global/cscratch1/sd/kisner/desi/fiberassign/input/sky_large.fits")
    #
    #     # Create a hierarchical triangle mesh lookup of the targets positions
    #     tree = TargetTree(tgs)
    #
    #     # Compute the targets available to each fiber for each tile.
    #     tgsavail = TargetsAvailable(tgs, tiles, tree)
    #
    #     # Free the tree
    #     del tree
    #
    #     # Compute the fibers on all tiles available for each target and sky
    #     favail = FibersAvailable(tgsavail)
    #
    #     # Create assignment object
    #     asgn = Assignment(tgs, tgsavail, favail)
    #
    #     # First-pass assignment of science targets
    #     asgn.assign_unused(TARGET_TYPE_SCIENCE)
    #
    #     # outfile = "out-test_full-firstpass"
    #     # write_assignment_ascii(tiles, asgn, outroot=outfile, single=True)
    #
    #     # Redistribute science targets
    #     asgn.redistribute_science()
    #
    #     # outfile = "out-test_full-redist"
    #     # write_assignment_ascii(tiles, asgn, outroot=outfile, single=True)
    #
    #     # Assign standards to unused fibers
    #     asgn.assign_unused(TARGET_TYPE_STANDARD, 10)
    #
    #     # Assign sky to unused fibers
    #     asgn.assign_unused(TARGET_TYPE_SKY, 40)
    #     #
    #     # # outfile = "out-test_full-stdsky_unused"
    #     # # write_assignment_ascii(tiles, asgn, outroot=outfile, single=True)
    #     #
    #     # # Force assignment of sufficient standards
    #     # asgn.assign_force(TARGET_TYPE_STANDARD, 10)
    #     #
    #     # # Force assignment of sufficient standards
    #     # asgn.assign_force(TARGET_TYPE_SKY, 40)
    #     #
    #     # # outfile = "out-test_full-stdsky_force"
    #     # # write_assignment_ascii(tiles, asgn, outroot=outfile, single=True)
    #     #
    #     # # Assign sky to unused fibers
    #     # asgn.assign_unused(TARGET_TYPE_SKY)
    #     #
    #     # # Assign safe location to unused fibers
    #     # asgn.assign_unused(TARGET_TYPE_SAFE)
    #     #
    #     # # outfile = "out-test_full-final"
    #     # # write_assignment_ascii(tiles, asgn, outroot=outfile, single=True)
    #
    #     return
