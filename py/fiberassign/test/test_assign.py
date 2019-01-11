"""
Test fiberassign target operations.
"""
import os

import shutil

import unittest

import json

import numpy as np

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles

from fiberassign.targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                                 TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                                 Targets, TargetsAvailable, TargetTree,
                                 FibersAvailable, load_target_file)

from fiberassign.assign import (Assignment, write_assignment_fits,
                                write_assignment_ascii)

from fiberassign.qa import qa_tiles

from fiberassign.vis import plot_tiles, plot_qa

from .simulate import (sim_data_dir, sim_tiles, sim_targets)


class TestAssign(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_io(self):
        input_mtl = os.path.join(sim_data_dir(), "mtl.fits")
        input_std = os.path.join(sim_data_dir(), "standards.fits")
        input_sky = os.path.join(sim_data_dir(), "sky.fits")
        nscience = sim_targets(input_mtl, TARGET_TYPE_SCIENCE, 0)
        nstd = sim_targets(input_std, TARGET_TYPE_STANDARD, nscience)
        nsky = sim_targets(input_sky, TARGET_TYPE_SKY, (nscience + nstd))

        tgs = Targets()
        load_target_file(tgs, input_mtl)
        load_target_file(tgs, input_std)
        load_target_file(tgs, input_sky)

        # Create a hierarchical triangle mesh lookup of the targets positions
        tree = TargetTree(tgs, 0.01)

        # Compute the targets available to each fiber for each tile.
        hw = load_hardware()
        tfile = os.path.join(sim_data_dir(), "footprint.fits")
        sim_tiles(tfile)
        tiles = load_tiles(hw, tiles_file=tfile)
        tgsavail = TargetsAvailable(tgs, tiles, tree)

        # Free the tree
        del tree

        # Compute the fibers on all tiles available for each target
        favail = FibersAvailable(tgsavail)

        # First pass assignment
        asgn = Assignment(tgs, tgsavail, favail)
        asgn.assign_unused(TARGET_TYPE_SCIENCE)

        # Write out
        write_assignment_ascii(tiles, asgn, out_dir=sim_data_dir(),
                               out_prefix="test_io_ascii")
        return

    def test_full(self):
        np.random.seed(123456789)
        input_mtl = os.path.join(sim_data_dir(), "mtl.fits")
        input_std = os.path.join(sim_data_dir(), "standards.fits")
        input_sky = os.path.join(sim_data_dir(), "sky.fits")
        nscience = sim_targets(input_mtl, TARGET_TYPE_SCIENCE, 0)
        nstd = sim_targets(input_std, TARGET_TYPE_STANDARD, nscience)
        nsky = sim_targets(input_sky, TARGET_TYPE_SKY, (nscience + nstd))

        tgs = Targets()
        load_target_file(tgs, input_mtl)
        load_target_file(tgs, input_std)
        load_target_file(tgs, input_sky)

        # Create a hierarchical triangle mesh lookup of the targets positions
        tree = TargetTree(tgs, 0.01)

        # Compute the targets available to each fiber for each tile.
        hw = load_hardware()
        tfile = os.path.join(sim_data_dir(), "footprint.fits")
        sim_tiles(tfile)
        tiles = load_tiles(hw, tiles_file=tfile)
        tgsavail = TargetsAvailable(tgs, tiles, tree)

        # Free the tree
        del tree

        # Compute the fibers on all tiles available for each target
        favail = FibersAvailable(tgsavail)

        # Create assignment object
        asgn = Assignment(tgs, tgsavail, favail)

        # First-pass assignment of science targets
        asgn.assign_unused(TARGET_TYPE_SCIENCE)

        # Redistribute science targets
        asgn.redistribute_science()

        # Assign standards, 10 per petal
        asgn.assign_unused(TARGET_TYPE_STANDARD, 10)
        asgn.assign_force(TARGET_TYPE_STANDARD, 10)

        # Assign sky to unused fibers, up to 40 per petal
        asgn.assign_unused(TARGET_TYPE_SKY, 40)
        asgn.assign_force(TARGET_TYPE_SKY, 40)

        # If there are any unassigned fibers, try to place them somewhere.
        asgn.assign_unused(TARGET_TYPE_SCIENCE)
        asgn.assign_unused(TARGET_TYPE_SKY)

        outdir = os.path.join(sim_data_dir(), "test_full")
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.makedirs(outdir)

        write_assignment_fits(tiles, asgn, out_dir=outdir,
                              out_prefix="fiberassign_")

        plotdir = "{}_plots".format(outdir)
        if os.path.isdir(plotdir):
            shutil.rmtree(plotdir)
        os.makedirs(plotdir)

        plot_tiles(hw, tiles, result_dir=outdir, result_prefix="fiberassign_",
                   plot_dir=plotdir, petals=[0])

        qa_tiles(hw, tiles, result_dir=outdir, result_prefix="fiberassign_")

        qadata = None
        with open(os.path.join(outdir, "qa.json"), "r") as f:
            qadata = json.load(f)

        for tile, props in qadata.items():
            self.assertEqual(4500, props["assign_science"])
            self.assertEqual(100, props["assign_std"])
            self.assertEqual(400, props["assign_sky"])

        plot_qa(qadata, "{}_qa".format(outdir), outformat="pdf",
                labels=True)

        return
