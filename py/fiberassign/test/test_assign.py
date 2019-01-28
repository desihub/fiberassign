"""
Test fiberassign target operations.
"""
import os

import shutil

import unittest

import json

import numpy as np

import fitsio

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles

from fiberassign.targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                                 TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                                 Targets, TargetsAvailable, TargetTree,
                                 FibersAvailable, load_target_file)

from fiberassign.assign import (Assignment, write_assignment_fits,
                                write_assignment_ascii, merge_results,
                                read_assignment_fits_tile)

from fiberassign.qa import qa_tiles

from fiberassign.vis import plot_tiles, plot_qa

from .simulate import (test_subdir_create, sim_tiles, sim_targets, sim_status)


class TestAssign(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_io(self):
        test_dir = test_subdir_create("assign_test_io")
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

        # Create a hierarchical triangle mesh lookup of the targets positions
        tree = TargetTree(tgs, 0.01)

        # Compute the targets available to each fiber for each tile.
        fstatus = os.path.join(test_dir, "fiberstatus.ecsv")
        sim_status(fstatus)
        hw = load_hardware(status_file=fstatus)
        tfile = os.path.join(test_dir, "footprint.fits")
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

        # Write out, merge, read back in and verify

        write_assignment_ascii(tiles, asgn, out_dir=test_dir,
                               out_prefix="test_io_ascii_")

        write_assignment_fits(tiles, asgn, out_dir=test_dir,
                              out_prefix="basic_", all_targets=False)

        write_assignment_fits(tiles, asgn, out_dir=test_dir,
                              out_prefix="full_", all_targets=True)

        plot_tiles(hw, tiles, result_dir=test_dir,
                   result_prefix="basic_", plot_dir=test_dir,
                   plot_prefix="basic_",
                   result_split_dir=False, petals=[0],
                   serial=True)

        plot_tiles(hw, tiles, result_dir=test_dir,
                   result_prefix="full_", plot_dir=test_dir,
                   plot_prefix="full_",
                   result_split_dir=False, petals=[0],
                   serial=True)

        target_files = [
            input_mtl,
            input_sky,
            input_std
        ]
        tile_ids = list(tiles.id)

        merge_results(target_files, tile_ids, result_dir=test_dir,
                      result_prefix="basic_", out_dir=test_dir,
                      out_prefix="basic_tile-", copy_fba=False)

        merge_results(target_files, tile_ids, result_dir=test_dir,
                      result_prefix="full_", out_dir=test_dir,
                      out_prefix="full_tile-", copy_fba=False)

        # Here we test reading with the standard reading function

        for tid in tile_ids:
            tdata = asgn.tile_fiber_target(tid)
            avail = tgsavail.tile_data(tid)
            # Check basic format
            infile = os.path.join(test_dir,
                                  "basic_tile-{:05d}.fits".format(tid))
            inhead, fiber_data, targets_data, avail_data, gfa_targets = \
                read_assignment_fits_tile((tid, infile))
            for fid, tgid, tgra, tgdec in zip(
                    fiber_data["FIBER"],
                    fiber_data["TARGETID"],
                    fiber_data["TARGET_RA"],
                    fiber_data["TARGET_DEC"]):
                if tgid >= 0:
                    self.assertEqual(tgid, tdata[fid])
                    props = tgs.get(tgid)
                    self.assertEqual(tgra, props.ra)
                    self.assertEqual(tgdec, props.dec)

            # Check full format
            infile = os.path.join(test_dir,
                                  "full_tile-{:05d}.fits".format(tid))
            inhead, fiber_data, targets_data, avail_data, gfa_targets = \
                read_assignment_fits_tile((tid, infile))
            for fid, tgid, tgra, tgdec in zip(
                    fiber_data["FIBER"],
                    fiber_data["TARGETID"],
                    fiber_data["TARGET_RA"],
                    fiber_data["TARGET_DEC"]):
                if tgid >= 0:
                    self.assertEqual(tgid, tdata[fid])
                    props = tgs.get(tgid)
                    self.assertEqual(tgra, props.ra)
                    self.assertEqual(tgdec, props.dec)

        # Now read the files directly with fitsio and verify against the input
        # target data.

        for tid in tile_ids:
            tdata = asgn.tile_fiber_target(tid)
            avail = tgsavail.tile_data(tid)
            # Check basic format
            infile = os.path.join(test_dir,
                                  "basic_tile-{:05d}.fits".format(tid))
            fdata = fitsio.FITS(infile, "r")
            fassign = fdata["FIBERASSIGN"].read()
            ftargets = fdata["TARGETS"].read()
            for fid, tgid, tgra, tgdec, tgsub, tgprior, tgobs in zip(
                    fassign["FIBER"],
                    fassign["TARGETID"],
                    fassign["TARGET_RA"],
                    fassign["TARGET_DEC"],
                    fassign["SUBPRIORITY"],
                    fassign["PRIORITY"],
                    fassign["OBSCONDITIONS"]):
                if tgid >= 0:
                    self.assertEqual(tgid, tdata[fid])
                    props = tgs.get(tgid)
                    self.assertEqual(tgra, props.ra)
                    self.assertEqual(tgdec, props.dec)
                    self.assertEqual(tgsub, props.subpriority)
                    self.assertEqual(tgprior, props.priority)
                    self.assertEqual(tgobs, props.obscond)
            for tgid, tgra, tgdec, tgsub, tgprior, tgobs in zip(
                    ftargets["TARGETID"],
                    ftargets["RA"],
                    ftargets["DEC"],
                    ftargets["SUBPRIORITY"],
                    ftargets["PRIORITY"],
                    ftargets["OBSCONDITIONS"]):
                props = tgs.get(tgid)
                self.assertEqual(tgra, props.ra)
                self.assertEqual(tgdec, props.dec)
                self.assertEqual(tgsub, props.subpriority)
                self.assertEqual(tgprior, props.priority)
                self.assertEqual(tgobs, props.obscond)

            # Check full format
            infile = os.path.join(test_dir,
                                  "full_tile-{:05d}.fits".format(tid))
            fdata = fitsio.FITS(infile, "r")
            fassign = fdata["FIBERASSIGN"].read()
            ftargets = fdata["TARGETS"].read()
            for fid, tgid, tgra, tgdec, tgsub, tgprior, tgobs in zip(
                    fassign["FIBER"],
                    fassign["TARGETID"],
                    fassign["TARGET_RA"],
                    fassign["TARGET_DEC"],
                    fassign["SUBPRIORITY"],
                    fassign["PRIORITY"],
                    fassign["OBSCONDITIONS"]):
                if tgid >= 0:
                    self.assertEqual(tgid, tdata[fid])
                    props = tgs.get(tgid)
                    self.assertEqual(tgra, props.ra)
                    self.assertEqual(tgdec, props.dec)
                    self.assertEqual(tgsub, props.subpriority)
                    self.assertEqual(tgprior, props.priority)
                    self.assertEqual(tgobs, props.obscond)
            for tgid, tgra, tgdec, tgsub, tgprior, tgobs in zip(
                    ftargets["TARGETID"],
                    ftargets["RA"],
                    ftargets["DEC"],
                    ftargets["SUBPRIORITY"],
                    ftargets["PRIORITY"],
                    ftargets["OBSCONDITIONS"]):
                props = tgs.get(tgid)
                self.assertEqual(tgra, props.ra)
                self.assertEqual(tgdec, props.dec)
                self.assertEqual(tgsub, props.subpriority)
                self.assertEqual(tgprior, props.priority)
                self.assertEqual(tgobs, props.obscond)

        plot_tiles(hw, tiles, result_dir=test_dir,
                   result_prefix="basic_tile-", plot_dir=test_dir,
                   plot_prefix="basic_tile-",
                   result_split_dir=False, petals=[0],
                   serial=True)

        plot_tiles(hw, tiles, result_dir=test_dir,
                   result_prefix="full_tile-", plot_dir=test_dir,
                   plot_prefix="full_tile-",
                   result_split_dir=False, petals=[0],
                   serial=True)

        return

    def test_full(self):
        test_dir = test_subdir_create("assign_test_full")
        np.random.seed(123456789)
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

        # Create a hierarchical triangle mesh lookup of the targets positions
        tree = TargetTree(tgs, 0.01)

        # Read hardware properties
        fstatus = os.path.join(test_dir, "fiberstatus.ecsv")
        sim_status(fstatus)
        hw = load_hardware(status_file=fstatus)
        tfile = os.path.join(test_dir, "footprint.fits")
        sim_tiles(tfile)
        tiles = load_tiles(hw, tiles_file=tfile)

        # Compute the targets available to each fiber for each tile.
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

        write_assignment_fits(tiles, asgn, out_dir=test_dir, all_targets=True)

        plot_tiles(hw, tiles, result_dir=test_dir, plot_dir=test_dir,
                   petals=[0])

        qa_tiles(hw, tiles, result_dir=test_dir)

        qadata = None
        with open(os.path.join(test_dir, "qa.json"), "r") as f:
            qadata = json.load(f)

        for tile, props in qadata.items():
            self.assertEqual(4495, props["assign_science"])
            self.assertEqual(100, props["assign_std"])
            self.assertEqual(400, props["assign_sky"])

        plot_qa(qadata, os.path.join(test_dir, "qa"), outformat="pdf",
                labels=True)

        return
