"""
Test fiberassign target operations.
"""
import os

import shutil

import unittest

from datetime import datetime

import json

import numpy as np

import fitsio

import desimodel

from fiberassign.utils import option_list, GlobalTimers

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles, Tiles

from fiberassign.targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                                 TARGET_TYPE_SUPPSKY,
                                 TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                                 Targets, TargetsAvailable, TargetTree,
                                 LocationsAvailable, load_target_file)

from fiberassign.assign import (Assignment, write_assignment_fits,
                                write_assignment_ascii, merge_results,
                                read_assignment_fits_tile, run)

from fiberassign.qa import qa_tiles

from fiberassign.vis import plot_tiles, plot_qa

from fiberassign.scripts.assign import parse_assign, run_assign_full

from fiberassign.scripts.plot import parse_plot, run_plot

from fiberassign.scripts.qa import parse_qa, run_qa

from fiberassign.scripts.qa_plot import parse_plot_qa, run_plot_qa


from .simulate import (test_subdir_create, sim_tiles, sim_targets,
                       sim_focalplane, petal_rotation, test_assign_date)


class TestAssign(unittest.TestCase):

    def setUp(self):
        self.density_science = 5000
        self.density_standards = 5000
        self.density_sky = 100
        self.density_suppsky = 5000
        pass

    def tearDown(self):
        pass

    def test_io(self):
        test_dir = test_subdir_create("assign_test_io")
        input_mtl = os.path.join(test_dir, "mtl.fits")
        input_std = os.path.join(test_dir, "standards.fits")
        input_sky = os.path.join(test_dir, "sky.fits")
        input_suppsky = os.path.join(test_dir, "suppsky.fits")
        tgoff = 0
        nscience = sim_targets(
            input_mtl,
            TARGET_TYPE_SCIENCE,
            tgoff,
            density=self.density_science
        )
        tgoff += nscience
        nstd = sim_targets(
            input_std,
            TARGET_TYPE_STANDARD,
            tgoff,
            density=self.density_standards
        )
        tgoff += nstd
        nsky = sim_targets(
            input_sky,
            TARGET_TYPE_SKY,
            tgoff,
            density=self.density_sky
        )
        tgoff += nsky
        nsuppsky = sim_targets(
            input_suppsky,
            TARGET_TYPE_SUPPSKY,
            tgoff,
            density=self.density_suppsky
        )

        tgs = Targets()
        load_target_file(tgs, input_mtl)
        load_target_file(tgs, input_std)
        load_target_file(tgs, input_sky)
        load_target_file(tgs, input_suppsky)

        # Create a hierarchical triangle mesh lookup of the targets positions
        tree = TargetTree(tgs, 0.01)

        # Compute the targets available to each fiber for each tile.
        fp, exclude, state = sim_focalplane(rundate=test_assign_date)
        hw = load_hardware(focalplane=(fp, exclude, state))
        tfile = os.path.join(test_dir, "footprint.fits")
        sim_tiles(tfile)
        tiles = load_tiles(tiles_file=tfile)
        tgsavail = TargetsAvailable(hw, tgs, tiles, tree)

        # Free the tree
        del tree

        # Compute the fibers on all tiles available for each target
        favail = LocationsAvailable(tgsavail)

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

        plotpetals = [0]
        # plotpetals = None

        plot_tiles(hw, tiles, result_dir=test_dir,
                   result_prefix="basic_", plot_dir=test_dir,
                   plot_prefix="basic_",
                   result_split_dir=False, petals=plotpetals,
                   serial=True)

        plot_tiles(hw, tiles, result_dir=test_dir,
                   result_prefix="full_", plot_dir=test_dir,
                   plot_prefix="full_",
                   result_split_dir=False, petals=plotpetals,
                   serial=True)

        target_files = [
            input_mtl,
            input_sky,
            input_std
        ]
        tile_ids = list(tiles.id)

        merge_results(target_files, list(), tile_ids, result_dir=test_dir,
                      result_prefix="basic_", out_dir=test_dir,
                      out_prefix="basic_tile-", copy_fba=False)

        merge_results(target_files, list(), tile_ids, result_dir=test_dir,
                      result_prefix="full_", out_dir=test_dir,
                      out_prefix="full_tile-", copy_fba=False)

        # Here we test reading with the standard reading function

        for tid in tile_ids:
            tdata = asgn.tile_location_target(tid)
            avail = tgsavail.tile_data(tid)
            # Check basic format
            infile = os.path.join(test_dir,
                                  "basic_tile-{:06d}.fits".format(tid))
            inhead, fiber_data, targets_data, avail_data, gfa_targets = \
                read_assignment_fits_tile((tid, infile))

            for lid, tgid, tgra, tgdec in zip(
                    fiber_data["LOCATION"],
                    fiber_data["TARGETID"],
                    fiber_data["TARGET_RA"],
                    fiber_data["TARGET_DEC"]):
                if tgid >= 0:
                    self.assertEqual(tgid, tdata[lid])
                    props = tgs.get(tgid)
                    self.assertEqual(tgra, props.ra)
                    self.assertEqual(tgdec, props.dec)

            # Check full format
            infile = os.path.join(test_dir,
                                  "full_tile-{:06d}.fits".format(tid))
            inhead, fiber_data, targets_data, avail_data, gfa_targets = \
                read_assignment_fits_tile((tid, infile))
            for lid, tgid, tgra, tgdec in zip(
                    fiber_data["LOCATION"],
                    fiber_data["TARGETID"],
                    fiber_data["TARGET_RA"],
                    fiber_data["TARGET_DEC"]):
                if tgid >= 0:
                    self.assertEqual(tgid, tdata[lid])
                    props = tgs.get(tgid)
                    self.assertEqual(tgra, props.ra)
                    self.assertEqual(tgdec, props.dec)

        # Now read the files directly with fitsio and verify against the input
        # target data.

        for tid in tile_ids:
            tdata = asgn.tile_location_target(tid)
            avail = tgsavail.tile_data(tid)
            # Check basic format
            infile = os.path.join(test_dir,
                                  "basic_tile-{:06d}.fits".format(tid))
            fdata = fitsio.FITS(infile, "r")
            fassign = fdata["FIBERASSIGN"].read()
            ftargets = fdata["TARGETS"].read()
            for lid, tgid, tgra, tgdec, tgsub, tgprior, tgobs in zip(
                    fassign["LOCATION"],
                    fassign["TARGETID"],
                    fassign["TARGET_RA"],
                    fassign["TARGET_DEC"],
                    fassign["SUBPRIORITY"],
                    fassign["PRIORITY"],
                    fassign["OBSCONDITIONS"]):
                if tgid >= 0:
                    self.assertEqual(tgid, tdata[lid])
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
                                  "full_tile-{:06d}.fits".format(tid))
            fdata = fitsio.FITS(infile, "r")
            fassign = fdata["FIBERASSIGN"].read()
            ftargets = fdata["TARGETS"].read()
            for lid, tgid, tgra, tgdec, tgsub, tgprior, tgobs in zip(
                    fassign["LOCATION"],
                    fassign["TARGETID"],
                    fassign["TARGET_RA"],
                    fassign["TARGET_DEC"],
                    fassign["SUBPRIORITY"],
                    fassign["PRIORITY"],
                    fassign["OBSCONDITIONS"]):
                if tgid >= 0:
                    self.assertEqual(tgid, tdata[lid])
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
                   result_split_dir=False, petals=plotpetals,
                   serial=True)

        plot_tiles(hw, tiles, result_dir=test_dir,
                   result_prefix="full_tile-", plot_dir=test_dir,
                   plot_prefix="full_tile-",
                   result_split_dir=False, petals=plotpetals,
                   serial=True)
        return

    def test_full(self):
        test_dir = test_subdir_create("assign_test_full")
        np.random.seed(123456789)
        input_mtl = os.path.join(test_dir, "mtl.fits")
        input_std = os.path.join(test_dir, "standards.fits")
        input_sky = os.path.join(test_dir, "sky.fits")
        input_suppsky = os.path.join(test_dir, "suppsky.fits")
        tgoff = 0
        nscience = sim_targets(
            input_mtl,
            TARGET_TYPE_SCIENCE,
            tgoff,
            density=self.density_science
        )
        tgoff += nscience
        nstd = sim_targets(
            input_std,
            TARGET_TYPE_STANDARD,
            tgoff,
            density=self.density_standards
        )
        tgoff += nstd
        nsky = sim_targets(
            input_sky,
            TARGET_TYPE_SKY,
            tgoff,
            density=self.density_sky
        )
        tgoff += nsky
        nsuppsky = sim_targets(
            input_suppsky,
            TARGET_TYPE_SUPPSKY,
            tgoff,
            density=self.density_suppsky
        )

        tgs = Targets()
        load_target_file(tgs, input_mtl)
        load_target_file(tgs, input_std)
        load_target_file(tgs, input_sky)
        load_target_file(tgs, input_suppsky)

        # Create a hierarchical triangle mesh lookup of the targets positions
        tree = TargetTree(tgs, 0.01)

        # Read hardware properties
        fp, exclude, state = sim_focalplane(rundate=test_assign_date)
        hw = load_hardware(focalplane=(fp, exclude, state))
        tfile = os.path.join(test_dir, "footprint.fits")
        sim_tiles(tfile)
        tiles = load_tiles(tiles_file=tfile)

        # Compute the targets available to each fiber for each tile.
        tgsavail = TargetsAvailable(hw, tgs, tiles, tree)

        # Free the tree
        del tree

        # Compute the fibers on all tiles available for each target
        favail = LocationsAvailable(tgsavail)

        # Create assignment object
        asgn = Assignment(tgs, tgsavail, favail)

        run(asgn)

        write_assignment_fits(tiles, asgn, out_dir=test_dir, all_targets=True)

        plotpetals = [0]
        #plotpetals = None
        plot_tiles(hw, tiles, result_dir=test_dir, plot_dir=test_dir,
                   result_prefix="fba-",
                   real_shapes=True, petals=plotpetals, serial=True)

        qa_tiles(hw, tiles, result_dir=test_dir)

        qadata = None
        with open(os.path.join(test_dir, "qa.json"), "r") as f:
            qadata = json.load(f)

        for tile, props in qadata.items():
            self.assertTrue(props["assign_science"] >= 4485)
            self.assertEqual(100, props["assign_std"])
            self.assertTrue(
                (props["assign_sky"] + props["assign_suppsky"]) >= 400
            )

        plot_qa(qadata, os.path.join(test_dir, "qa"), outformat="pdf",
                labels=True)
        return

    def test_cli(self):
        test_dir = test_subdir_create("assign_test_cli")
        np.random.seed(123456789)
        input_mtl = os.path.join(test_dir, "mtl.fits")
        input_std = os.path.join(test_dir, "standards.fits")
        input_sky = os.path.join(test_dir, "sky.fits")
        input_suppsky = os.path.join(test_dir, "suppsky.fits")
        tgoff = 0
        nscience = sim_targets(
            input_mtl,
            TARGET_TYPE_SCIENCE,
            tgoff,
            density=self.density_science
        )
        tgoff += nscience
        nstd = sim_targets(
            input_std,
            TARGET_TYPE_STANDARD,
            tgoff,
            density=self.density_standards
        )
        tgoff += nstd
        nsky = sim_targets(
            input_sky,
            TARGET_TYPE_SKY,
            tgoff,
            density=self.density_sky
        )
        tgoff += nsky
        nsuppsky = sim_targets(
            input_suppsky,
            TARGET_TYPE_SUPPSKY,
            tgoff,
            density=self.density_suppsky
        )

        tfile = os.path.join(test_dir, "footprint.fits")
        sim_tiles(tfile)

        opts = {
            "targets": [input_mtl, input_std, input_sky, input_suppsky],
            "dir": test_dir,
            "footprint": tfile,
            "standards_per_petal": 10,
            "sky_per_petal": 40,
            "overwrite": True,
            "rundate": test_assign_date
        }
        optlist = option_list(opts)
        args = parse_assign(optlist)
        run_assign_full(args)

        plotpetals = "0"
        #plotpetals = "0,1,2,3,4,5,6,7,8,9"
        opts = {
            "footprint": tfile,
            "dir": test_dir,
            "petals": plotpetals,
            "serial": True,
            "rundate": test_assign_date
        }
        optlist = option_list(opts)
        args = parse_plot(optlist)
        run_plot(args)

        opts = {
            "dir": test_dir
        }
        optlist = option_list(opts)
        args = parse_qa(optlist)
        run_qa(args)

        opts = {
            "qafile": os.path.join(test_dir, "qa.json")
        }
        optlist = option_list(opts)
        args = parse_plot_qa(optlist)
        run_plot_qa(args)

        with open(os.path.join(test_dir, "qa.json"), "r") as f:
            qadata = json.load(f)

        for tile, props in qadata.items():
            self.assertTrue(props["assign_science"] >= 4490)
            self.assertEqual(100, props["assign_std"])
            self.assertTrue(
                (props["assign_sky"] + props["assign_suppsky"]) >= 400
            )
        return

    def test_fieldrot(self):
        test_dir = test_subdir_create("assign_test_fieldrot")
        np.random.seed(123456789)
        input_mtl = os.path.join(test_dir, "mtl.fits")
        input_std = os.path.join(test_dir, "standards.fits")
        input_sky = os.path.join(test_dir, "sky.fits")
        input_suppsky = os.path.join(test_dir, "suppsky.fits")
        tgoff = 0
        nscience = sim_targets(
            input_mtl,
            TARGET_TYPE_SCIENCE,
            tgoff,
            density=self.density_science
        )
        tgoff += nscience
        nstd = sim_targets(
            input_std,
            TARGET_TYPE_STANDARD,
            tgoff,
            density=self.density_standards
        )
        tgoff += nstd
        nsky = sim_targets(
            input_sky,
            TARGET_TYPE_SKY,
            tgoff,
            density=self.density_sky
        )
        tgoff += nsky
        nsuppsky = sim_targets(
            input_suppsky,
            TARGET_TYPE_SUPPSKY,
            tgoff,
            density=self.density_suppsky
        )

        # Simulate the tiles
        tfile = os.path.join(test_dir, "footprint.fits")
        sim_tiles(tfile)

        # petal mapping
        rotator = petal_rotation(1, reverse=False)

        rots = [0, 36]

        tile_ids = None

        for rt in rots:
            odir = "theta_{:02d}".format(rt)

            tgs = Targets()
            load_target_file(tgs, input_mtl)
            load_target_file(tgs, input_std)
            load_target_file(tgs, input_sky)
            load_target_file(tgs, input_suppsky)

            # Create a hierarchical triangle mesh lookup of the targets
            # positions
            tree = TargetTree(tgs, 0.01)

            # Manually override the field rotation
            tiles = load_tiles(tiles_file=tfile, obstheta=float(rt))
            if tile_ids is None:
                tile_ids = list(tiles.id)

            # Simulate a fake focalplane
            fp, exclude, state = sim_focalplane(rundate=test_assign_date, fakepos=True)

            # Load the focalplane
            hw = load_hardware(focalplane=(fp, exclude, state))

            # Compute the targets available to each fiber for each tile.
            tgsavail = TargetsAvailable(hw, tgs, tiles, tree)

            # Compute the fibers on all tiles available for each target
            favail = LocationsAvailable(tgsavail)

            # Create assignment object
            asgn = Assignment(tgs, tgsavail, favail)

            # First-pass assignment of science targets
            asgn.assign_unused(TARGET_TYPE_SCIENCE)

            out = os.path.join(test_dir, odir)

            write_assignment_fits(tiles, asgn, out_dir=out, all_targets=True)

            ppet = 6
            if odir == "theta_36":
                ppet = rotator[6]
            plot_tiles(hw, tiles, result_dir=out, plot_dir=out,
                       real_shapes=True, petals=[ppet], serial=True)

            # Explicitly free everything
            del asgn
            del favail
            del tgsavail
            del hw
            del tiles
            del tree
            del tgs

        # For each tile, compare the assignment output and verify that they
        # agree with a one-petal rotation.

        # NOTE:  The comparison below will NOT pass, since we are still
        # Sorting by highest priority available target and then (in case
        # of a tie) by fiber ID.  See line 333 of assign.cpp.  Re-enable this
        # test after that is changed to sort by location in case of a tie.

        # for tl in tile_ids:
        #     orig_path = os.path.join(
        #         test_dir, "theta_00", "fiberassign_{:06d}.fits".format(tl)
        #     )
        #     orig_header, orig_data, _, _, _ = \
        #         read_assignment_fits_tile((tl, orig_path))
        #     rot_path = os.path.join(
        #         test_dir, "theta_36", "fiberassign_{:06d}.fits".format(tl)
        #     )
        #     rot_header, rot_data, _, _, _ = \
        #         read_assignment_fits_tile((tl, rot_path))
        #     comppath = os.path.join(
        #         test_dir, "comp_00-36_{:06d}.txt".format(tl)
        #     )
        #     with open(comppath, "w") as fc:
        #         for dev, petal, tg in zip(
        #                 orig_data["DEVICE_LOC"], orig_data["PETAL_LOC"],
        #                 orig_data["TARGETID"]
        #         ):
        #             for newdev, newpetal, newtg in zip(
        #                     rot_data["DEVICE_LOC"], rot_data["PETAL_LOC"],
        #                     rot_data["TARGETID"]
        #             ):
        #                 rpet = rotator[newpetal]
        #                 if (newdev == dev) and (rpet == petal):
        #                     fc.write(
        #                         "{}, {} = {} : {}, {} = {}\n"
        #                         .format(petal, dev, tg, rpet, newdev, newtg)
        #                     )
        #                     # self.assertEqual(newtg, tg)

        return


def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
