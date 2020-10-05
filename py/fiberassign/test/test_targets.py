"""
Test fiberassign target operations.
"""

import os

import unittest

import numpy as np

from desitarget.targetmask import desi_mask

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles

from fiberassign.targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                                 TARGET_TYPE_SUPPSKY,
                                 TARGET_TYPE_STANDARD, load_target_file,
                                 desi_target_type, default_main_sciencemask,
                                 default_main_skymask, default_main_stdmask,
                                 default_main_suppskymask,
                                 default_main_safemask,
                                 default_main_excludemask,
                                 Targets, TargetTree, TargetsAvailable,
                                 LocationsAvailable)

from .simulate import (test_subdir_create, sim_tiles, sim_targets, test_assign_date)


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
        input_suppsky = os.path.join(test_dir, "suppsky.fits")
        tgoff = 0
        nscience = sim_targets(input_mtl, TARGET_TYPE_SCIENCE, tgoff)
        tgoff += nscience
        nstd = sim_targets(input_std, TARGET_TYPE_STANDARD, tgoff)
        tgoff += nstd
        nsky = sim_targets(input_sky, TARGET_TYPE_SKY, tgoff)
        tgoff += nsky
        nsuppsky = sim_targets(input_suppsky, TARGET_TYPE_SUPPSKY, tgoff)

        tgs = Targets()
        load_target_file(tgs, input_mtl)
        load_target_file(tgs, input_std)
        load_target_file(tgs, input_sky)
        load_target_file(tgs, input_suppsky)
        print(tgs)

        # Test access
        ids = tgs.ids()
        tt = tgs.get(ids[0])
        tt.ra += 1.0e-5
        tt.dec += 1.0e-5
        tt.subpriority = 0.99

        # Create a hierarchical triangle mesh lookup of the targets positions
        tree = TargetTree(tgs, 0.01)

        # Compute the targets available to each fiber for each tile.
        hw = load_hardware()
        tfile = os.path.join(test_dir, "footprint.fits")
        sim_tiles(tfile)
        tiles = load_tiles(tiles_file=tfile)
        tgsavail = TargetsAvailable(hw, tgs, tiles, tree)

        # Free the tree
        del tree

        # Compute the fibers on all tiles available for each target
        favail = LocationsAvailable(tgsavail)

        return

    def test_target_type(self):
        """
        test fiberassign.targets.desi_target_type()
        """
        # Array inputs
        desi_target = [
            desi_mask["ELG"].mask,
            desi_mask["STD_FAINT"].mask,
            desi_mask["SKY"].mask,
            desi_mask["SUPP_SKY"].mask,
            desi_mask["IN_BRIGHT_OBJECT"].mask,
            ]
        fbatype = np.array([
            TARGET_TYPE_SCIENCE,
            TARGET_TYPE_STANDARD,
            TARGET_TYPE_SKY,
            TARGET_TYPE_SUPPSKY,
            0
            ])
        result = desi_target_type(
            desi_target, default_main_sciencemask(), default_main_stdmask(),
            default_main_skymask(), default_main_suppskymask(),
            default_main_safemask(), default_main_excludemask())
        self.assertTrue(np.all(result == fbatype))

        # Scalar inputs
        for i in range(len(desi_target)):
            result = desi_target_type(
                desi_target[i],
                default_main_sciencemask(), default_main_stdmask(),
                default_main_skymask(), default_main_suppskymask(),
                default_main_safemask(), default_main_excludemask())
            self.assertEqual(result, fbatype[i])

        # Does excludemask work?
        mask = desi_mask["ELG"].mask
        result = desi_target_type(
            mask, default_main_sciencemask(), default_main_stdmask(),
            default_main_skymask(), default_main_suppskymask(),
            default_main_safemask(), default_main_excludemask())
        self.assertEqual(result, TARGET_TYPE_SCIENCE)

        result = desi_target_type(
            mask, default_main_sciencemask(), default_main_stdmask(),
            default_main_skymask(), default_main_suppskymask(),
            default_main_safemask(), mask)
        self.assertEqual(result, 0)

        result = desi_target_type(
            [mask, mask], default_main_sciencemask(), default_main_stdmask(),
            default_main_skymask(), default_main_suppskymask(),
            default_main_safemask(), mask)
        self.assertTrue(not np.any(result))


def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m fiberassign.test.test_targets
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
