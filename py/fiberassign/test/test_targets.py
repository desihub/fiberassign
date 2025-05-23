"""
Test fiberassign target operations.
"""

import os

import unittest

import numpy as np

from desitarget.targetmask import desi_mask, mws_mask

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
                                 default_main_gaia_stdmask,
                                 Targets, TargetsAvailable,
                                 LocationsAvailable, targets_in_tiles, create_tagalong)

from .simulate import sim_data_subdir_create, sim_tiles, sim_targets


class TestTargets(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_available(self):
        test_dir = sim_data_subdir_create("targets_test_available")
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
        tagalong = create_tagalong(plate_radec=False)
        load_target_file(tgs, tagalong, input_mtl)
        load_target_file(tgs, tagalong, input_std)
        load_target_file(tgs, tagalong, input_sky)
        load_target_file(tgs, tagalong, input_suppsky)
        print(tgs)

        # Test access
        ids = tgs.ids()
        tt = tgs.get(ids[0])
        tt.subpriority = 0.99

        # Compute the targets available to each fiber for each tile.
        hw = load_hardware()
        tfile = os.path.join(test_dir, "footprint.fits")
        sim_tiles(tfile)
        tiles = load_tiles(tiles_file=tfile)
        # Precompute target positions
        tile_targetids, tile_x, tile_y, tile_xy_cs5 = targets_in_tiles(hw, tgs, tiles, tagalong)
        tgsavail = TargetsAvailable(hw, tiles, tile_targetids, tile_x, tile_y)

        # Free the tree
        del tile_targetids, tile_x, tile_y

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
            0,
            desi_mask["SKY"].mask,
            desi_mask["SUPP_SKY"].mask,
            desi_mask["IN_BRIGHT_OBJECT"].mask,
            ]
        mws_target = [
            0,
            0,
            mws_mask["GAIA_STD_FAINT"].mask,
            0,
            0,
            0,
        ]
        fbatype = np.array([
            TARGET_TYPE_SCIENCE,
            TARGET_TYPE_STANDARD,
            TARGET_TYPE_STANDARD,
            TARGET_TYPE_SKY,
            TARGET_TYPE_SUPPSKY,
            0
            ])
        result = desi_target_type(
            desi_target, mws_target, default_main_sciencemask(), default_main_stdmask(),
            default_main_skymask(), default_main_suppskymask(),
            default_main_safemask(), default_main_excludemask(),
            default_main_gaia_stdmask())
        self.assertTrue(np.all(result == fbatype))

        # Scalar inputs
        for i in range(len(desi_target)):
            result = desi_target_type(
                desi_target[i], mws_target[i],
                default_main_sciencemask(), default_main_stdmask(),
                default_main_skymask(), default_main_suppskymask(),
                default_main_safemask(), default_main_excludemask(),
                default_main_gaia_stdmask())
            self.assertEqual(result, fbatype[i])

        # Does excludemask work?
        mask = desi_mask["ELG"].mask
        mask2 = 0 # MWS_TARGET
        result = desi_target_type(
            mask, mask2, default_main_sciencemask(), default_main_stdmask(),
            default_main_skymask(), default_main_suppskymask(),
            default_main_safemask(), default_main_excludemask(),
            default_main_gaia_stdmask())
        self.assertEqual(result, TARGET_TYPE_SCIENCE)

        result = desi_target_type(
            mask, mask2, default_main_sciencemask(), default_main_stdmask(),
            default_main_skymask(), default_main_suppskymask(),
            default_main_safemask(), mask,
            default_main_gaia_stdmask())
        self.assertEqual(result, 0)

        result = desi_target_type(
            [mask, mask], [mask2, mask2], default_main_sciencemask(), default_main_stdmask(),
            default_main_skymask(), default_main_suppskymask(),
            default_main_safemask(), mask, default_main_gaia_stdmask())
        self.assertTrue(not np.any(result))


def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m fiberassign.test.test_targets
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
