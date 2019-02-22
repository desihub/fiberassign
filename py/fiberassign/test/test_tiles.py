"""
Test fiberassign tile operations.
"""

import os

import unittest

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles

from .simulate import test_subdir_create, sim_tiles


class TestTiles(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read(self):
        test_dir = test_subdir_create("tiles_test_read")
        hw = load_hardware()
        stiles = [1165, 18465, 16870]
        tfile = os.path.join(test_dir, "footprint.fits")
        sim_tiles(tfile)
        tls = load_tiles(tiles_file=tfile, select=stiles)
        print(tls)
        indx = 0
        for st in stiles:
            self.assertEqual(tls.order[st], indx)
            indx += 1
        return
