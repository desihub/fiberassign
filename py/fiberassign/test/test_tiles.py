"""
Test fiberassign tile operations.
"""

import os

import unittest

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles

from .simulate import sim_data_dir, sim_tiles


class TestTiles(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read(self):
        hw = load_hardware()
        tfile = os.path.join(sim_data_dir(), "footprint.fits")
        sim_tiles(tfile)
        tls = load_tiles(hw, tiles_file=tfile)
        print(tls)
        return
