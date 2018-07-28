"""
Test fiberassign tile operations.
"""

import os
import unittest

from pkg_resources import resource_filename

from fiberassign._internal import Tiles

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles


class TestTiles(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    #
    # def test_read(self):
    #     hw = load_hardware()
    #     tls = load_tiles(hw)
    #     print(tls)
    #     return
