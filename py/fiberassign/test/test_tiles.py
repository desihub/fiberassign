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

    # def test_read(self):
    #     test_dir = test_subdir_create("tiles_test_read")
    #     hw = load_hardware()
    #     tfile = os.path.join(test_dir, "footprint.fits")
    #     sfile = os.path.join(test_dir, "footprint_keep.txt")
    #     sim_tiles(tfile, selectfile=sfile)
    #     stiles = list()
    #     with open(sfile, "r") as f:
    #         for line in f:
    #             # Try to convert the first column to an integer.
    #             try:
    #                 stiles.append(int(line.split()[0]))
    #             except ValueError:
    #                 pass
    #     tls = load_tiles(tiles_file=tfile, select=stiles)
    #     print(tls)
    #     indx = 0
    #     for st in stiles:
    #         self.assertEqual(tls.order[st], indx)
    #         indx += 1
    #     return


def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
