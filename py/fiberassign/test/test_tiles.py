"""
Test fiberassign tile operations.
"""

import os

import unittest

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles

from .simulate import (
    test_subdir_create,
    sim_tiles,
    test_assign_date
)

from astropy.table import Table
from astropy.time import Time


class TestTiles(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read(self):
        test_dir = test_subdir_create("tiles_test_read")
        print('test_dir', test_dir)
        hw = load_hardware()
        tfile = os.path.join(test_dir, "footprint.fits")
        sfile = os.path.join(test_dir, "footprint_keep.txt")
        sim_tiles(tfile, selectfile=sfile)
        stiles = list()
        with open(sfile, "r") as f:
            for line in f:
                # Try to convert the first column to an integer.
                try:
                    stiles.append(int(line.split()[0]))
                except ValueError:
                    pass
        tls = load_tiles(tiles_file=tfile, select=stiles)
        print(tls)
        indx = 0
        for st in stiles:
            self.assertEqual(tls.order[st], indx)
            indx += 1

        #- Syntax / coverage tests
        tls = load_tiles(tiles_file=tfile, select=stiles, obstime='2020-10-20')
        self.assertEqual(tls.obstime[0], Time('2020-10-20T00:00:00').isot)
        tls = load_tiles(tiles_file=tfile, select=stiles, obstime='2020-10-20T10:20:30')
        self.assertEqual(tls.obstime[0], Time('2020-10-20T10:20:30').isot)

        #- including obsdate in tile file
        tiles = Table.read(tfile)
        tiles['OBSDATE'] = '2020-11-22'    # creates full column
        tiles['OBSDATE'][1] = '2020-11-23' # override second entry
        tfile = os.path.join(test_dir, "footprint-obsdate.fits")
        tiles.write(tfile)
        tls = load_tiles(tiles_file=tfile, select=stiles)
        self.assertEqual(tls.obstime[0], Time('2020-11-22').isot)
        self.assertEqual(tls.obstime[1], Time('2020-11-23').isot)

        return


def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
