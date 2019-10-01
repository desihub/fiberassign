"""
Test fiberassign target operations.
"""

import unittest

import numpy as np

from fiberassign.utils import Timer

from fiberassign.hardware import load_hardware


class TestHardware(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read(self):
        hw = load_hardware()
        print(hw)
        return

    # def test_collision_xy(self):
    #     hw = load_hardware()
    #     center_mm = hw.loc_pos_xy_mm
    #     locs = hw.locations
    #     nrot = 100
    #     rotrad = 0.5
    #     rotincr = 2 * np.pi / nrot
    #     tm = Timer()
    #     tm.start()
    #     for rot in range(nrot):
    #         xoff = rotrad * np.cos(rot * rotincr)
    #         yoff = rotrad * np.sin(rot * rotincr)
    #         xy = [(center_mm[p][0] + xoff, center_mm[p][1] + yoff)
    #               for p in locs]
    #         result = hw.check_collisions_xy(locs, xy, 0)
    #     tm.stop()
    #     tm.report("check_collisions_xy 100 configurations")
    #     return
    #
    # def test_collision_thetaphi(self):
    #     hw = load_hardware()
    #     locs = hw.locations
    #     ntheta = 10
    #     nphi = 10
    #     thetaincr = 2 * np.pi / ntheta
    #     phiincr = np.pi / nphi
    #     tm = Timer()
    #     tm.start()
    #     for thetarot in range(ntheta):
    #         for phirot in range(nphi):
    #             theta = [thetarot * thetaincr for x in locs]
    #             phi = [phirot * phiincr for x in locs]
    #             result = hw.check_collisions_thetaphi(locs, theta, phi, 0)
    #     tm.stop()
    #     tm.report("check_collisions_thetaphi 100 configurations")
    #     return

def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
