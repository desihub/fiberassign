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

    def test_collision_xy(self):
        hw = load_hardware()
        center_mm = hw.center_mm
        fiber_id = hw.fiber_id
        nrot = 100
        rotrad = 0.5
        rotincr = 2 * np.pi / nrot
        tm = Timer()
        tm.start()
        for rot in range(nrot):
            xoff = rotrad * np.cos(rot * rotincr)
            yoff = rotrad * np.sin(rot * rotincr)
            xy = [(center_mm[p][0] + xoff, center_mm[p][1] + yoff)
                  for p in fiber_id]
            result = hw.check_collisions_xy(fiber_id, xy)
        tm.stop()
        tm.report("check_collisions_xy 100 configurations")
        return

    def test_collision_thetaphi(self):
        hw = load_hardware()
        fiber_id = hw.fiber_id
        ntheta = 10
        nphi = 10
        thetaincr = 2 * np.pi / ntheta
        phiincr = np.pi / nphi
        tm = Timer()
        tm.start()
        for thetarot in range(ntheta):
            for phirot in range(nphi):
                theta = [thetarot * thetaincr for x in fiber_id]
                phi = [phirot * phiincr for x in fiber_id]
                result = hw.check_collisions_thetaphi(fiber_id, theta, phi)
        tm.stop()
        tm.report("check_collisions_thetaphi 100 configurations")
        return
