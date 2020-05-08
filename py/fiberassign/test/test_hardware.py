"""
Test fiberassign target operations.
"""

import unittest

from datetime import datetime

import numpy as np

import desimodel.io as dmio

from fiberassign.utils import Timer

from fiberassign.hardware import load_hardware

from .simulate import test_assign_date


class TestHardware(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read(self):
        hw = load_hardware(rundate=test_assign_date)
        print(hw)
        locs = hw.locations
        cs5 = hw.loc_pos_cs5_mm
        curved = hw.loc_pos_curved_mm
        # for p in locs:
        #     print("{}: curved = [{}, {}] cs5 = [{}, {}]".format(
        #         p, curved[p][0], curved[p][1], cs5[p][0], cs5[p][1]
        #     ))
        return

    def test_collision_xy(self):
        hw = load_hardware(rundate=test_assign_date)
        center_mm = hw.loc_pos_curved_mm
        locs = hw.locations
        nrot = 100
        rotrad = 0.5
        rotincr = 2 * np.pi / nrot
        tm = Timer()
        tm.start()
        for rot in range(nrot):
            xoff = rotrad * np.cos(rot * rotincr)
            yoff = rotrad * np.sin(rot * rotincr)
            xy = [(center_mm[p][0] + xoff, center_mm[p][1] + yoff)
                  for p in locs]
            result = hw.check_collisions_xy(locs, xy, 0)
        tm.stop()
        tm.report("check_collisions_xy 100 configurations")
        return

    def test_collision_thetaphi(self):
        hw = load_hardware(rundate=test_assign_date)
        locs = hw.locations
        ntheta = 10
        nphi = 10
        thetaincr = 2 * np.pi / ntheta
        phiincr = np.pi / nphi
        tm = Timer()
        tm.start()
        for thetarot in range(ntheta):
            for phirot in range(nphi):
                theta = [thetarot * thetaincr for x in locs]
                phi = [phirot * phiincr for x in locs]
                result = hw.check_collisions_thetaphi(locs, theta, phi, 0)
        tm.stop()
        tm.report("check_collisions_thetaphi 100 configurations")
        return

    def test_thetaphi_range(self):
        # Function to test that all positioners can reach a circle of
        # targets at fixed distance from their center.
        def check_reachable(hrdw, radius, increments, log_fail=True):
            centers = hw.loc_pos_curved_mm
            theta_arms = hw.loc_theta_arm
            phi_arms = hw.loc_phi_arm
            theta_mins = hw.loc_theta_min
            theta_maxs = hw.loc_theta_max
            theta_offsets = hw.loc_theta_offset
            phi_mins = hw.loc_phi_min
            phi_maxs = hw.loc_phi_max
            phi_offsets = hw.loc_phi_offset

            n_failed = 0
            for loc in hw.locations:
                center = centers[loc]
                theta_arm = theta_arms[loc]
                phi_arm = phi_arms[loc]
                theta_min = theta_mins[loc]
                theta_max = theta_maxs[loc]
                theta_offset = theta_offsets[loc]
                phi_min = phi_mins[loc]
                phi_max = phi_maxs[loc]
                phi_offset = phi_offsets[loc]
                ang = np.arange(increments) / (2 * np.pi)
                test_x = radius * np.cos(ang) + center[0]
                test_y = radius * np.sin(ang) + center[1]
                for xy in zip(test_x, test_y):
                    result = hrdw.xy_to_thetaphi(
                        center,
                        xy,
                        theta_arm,
                        phi_arm,
                        theta_offset,
                        phi_offset,
                        theta_min,
                        phi_min,
                        theta_max,
                        phi_max
                    )
                    if result[0] is None or result[1] is None:
                        if log_fail:
                            print(
                                "loc {} at ({}, {}) cannot reach ({}, {})".format(
                                    loc, center[0], center[1], xy[0], xy[1]
                                ), flush=True
                            )
                        n_failed += 1
                        break
                    else:
                        if not log_fail:
                            # log success instead
                            print(
                                "loc {} at ({}, {}) to ({}, {}) with ({}, {})".format(
                                    loc, center[0], center[1],
                                    xy[0], xy[1],
                                    result[0]*180.0/np.pi, result[1]*180.0/np.pi
                                ), flush=True
                            )
            return n_failed

        # Test nominal focalplane
        hw = load_hardware(rundate=test_assign_date)
        failed = check_reachable(hw, 3.0, 100)
        if (failed > 0):
            print(
                "{} positioners failed to reach ring at 3mm from center".format(
                    failed
                ),
                flush=True
            )
            self.assertTrue(False)

        # Now we are going to artificially restrict the phi angle range and test that
        # we cannot access the outer areas of the patrol radius.
        runtime = datetime.strptime(test_assign_date, "%Y-%m-%dT%H:%M:%S")
        fp, exclude, state, tmstr = dmio.load_focalplane(runtime)

        limit_radius = 2.0
        open_limit = 2.0 * np.arcsin(0.5 * limit_radius / 3.0)
        phi_limit_min = (np.pi - open_limit) * 180.0 / np.pi

        new_min = phi_limit_min - np.array(fp["OFFSET_P"])
        fp["MIN_P"][:] = new_min

        hw = load_hardware(focalplane=(fp, exclude, state))
        failed = check_reachable(hw, 3.0, 100, log_fail=False)
        if (failed != len(hw.locations)):
            print(
                "{} positioners reached 3mm from center, despite restricted phi".format(
                    len(hw.locations) - failed
                ),
                flush=True
            )
            self.assertTrue(False)
        return

def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
