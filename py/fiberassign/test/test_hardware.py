"""
Test fiberassign target operations.
"""

import unittest

from datetime import datetime

import numpy as np

import desimodel.io as dmio

from fiberassign.utils import Timer

from fiberassign.hardware import load_hardware

from .simulate import sim_assign_date


class TestHardware(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read(self):
        hw = load_hardware(rundate=sim_assign_date)
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
        hw = load_hardware(rundate=sim_assign_date)
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
        hw = load_hardware(rundate=sim_assign_date)
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
        hw = load_hardware(rundate=sim_assign_date)
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
        runtime = datetime.strptime(sim_assign_date, "%Y-%m-%dT%H:%M:%S%z")
        fp, exclude, state, tmstr = dmio.load_focalplane(runtime)

        # make a copy so that we aren't modifying the desimodel cache
        fp = fp.copy()

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

    def test_thetaphi_xy(self):
        # Test round trip consistency.
        def check_positioner(hrdw, radius, increments, log_fail=True):
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
                    thetaphi = hrdw.xy_to_thetaphi(
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
                    if thetaphi[0] is None or thetaphi[1] is None:
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
                                    thetaphi[0]*180.0/np.pi, thetaphi[1]*180.0/np.pi
                                ), flush=True
                            )
                    result = hrdw.thetaphi_to_xy(
                        center,
                        thetaphi[0],
                        thetaphi[1],
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
                                "loc {} at ({}, {}) invalid angles ({}, {})".format(
                                    loc, center[0], center[1],
                                    thetaphi[0]*180.0/np.pi, thetaphi[1]*180.0/np.pi
                                ), flush=True
                            )
                        n_failed += 1
                        break
                    else:
                        if not log_fail:
                            # log success instead
                            print(
                                "loc {} at ({}, {}) angles ({}, {}) to ({}, {})".format(
                                    loc, center[0], center[1],
                                    thetaphi[0]*180.0/np.pi, thetaphi[1]*180.0/np.pi,
                                    result[0], result[1]
                                ), flush=True
                            )
                    if not np.allclose([xy[0], xy[1]], [result[0], result[1]]):
                        print(
                            "loc {} at ({}, {}) failed roundtrip ({}, {}) != ({}, {})".format(
                                loc, center[0], center[1],
                                xy[0], xy[1], result[0], result[1]
                            ), flush=True
                        )
                        n_failed += 1
            return n_failed

        # Test nominal focalplane
        hw = load_hardware(rundate=sim_assign_date)
        failed = check_positioner(hw, 3.0, 100)
        if (failed > 0):
            print(
                "{} positioners failed X/Y roundtrip at 3mm from center".format(
                    failed
                ),
                flush=True
            )
            self.assertTrue(False)

        return
    
    def run_fba_compile_debug(self, tileid, tid, fiber):
        import os
        import fitsio
        import numpy as np
        from fiberassign.scripts.assign import parse_assign, run_assign_init
        from fiberassign.targets import targets_in_tiles, TargetsAvailable

        # tile settings
        tileidpad = "{:06d}".format(tileid)
        fafn = os.path.join(os.getenv("DESI_ROOT"), "target", "fiberassign", "tiles", "trunk", tileidpad[:3], "fiberassign-{}.fits.gz".format(tileidpad))
        hdr = fitsio.read_header(fafn, 0)
        rundate, fa_ha = hdr["RUNDATE"], hdr["FA_HA"]

        # some intermediate files
        fadir = os.path.join(os.getenv("DESI_ROOT"), "survey", "fiberassign", "main")
        tilesfn = os.path.join(fadir, tileidpad[:3], "{}-tiles.fits".format(tileidpad))
        targfn = os.path.join(fadir, tileidpad[:3], "{}-targ.fits".format(tileidpad))

        # is this tid is indeed reachable in the kpno file
        d = fitsio.read(fafn, "POTENTIAL_ASSIGNMENTS")
        kpno = ((d["TARGETID"] == tid) & (d["FIBER"] == fiber)).sum() == 1

        # compute the targets reachable by the fibers
        opts = ["--targets", targfn, "--footprint", tilesfn, "--rundate", rundate, "--ha", str(fa_ha)]
        ag = parse_assign(opts)
        hw, tiles, tgs, tagalong = run_assign_init(ag, plate_radec=True)
        tile_targetids, tile_x, tile_y, tile_xy_cs5 = targets_in_tiles(hw, tgs, tiles, tagalong)

        # cut on our tid for speed up
        ii = np.where(tile_targetids[tileid] == tid)[0]
        tile_targetids[tileid] = tile_targetids[tileid][ii]
        tile_x[tileid] = tile_x[tileid][ii]
        tile_y[tileid] = tile_y[tileid][ii]

        # compute available targets
        tgsavail = TargetsAvailable(hw, tiles, tile_targetids, tile_x, tile_y)

        # check for our tid
        loc = [loc for loc in hw.locations if hw.loc_fiber[loc] == fiber][0]
        nersc = tid in tgsavail.tile_data(tileid)[loc]
        return nersc

    def test_hardware_fiberassign_behavior_abs(self):
        """Test the behavior of the bug found in issue #258 where ::abs in
        hardware.cpp returned an int for r_min instead of a double in older
        compiler versions.  This is fixed to use fabs instead but the old
        buggy behavior is replicated by setting FIBERASSIGN_BEHAVIOR=0.
        """

        import os

        # non-reproducible case
        # that tid is reachable by that fiber at kpno, not for nersc rerun with compiling fiberassign with godesi/main
        tileid, tid, fiber = 3857, 39633537837566766, 1761

        # set DESIMODEL here, to avoid having to set it in the configuration
        os.environ["DESIMODEL"] = "/global/common/software/desi/perlmutter/desiconda/current/code/desimodel/main"
        if 'FIBERASSIGN_BEHAVIOR' in os.environ:
            current_behavior = os.environ['FIBERASSIGN_BEHAVIOR']
        else:
            current_behavior = -1

        print ("Testing 'buggy' behavior")
        os.environ['FIBERASSIGN_BEHAVIOR'] = "0"
        assert(self.run_fba_compile_debug(tileid, tid, fiber) is True)

        print ("Testing 'correct' behavior")
        os.environ['FIBERASSIGN_BEHAVIOR'] = "1"
        assert(self.run_fba_compile_debug(tileid, tid, fiber) is False)

        #reset behavior
        if (current_behavior != -1):
            #Unset is equivalent to 1 so leave it set
            os.environ['FIBERASSIGN_BEHAVIOR'] = current_behavior

        return True


def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
