"""
Test fiberassign target operations.
"""
import os
import subprocess
import re
import shutil
import unittest
from datetime import datetime
import json
import glob

import numpy as np

import fitsio

import desimodel

import fiberassign

from fiberassign.utils import option_list, GlobalTimers

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles, Tiles

from fiberassign.targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                                 TARGET_TYPE_SUPPSKY,
                                 TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                                 Targets, TargetsAvailable,
                                 LocationsAvailable, load_target_file, targets_in_tiles, create_tagalong)

from fiberassign.assign import (Assignment, write_assignment_fits,
                                write_assignment_ascii, merge_results)

from fiberassign.qa import qa_tiles, qa_targets

from fiberassign.vis import plot_tiles, plot_qa, set_matplotlib_pdf_backend

from fiberassign.scripts.assign import parse_assign, run_assign_full

from fiberassign.scripts.plot import parse_plot, run_plot

from fiberassign.scripts.qa import parse_qa, run_qa

from fiberassign.scripts.qa_plot import parse_plot_qa, run_plot_qa


from .simulate import (sim_data_subdir_create, sim_tiles, sim_targets,
                       sim_focalplane, petal_rotation, sim_assign_date)


class TestQA(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Find the location of scripts.  First try the case where we are running
        # tests from the top level of the source tree.
        cls.topDir = os.path.dirname( # top-level
            os.path.dirname( # build/
                os.path.dirname( # lib.arch/
                    os.path.dirname( # fiberassign/
                        os.path.dirname(os.path.abspath(__file__)) # test/
                        )
                    )
                )
            )
        cls.binDir = os.path.join(cls.topDir, "bin")
        if not os.path.isdir(cls.binDir):
            # We are running from some other directory from an installed package
            cls.topDir = os.path.dirname( # top-level
                os.path.dirname( # lib/
                    os.path.dirname( # python3.x/
                        os.path.dirname( # site-packages/
                            os.path.dirname( # egg/
                                os.path.dirname( # fiberassign/
                                    os.path.dirname(os.path.abspath(__file__)) # test/
                                )
                            )
                        )
                    )
                )
            )
            cls.binDir = os.path.join(cls.topDir, "bin")

    def setUp(self):
        self.density_science = 5000
        self.density_standards = 5000
        self.density_sky = 10
        self.density_suppsky = 5000
        pass

    def tearDown(self):
        pass

    def test_science(self):
        set_matplotlib_pdf_backend()
        import matplotlib.pyplot as plt
        test_dir = sim_data_subdir_create("qa_test_science")
        log_file = os.path.join(test_dir, "log.txt")

        np.random.seed(123456789)
        input_mtl = os.path.join(test_dir, "mtl.fits")
        # For this test, we will use just 2 science target classes, in order to verify
        # we get approximately the correct distribution
        sdist = [
            (3000, 1, 0.25, "QSO"),
            (2000, 1, 0.75, "ELG")
        ]
        nscience = sim_targets(
            input_mtl,
            TARGET_TYPE_SCIENCE,
            0,
            density=self.density_science,
            science_frac=sdist
        )

        log_msg = "Simulated {} science targets\n".format(nscience)

        tgs = Targets()
        tagalong = create_tagalong(plate_radec=False)
        load_target_file(tgs, tagalong, input_mtl)

        # Read hardware properties
        fp, exclude, state = sim_focalplane(rundate=sim_assign_date)
        hw = load_hardware(focalplane=(fp, exclude, state))
        tfile = os.path.join(test_dir, "footprint.fits")
        sim_tiles(tfile)
        tiles = load_tiles(tiles_file=tfile)

        # Precompute target positions
        tile_targetids, tile_x, tile_y, tile_xy_cs5 = targets_in_tiles(hw, tgs, tiles, tagalong)
        # Compute the targets available to each fiber for each tile.
        tgsavail = TargetsAvailable(hw, tiles, tile_targetids, tile_x, tile_y)

        # Compute the fibers on all tiles available for each target
        favail = LocationsAvailable(tgsavail)

        # Pass empty map of STUCK positioners that land on good sky
        stucksky = {}

        # Create assignment object
        asgn = Assignment(tgs, tgsavail, favail, stucksky)

        # First-pass assignment of science targets
        asgn.assign_unused(TARGET_TYPE_SCIENCE)

        # Redistribute
        asgn.redistribute_science()

        write_assignment_fits(tiles, tagalong, asgn, out_dir=test_dir, all_targets=True,
                              tile_xy_cs5=tile_xy_cs5)

        tile_ids = list(tiles.id)

        merge_results(
            [input_mtl], list(), tile_ids, result_dir=test_dir, copy_fba=False
        )

        # FIXME:  In order to use the qa_targets function, we need to know the
        # starting requested number of observations (NUMOBS_INIT).  Then we can use
        # that value for each target and compare to the number actually assigned.
        # However, the NUMOBS_INIT column was removed from the merged TARGET table.
        # If we are ever able to reach consensus on restoring that column, then we
        # can re-enable these tests below.
        #
        # qa_targets(
        #     hw,
        #     tiles,
        #     result_dir=test_dir,
        #     result_prefix="fiberassign-"
        # )
        #
        # # Load the target catalog so that we have access to the target properties
        #
        # fd = fitsio.FITS(input_mtl, "r")
        # scidata = np.array(np.sort(fd[1].read(), order="TARGETID"))
        # fd.close()
        # del fd
        #
        # # How many possible positioner assignments did we have?
        # nassign = 5000 * len(tile_ids)
        #
        # possible = dict()
        # achieved = dict()
        #
        # namepat = re.compile(r".*/qa_target_count_(.*)_init-(.*)\.fits")
        # for qafile in glob.glob("{}/qa_target_count_*.fits".format(test_dir)):
        #     namemat = namepat.match(qafile)
        #     name = namemat.group(1)
        #     obs = int(namemat.group(2))
        #     if obs == 0:
        #         continue
        #     fd = fitsio.FITS(qafile, "r")
        #     fdata = fd["COUNTS"].read()
        #     # Sort by target ID so we can select easily
        #     fdata = np.sort(fdata, order="TARGETID")
        #     tgid = np.array(fdata["TARGETID"])
        #     counts = np.array(fdata["NUMOBS_DONE"])
        #     avail = np.array(fdata["NUMOBS_AVAIL"])
        #     del fdata
        #     fd.close()
        #
        #     # Select target properties.  BOTH TARGET LISTS MUST BE SORTED.
        #     rows = np.where(np.isin(scidata["TARGETID"], tgid, assume_unique=True))[0]
        #
        #     ra = np.array(scidata["RA"][rows])
        #     dec = np.array(scidata["DEC"][rows])
        #     dtarget = np.array(scidata["DESI_TARGET"][rows])
        #     init = np.array(scidata["NUMOBS_INIT"][rows])
        #
        #     requested = obs * np.ones_like(avail)
        #
        #     under = np.where(avail < requested)[0]
        #     over = np.where(avail > requested)[0]
        #
        #     limavail = np.array(avail)
        #     limavail[over] = obs
        #
        #     deficit = np.zeros(len(limavail), dtype=np.int)
        #
        #     deficit[:] = limavail - counts
        #     deficit[avail == 0] = 0
        #
        #     possible[name] = np.sum(limavail)
        #     achieved[name] = np.sum(counts)
        #
        #     log_msg += "{}-{}:\n".format(name, obs)
        #
        #     pindx = np.where(deficit > 0)[0]
        #     poor_tgid = tgid[pindx]
        #     poor_dtarget = dtarget[pindx]
        #     log_msg += "  Deficit > 0: {}\n".format(len(poor_tgid))
        #     poor_ra = ra[pindx]
        #     poor_dec = dec[pindx]
        #     poor_deficit = deficit[pindx]
        #
        #     # Plot Target availability
        #     # Commented out by default, since in the case of high target density
        #     # needed for maximizing assignments, there are far more targets than
        #     # the number of available fiber placements.
        #
        #     # marksize = 4 * np.ones_like(deficit)
        #     #
        #     # fig = plt.figure(figsize=(12, 12))
        #     # ax = fig.add_subplot(1, 1, 1)
        #     # ax.scatter(ra, dec, s=2, c="black", marker="o")
        #     # for pt, pr, pd, pdef in zip(poor_tgid, poor_ra, poor_dec, poor_deficit):
        #     #     ploc = plt.Circle(
        #     #         (pr, pd), radius=(0.05*pdef), fc="none", ec="red"
        #     #     )
        #     #     ax.add_artist(ploc)
        #     # ax.set_xlabel("RA", fontsize="large")
        #     # ax.set_ylabel("DEC", fontsize="large")
        #     # ax.set_title(
        #     #     "Target \"{}\": (min(avail, requested) - counts) > 0".format(
        #     #         name, obs
        #     #     )
        #     # )
        #     # #ax.legend(handles=lg, framealpha=1.0, loc="upper right")
        #     # plt.savefig(os.path.join(test_dir, "{}-{}_deficit.pdf".format(name, obs)), dpi=300, format="pdf")
        #
        # log_msg += \
        #     "Assigned {} tiles for total of {} possible target observations\n".format(
        #         len(tile_ids), nassign
        #     )
        # ach = 0
        # for nm in possible.keys():
        #     ach += achieved[nm]
        #     log_msg += \
        #         "  type {} had {} possible target obs and achieved {}\n".format(
        #             nm, possible[nm], achieved[nm]
        #         )
        # frac = 100.0 * ach / nassign
        # log_msg += \
        #     "  {} / {} = {:0.2f}% of fibers were assigned\n".format(
        #         ach, nassign, frac
        #     )
        # for nm in possible.keys():
        #     log_msg += \
        #         "  type {} had {:0.2f}% of achieved observations\n".format(
        #             nm, achieved[nm] / ach
        #         )
        # with open(log_file, "w") as f:
        #     f.write(log_msg)
        #
        # self.assertGreaterEqual(frac,  99.0)

        # Test if qa-fiberassign script runs without crashing
        script = os.path.join(self.binDir, "qa-fiberassign")
        if os.path.exists(script):
            fafiles = glob.glob(f"{test_dir}/fiberassign-*.fits")
            cmd = "{} --targets {}".format(script, " ".join(fafiles))
            err = subprocess.call(cmd.split())
            self.assertEqual(err, 0, f"FAILED ({err}): {cmd}")
        else:
            print(f"ERROR: didn't find {script}")


def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
