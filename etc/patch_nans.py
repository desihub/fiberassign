from astropy.io import fits # use astropy.io.fits here because it allows us to implace update a fits file.
import numpy as np
from pathlib import Path
from fiberassign.hardware import load_hardware, xy2radec, xy2cs5
from types import SimpleNamespace

import argparse
from datetime import datetime, timezone

parser = argparse.ArgumentParser()
parser.add_argument("input", type=str, help="base directory of copied fiberassign files to process")
args = parser.parse_args()
# We will mimic the loop that sets the RA/DEC values at save time
# using the run parameters saved in the original fiberassign files.

orig_loc = Path(args.input)

for fba_file in sorted(orig_loc.glob("*/fiberassign*fits.gz")):
    with fits.open(fba_file, "update") as h:
        tbl = h["FIBERASSIGN"].data
        # Table of objects that need to be patched.
        patch_tbl = tbl[np.isnan(tbl["TARGET_RA"]) | np.isnan(tbl["TARGET_DEC"])]

        print(f"{fba_file}, {len(patch_tbl)} to update")
        if len(patch_tbl) == 0: continue

        header = h[0].header
        runtime = header["FA_RUN"]

        tile_ra = header["TILERA"]
        tile_dec = header["TILEDEC"]
        tile_obstime = header["FA_PLAN"]
        tile_obsha = header["FA_HA"]
        tile_obstheta = header["FIELDROT"]

        hw = load_hardware(None, rundate=runtime)

        # To be 100% clear, I don't exactly understand why we pull those out into
        # new object. Again, I'm just mimicing the loop that sets these
        # values in fiberassign when the output is saved.
        hwduck = SimpleNamespace()
        for attribute in [
                'state', 'loc_pos_curved_mm', 'loc_theta_pos',
                'loc_theta_offset', 'loc_phi_pos', 'loc_phi_offset',
                'loc_theta_arm', 'loc_phi_arm', 'loc_theta_min',
                'loc_phi_min', 'loc_theta_max', 'loc_phi_max']:
            setattr(hwduck, attribute, getattr(hw, attribute))

        # Extract the x, y coordinate of the stuck positioners from the
        # hardware state using their theta and phi positions.
        all_x = np.zeros(len(patch_tbl))
        all_y = np.zeros(len(patch_tbl))

        for i, loc in enumerate(patch_tbl["LOCATION"]):
            xy = hw.thetaphi_to_xy(
                                hwduck.loc_pos_curved_mm[loc],
                                hwduck.loc_theta_pos[loc] + hwduck.loc_theta_offset[loc],
                                hwduck.loc_phi_pos  [loc] + hwduck.loc_phi_offset  [loc],
                                hwduck.loc_theta_arm[loc],
                                hwduck.loc_phi_arm[loc],
                                hwduck.loc_theta_offset[loc],
                                hwduck.loc_phi_offset[loc],
                                hwduck.loc_theta_min[loc],
                                hwduck.loc_phi_min[loc],
                                hwduck.loc_theta_max[loc],
                                hwduck.loc_phi_max[loc],
                                ignore_range=True
                            )
            all_x[i] = xy[0]
            all_y[i] = xy[1]

        # Convert those xy positions to the ra and dec positions.
        ra,dec = xy2radec(hw, tile_ra, tile_dec, tile_obstime, tile_obstheta, tile_obsha,
                    all_x, all_y, False, 0)

        patch_tbl["TARGET_RA"] = ra
        patch_tbl["TARGET_DEC"] = dec

        # Patch these to match, because they also got set to nan originally
        fiber_x, fiber_y = xy2cs5(all_x, all_y)

        patch_tbl["FIBERASSIGN_X"] = fiber_x
        patch_tbl["FIBERASSIGN_Y"] = fiber_y

        # Copy the rows from the patch table into the full table.
        # Remember we onl changed the values of TARGET_RAm DEC and FIBERASSIGN X, Y
        # in the patch table, so all other values should remain unchanged.
        for i, loc in enumerate(patch_tbl["LOCATION"]):
            tbl[tbl["LOCATION"] == loc] = patch_tbl[i]

        time_string = datetime.isoformat(datetime.now(timezone.utc), timespec="seconds")
        header.add_comment(f"NaN RA/DEC Patched on {time_string}")

    # According to the docs when we exit the context manager the changes should
    # flush to the file....