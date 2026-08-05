#!/usr/bin/env python

# Script that patches mistaken TARGETIDs in calibration tile 83063.
# This is necessary due to an issue where the input targets for that tile
# was accidentally regenerated, see https://github.com/desihub/desitarget/issues/892
# Example call:
# python patch_83063.py /path/to/fiberassign-083063.fits.gz $DESI_SURVEYOPS/tertiary/
# This will patch whatever `fiberassign-083063.fits.gz` file you have with the tertiary
# program details in $DESI_SURVEYOPS/tertiary/.


import numpy as np
from astropy.table import Table, join
from astropy.io import fits

from desitarget.targets import decode_targetid

import argparse
from datetime import datetime, timezone
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("tile_file", type=str, help="Broken tile file to patch.")
parser.add_argument("tertiary", type=str, help="where the tertiary data is stored")
args = parser.parse_args()

tile_loc = Path(args.tile_file)
tertiary_dir = Path(args.tertiary)
# fa_saved = tile_dir / "fiberassign-083063.fits.gz"

old_targ = Table.read(tertiary_dir / "0008" / "tertiary-targets-0008.fits")
wrong_targ = Table.read(tertiary_dir / "0008" / "tertiary-targets-0008-83063.fits")

# Sanity check
assert np.all(np.isin(wrong_targ["ORIG_TARGETID"], old_targ["ORIG_TARGETID"])), "Some wrong targetids don't have a match in the correct file!"

cols = ["TARGETID", "ORIG_TARGETID"]
mapping = join(old_targ[cols], wrong_targ[cols], keys="ORIG_TARGETID", table_names=["CORRECT", "WRONG"])

with fits.open(tile_loc, mode="update") as h:
    # header = h[0].header
    time_string = datetime.isoformat(datetime.now(timezone.utc), timespec="seconds")
    for i, hdu in enumerate(h):
        print(f"Patching {hdu.name}...")
        tbl = hdu.data
        header = hdu.header

        if i == 0:
            if "TARGETIDs Patched" in header["COMMENT"][-1]:
                print("This file has already been patched! Aborting!")
                break
            header.add_comment(f"TARGETIDs Patched on {time_string}")

        if tbl is not None:
            for i, row in enumerate(tbl):
                tid = row["TARGETID"]
                if tid < 0: continue # Skip the stucks.

                rs = decode_targetid(tid)[2]
                if rs == 8888: # These are the ones that got new targetids in the fba_calibration pipeline.
                    idx = np.where(mapping["TARGETID_WRONG"] == tid)[0][0]
                    row["TARGETID"] = mapping[idx]["TARGETID_CORRECT"]

