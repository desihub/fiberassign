#!/usr/bin/env python

import os
import subprocess
from time import time, sleep
import numpy as np
from astropy.time import Time
from astropy.table import Table, vstack
from fiberassign.utils import get_revs_dates, get_rev_fiberassign_changes
from desiutil.log import get_logger
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

log = get_logger()


def parse():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description="Get the list of svn revisions for the fiberassign files for a set of tiles",
        epilog="*** No parallelization as this can bring down the desi server ! ***\n",
    )
    parser.add_argument("--outfn", help="output file", type=str, default=None)
    parser.add_argument(
        "--tilesfn",
        help="full path to prod tiles file, like $DESI_ROOT/spectro/redux/iron/tiles-iron.csv",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--svndir",
        help="svn folder",
        type=str,
        default=os.path.join(os.getenv("DESI_TARGET"), "fiberassign", "tiles", "trunk"),
    )
    args = parser.parse_args()
    for kwargs in args._get_kwargs():
        log.info(kwargs)
    return args


def main():

    # AR
    all_start = time()
    log.info("start\tTIMESTAMP={}".format(Time.now().isot))
    args = parse()

    # AR read the tiles file
    tiles = Table.read(args.tilesfn)
    log.info("looking at {} tiles".format(len(tiles)))

    # AR first identify the relevant subdirs folders (for speed-up)
    subdirs = np.unique(["{:06d}".format(tileid)[:3] for tileid in tiles["TILEID"]])
    subdirs = ",".join(subdirs)
    log.info(
        "looking at the following {} subdirs: {}".format(
            len(subdirs.split(",")), subdirs
        )
    )

    # AR create the unique list of all revisions
    revs, _ = get_revs_dates(args.svndir, subdirs)
    log.info("found {} revisions".format(revs.size))

    # AR get all changes to fiberassign files
    # AR (there could be duplicated files if they appear in different revisions)
    ds = []
    for rev in revs:
        d = get_rev_fiberassign_changes(args.svndir, rev, subdirs=subdirs)
        ds.append(d)
        # AR protect from overloading the server
        sleep(0.5)

    # AR stack
    d = vstack(ds)
    d.meta["TILESFN"] = args.tilesfn

    # AR cut on tiles
    sel = np.in1d(d["TILEID"], tiles["TILEID"])
    d = d[sel]
    # AR write
    if os.path.splitext(args.outfn)[1] == ".asc":
        d.write(args.outfn, format="ascii.commented_header")
    else:
        d.write(args.outfn)

    # AR time report
    all_dt = time() - all_start
    log.info("done\tTIMESTAMP={}\t(took {:.1f}s)".format(Time.now().isot, all_dt))


if __name__ == "__main__":
    main()
