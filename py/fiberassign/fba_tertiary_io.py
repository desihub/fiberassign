#!/usr/bin/env python

import os
import numpy as np
from astropy.table import Table
from time import time
from datetime import datetime, timezone
from astropy.time import Time
from desitarget.targetmask import desi_mask
from desitarget.targetmask import obsconditions
from desitarget.targets import encode_targetid, decode_targetid
from desitarget.targetmask import desi_mask, scnd_mask
from desitarget.mtl import get_utc_date
from desimodel.footprint import is_point_in_desi
from fiberassign.utils import assert_isoformat_utc
from fiberassign.utils import Logger

log = Logger.get()

# AR release is hard-coded to 8888 for tertiary programs
release = 8888

# AR default values
default = {
    "TIME_MJD_BEGIN": "2020-01-01T00:00:00+00:00",
    "TIME_MJD_END": "2120-01-01T00:00:00+00:00",
    "TOO_PRIO": "HI", # for TOO_PRIO and SCND_TARGET
}

# AR allowed values for some of the required header keywords
allow_hdrkeys = {
    "OBSCONDS" : ["BRIGHT", "DARK"], # AR not allowing BACKUP
    "SBPROF" : ["ELG", "BGS", "PSF", "FLT"],
}

# AR required columns from {args.targdir}/tertiary-targets-{args.prognum}.fits
req_keys = [
    "RA",
    "DEC",
    "PMRA",
    "PMDEC",
    "REF_EPOCH",
    "NUMOBS_INIT",
    "PRIORITY_INIT",
    "PRIORITY_DONE",
    "CHECKER",
]

# AR required keywords in the header
req_hdrkeys = [
    "FAPRGRM",
    "OBSCONDS",
    "SBPROF",
    "GOALTIME",
]

# AR reproducibility for SUBPRIORITY
np_rand_seed = 1234
np.random.seed(np_rand_seed)


def get_targfn(targdir, prognum):
    return os.path.join(
        targdir, "tertiary-targets-{:04d}.fits".format(prognum)
    )

def get_toofn(targdir, prognum, tileid):
    return os.path.join(
        targdir, "ToO-{:04d}-{:06d}.ecsv".format(prognum, tileid)
    )


def get_mjd(utc_time_mjd):
    assert_isoformat_utc(utc_time_mjd)
    return Time(datetime.strptime(utc_time_mjd, "%Y-%m-%dT%H:%M:%S%z")).mjd


def assert_tertiary_targ(targ, targhdr):

    # AR check all required keys are there
    miss_keys = []
    for key in req_keys:
        if key not in targ.dtype.names:
            miss_keys.append(key)
    if len(miss_keys) > 0:
        msg = "{} keys are missing from {}; exiting".format(
            ",".join(miss_keys), targfn
        )
        log.error(msg)
        raise IOError(msg)

    # AR check all required header keywords are there
    miss_hdrkeys = []
    for hdrkey in req_hdrkeys:
        if hdrkey not in targhdr:
            miss_hdrkeys.append(hdrkey)
    if len(miss_hdrkeys) > 0:
        msg = "{} keys are missing from the header of {}; exiting".format(
            ",".join(miss_hdrkeys), targfn
        )
        log.error(msg)
        raise IOError(msg)

    # AR check that the required header keywords are as allowed
    for hdrkey in allow_hdrkeys:
        if targhdr[hdrkey] not in allow_hdrkeys[hdrkey]:
            msg = "hdr['{}']={} not in allowed values ({})".format(
                hdrkey, targhdr[hdrkey], allow_hdrkeys[hdrkey],
            )
            log.error(msg)
            raise IOError(msg)


def get_numobs_priority(too, priority_dones, prognum, previous_tileids=None, fadir=None):

    # AR initialize NUMOBS, NUMOBS_MORE, PRIORITY
    numobss = np.zeros(len(too), dtype=int)
    numobs_mores = too["NUMOBS_INIT"].copy()
    priorities = too["PRIORITY_INIT"].copy()

    # AR previous assignment?
    # AR we proceed as follows:
    # AR - each time a target is assigned, we do NUMOBS_MORE-=1 and NUMOBS+=1
    # AR - when NUMOBS_MORE=0, we change PRIORITY to PRIORITY_DONE
    if previous_tileids is None:
        log.info("previous_tileids is None -> setting NUMOBS=0, NUMOBS_MORE=NUMOBS_INIT, PRIORITY=PRIORITY_INIT")
    else:
        for tileid in previous_tileids.split(","):
            tileid = int(tileid)
            # AR first try args.fadir (case where tiles are designed on the same day)
            fn = os.path.join(
                fadir,
                "{}".format("{:06d}".format(tileid)[:3]),
                "fiberassign-{:06d}.fits.gz".format(tileid),
            )
            isfn = False
            if os.path.isfile(fn):
                isfn = True
            # AR then try the DESI_TARGET/fiberassign
            else:
                log.inf("no {}".format(fn))
                fn = os.path.join(
                    os.getenv("DESI_TARGET"),
                    "fiberassign",
                    "tiles",
                    "trunk",
                    "{}".format("{:06d}".format(tileid)[:3]),
                    "fiberassign-{:06d}.fits.gz".format(tileid),
                )
                if os.path.isfile(fn):
                    isfn = True
                else:
                    log.info("no {}".format(fn))
            if isfn:
                log.info("reading {}".format(fn))
                # AR identify assigned targets
                prev_d = Table.read(fn, "FIBERASSIGN")
                _, prev_prognum, prev_release, _, _, _ = decode_targetid(prev_d["TARGETID"])
                sel = np.in1d(prev_release, release) & np.in1d(prev_prognum, prognum)
                prev_d = prev_d[sel]
                sel = np.in1d(too["TARGETID"], prev_d["TARGETID"])
                # AR update columns
                numobss[sel] += 1
                numobs_mores[sel] -= 1
                numobs_mores = np.clip(numobs_mores, 0, None)
                log.info(
                    "{} targets were already assigned on TILEID={}".format(
                        sel.sum(), tileid,
                    )
                )
            else:
                log.warning(
                    "no fiberassign-{:06d}.fits.gz file found: skipping this tile".format(
                        tileid
                    )
                )

        # AR now update PRIORITY for NUMOBS_MORE = 0 targets
        sel = numobs_mores == 0
        priorities[sel] = priority_dones[sel]
        log.info("PRIORITY updated to PRIORITY_DONE for {} targets".format(sel.sum()))

    return numobss, numobs_mores, priorities


def create_tertiary_too(args):

    # AR output file
    toofn = get_toofn(args.targdir, args.prognum, args.tileid)
    if os.path.isfile(toofn):
        msg = "{} already exists! exiting".format(toofn)
        log.error(msg)
        raise IOError(msg)
    else:
        log.info("toofn = {}".format(toofn))

    # AR mjd_begin, mjd_end
    mjd_begin = get_mjd(args.utc_time_mjd_begin)
    mjd_end = get_mjd(args.utc_time_mjd_end)
    log.info("set MJD_BEGIN={}, MJD_END={}".format(mjd_begin, mjd_end))

    # AR read input file
    targfn = get_targfn(args.targdir, args.prognum)
    if not os.path.isfile(targfn):
        msg = "missing {}; exiting".format(targfn)
        log.error(msg)
        raise IOerror(msg)
    targ = Table.read(targfn)
    targhdr = targ.meta
    log.info("reading {} targets from {}".format(len(targ), targfn))

    # AR check input catalog
    assert_tertiary_targ(targ, targhdr)

    # AR scnd_mask_name
    if args.scnd_mask_name is None:
        args.scnd_mask_name = "{}_TOO_{}P".format(targhdr["OBSCONDS"], default["TOO_PRIO"])
        log.info("setting args.scnd_mask_name = '{}_TOO_{}P'".format(targhdr["OBSCONDS"], default["TOO_PRIO"]))

    # AR TARGETID (possibly overwrite existing one)
    # AR we use args.prognum as BRICKID
    # AR need to be done *before* cutting on the tile footprint
    if "TARGETID" in targ.dtype.names:
        log.warning("overwriting the original TARGETIDs!")
    targ["TARGETID"] = encode_targetid(
        release=release, brickid=args.prognum, objid=np.arange(len(targ))
    )

    # AR SCND_ORDER
    # AR record the row from the tertiary-targets-{args.prognum}.fits
    if "SCND_ORDER" in targ.dtype.names:
        log.warning("overwriting the original SCND_ORDERs!")
    targ["SCND_ORDER"] = np.arange(len(targ), dtype=int)


    # AR restrict to targets in the tile
    tiles = Table()
    tiles["RA"], tiles["DEC"] = [args.tilera], [args.tiledec]
    sel = is_point_in_desi(tiles, targ["RA"], targ["DEC"])
    log.info(
        "keep {}/{} targets in tile with (tilera,tiledec)={},{}".format(
            sel.sum(), len(targ), args.tilera, args.tiledec,
        )
    )
    targ = targ[sel]
    ntarg = len(targ)

    # AR first get the correct structure
    reffn = os.path.join(os.getenv("DESI_SURVEYOPS"), "mtl", "main", "ToO", "ToO.ecsv")
    dref = Table.read(reffn)

    # AR build output file
    too = Table()
    for key in dref.dtype.names:
        too[key] = np.zeros_like(dref[key], shape=(ntarg))
    for key in req_keys + ["TARGETID", "SCND_ORDER"]:
        too[key] = targ[key]
    too["DESI_TARGET"] = desi_mask[args.desi_mask_name]
    too["SCND_TARGET"] = scnd_mask[args.scnd_mask_name]
    too["OBSCONDITIONS"] = obsconditions[targhdr["OBSCONDS"]]
    too["TOO_TYPE"] = "TILE"
    too["TOO_PRIO"] = default["TOO_PRIO"]
    too["OCLAYER"] = targhdr["OBSCONDS"]
    too["MJD_BEGIN"] = mjd_begin
    too["MJD_END"] = mjd_end
    too["TIMESTAMP"] = get_utc_date("main")


    # AR SUBPRIORITY
    # AR creating it if not present
    if "SUBPRIORITY" in targ.dtype.names:
        log.info("SUBPRIORITY column present in {} -> propagating it".format(targfn))
        too["SUBPRIORITY"] = targ["SUBPRIORITY"]
    else:
        log.info("SUBPRIORITY column not present in {} -> creating it".format(targfn))
        too["SUBPRIORITY"] = np.random.uniform(size=ntarg)


    # AR NUMOBS, NUMOBS_MORE, PRIORITY
    numobss, numobs_mores, priorities = get_numobs_priority(
        too,
        targ["PRIORITY_DONE"],
        args.prognum,
        previous_tileids=args.previous_tileids,
        fadir=args.fadir,
    )
    too["NUMOBS"] = numobss
    too["NUMOBS_MORE"] = numobs_mores
    too["PRIORITY"] = priorities

    # AR store args (we exclude any None argument)
    tmparr = []
    for kwargs in args._get_kwargs():
        if kwargs[1] is not None:
            tmparr += ["--{} {}".format(kwargs[0], kwargs[1])]
    tooargs = " ".join(tmparr)
    log.info("TOOARGS = {}".format(tooargs))
    too.meta["TOOARGS"] = tooargs
    # AR store other infos
    too.meta["TARGFN"] = targfn
    for hdrkey in req_hdrkeys:
        too.meta[hdrkey] = targhdr[hdrkey]

    # AR write
    too.write(toofn)


if __name__ == "__main__":
    main()
