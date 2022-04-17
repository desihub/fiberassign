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
from fiberassign.utils import get_mjd
from fiberassign.utils import Logger

log = Logger.get()

# AR release is hard-coded to 8888 for tertiary programs
release = 8888

# AR default values
default = {
    "TIME_MJD_BEGIN": "2020-01-01T00:00:00+00:00",
    "TIME_MJD_END": "2120-01-01T00:00:00+00:00",
    "TOO_PRIO": "HI",  # for TOO_PRIO and SCND_TARGET
}

# AR allowed values for some of the required header keywords
allow_hdrkeys = {
    "OBSCONDS": ["BRIGHT", "DARK"],  # AR not allowing BACKUP
    "SBPROF": ["ELG", "BGS", "PSF", "FLT"],
}

# AR required columns from {args.targdir}/tertiary-targets-{args.prognum}.fits
req_keys = [
    "TARGETID",  # AR should be encode_targetid(release=8888, brickid=PROGNUM, objid=np.arange(len(d)))
    "RA",
    "DEC",
    "PMRA",
    "PMDEC",
    "REF_EPOCH",
    "TERTIARY_TARGET",
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


def get_targdir(prognum):
    """
    Get the fiducial folder name for a given tertiary PROGNUM.

    Args:
        prognum: the tertiary PROGNUM (int)

    Returns:
        os.path.join(os.getenv("DESI_SURVEYOPS"), "tertiary", "{:04d}".format(prognum)) (string)
    """
    return os.path.join(os.getenv("DESI_SURVEYOPS"), "tertiary", "{:04d}".format(prognum))


def get_targfn(prognum, targdir=None):
    """
    Get the full path to a tertiary PROGNUM input targets file.

    Args:
        prognum: the tertiary PROGNUM (int)
        targdir (optional, defaults to get_targdir()): folder name (string)

    Returns:
        os.path.join(targdir, "tertiary-targets-{:04d}.fits".format(prognum)) (string)
    """
    if targdir is None:
        targdir = get_targdir(prognum)
    return os.path.join(targdir, "tertiary-targets-{:04d}.fits".format(prognum))


def get_priofn(prognum, targdir=None):
    """
    Get the full path to a tertiary PROGNUM input priorities file.

    Args:
        prognum: the tertiary PROGNUM (int)
        targdir (optional, defaults to get_targdir()): folder name (string)

    Returns:
        os.path.join(targdir, "tertiary-priorities-{:04d}.ecsv".format(prognum)) (string)
    """
    if targdir is None:
        targdir = get_targdir(prognum)
    return os.path.join(targdir, "tertiary-priorities-{:04d}.ecsv".format(prognum))


def get_toofn(prognum, tileid, targdir=None):
    """
    Get the full path to a tertiary PROGNUM ToO target file for a given tile.

    Args:
        prognum: the tertiary PROGNUM (int)
        tileid: the tile TILEID (int)
        targdir (optional, defaults to get_targdir()): folder name (string)

    Returns:
        os.path.join(targdir, "ToO-{:04d}-{:06d}.ecsv".format(prognum, tileid)) (string)
    """
    if targdir is None:
        targdir = get_targdir(prognum)
    return os.path.join(targdir, "ToO-{:04d}-{:06d}.ecsv".format(prognum, tileid))


def read_targfn(targfn):
    """
    Read a tertiary PROGNUM input targets file.

    Args:
        targfn: full path to the tertiary PROGNUM input targets file (string)

    Returns:
        targ: a Table() structured array
        targhdr: the related header
    """
    if not os.path.isfile(targfn):
        msg = "missing {}; exiting".format(targfn)
        log.error(msg)
        raise IOerror(msg)
    targ = Table.read(targfn, "TARGETS")
    targhdr = targ.meta
    log.info("reading {} targets from {}".format(len(targ), targfn))
    return targ, targhdr


def read_priofn(priofn):
    """
    Read a tertiary PROGNUM input priorities file.

    Args:
        priofn: full path to the tertiary PROGNUM input priorities file (string)

    Returns:
        prio:  Table() structured array
    """
    if not os.path.isfile(priofn):
        msg = "missing {}; exiting".format(priofn)
        log.error(msg)
        raise IOError(msg)
    prio = Table.read(priofn)
    return prio


def assert_tertiary_targ(prognum, targfn):
    """
    Verifies that a tertiary PROGNUM input targets file is correctly formatted.

    Args:
        prognum: the tertiary PROGNUM (int)
        targfn: full path to the tertiary PROGNUM input targets file (string)

    Notes:
        Will raise an error if the file is not correctly formatted.
    """
    # AR read
    targ, targhdr = read_targfn(targfn)

    # AR check all required keys are there
    miss_keys = []
    for key in req_keys:
        if key not in targ.dtype.names:
            miss_keys.append(key)
    if len(miss_keys) > 0:
        msg = "{} keys are missing from {}; exiting".format(",".join(miss_keys), targfn)
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
                hdrkey,
                targhdr[hdrkey],
                allow_hdrkeys[hdrkey],
            )
            log.error(msg)
            raise IOError(msg)

    # AR check TARGETID
    targetids = encode_targetid(
        release=8888, brickid=prognum, objid=np.arange(len(targ))
    )
    sel = targ["TARGETID"] != targetids
    if sel.sum() > 0:
        msg = "TARGETID is not as expected, i.e. encode_targetid(release=8888, brickid={}, objid=np.arange(len(targ)))".format(
            prognum,
        )
        log.error(msg)
        raise IOError(msg)

    # AR check RA, DEC
    sel = (targ["RA"] >= 0) & (targ["RA"] < 360)
    sel &= (targ["DEC"] >= -90) & (targ["DEC"] <= 90)
    if sel.sum() != len(targ):
        msg = "{} targets do not verify 0 <= RA < 360 and -90 <= DEC <= 90".format(
            len(targ) - sel.sum()
        )
        log.error(msg)
        raise IOError(msg)

    # AR check PMRA, PMDEC, REF_EPOCH
    sel = (np.isfinite(targ["PMRA"])) & (np.isfinite(targ["PMDEC"]))
    if sel.sum() != len(targ):
        msg = "{} targets do have not finite PMRA or PMDEC".format(
            len(targ) - sel.sum()
        )
        log.error(msg)
        raise IOError(msg)
    sel = (np.isfinite(targ["REF_EPOCH"])) & (targ["REF_EPOCH"] > 0)
    if sel.sum() != len(targ):
        msg = "{} targets do have not finite REF_EPOCH or REF_EPOCH <= 0".format(
            len(targ) - sel.sum()
        )
        log.error(msg)
        raise IOError(msg)

    # AR check SUBPRIORITY
    # AR - 0 < SUBPRIORITY < 1 (not sure 0 and 1 have to be excluded.. in doubt do so)
    # AR - check uniqueness of values, which is expected for a random (only raises a warning)
    if "SUBPRIORITY" in targ.dtype.names:
        sel = (targ["SUBPRIORITY"] > 0) & (targ["SUBPRIORITY"] < 1)
        if sel.sum() != len(targ):
            msg = "{} targets do not have 0 < SUBPRIORITY < 1".format(
                len(targ) - sel.sum()
            )
            log.error(msg)
            raise IOError(msg)
        tmpn = np.unique(targ["SUBPRIORITY"]).size
        if tmpn != len(targ):
            msg = "SUBPRIORITY does not look like a random draw ({} unique values for {} rows)".format(
                tmpn, len(targ)
            )
            log.warning(msg)


    # AR warning if some columns are present but will be overwritten
    keys = np.array(targ.dtype.names)
    sel = np.in1d(
        keys,
        [
            "SCND_ORDER",
            "PRIORITY_INIT",
            "NUMOBS_INIT",
            "NUMOBS",
            "NUMOBS_MORE",
            "PRIORITY",
        ],
    )
    if sel.sum() > 0:
        log.warning(
            "The following columns present in the input file will be overwritten: {}".format(
                keys[sel]
            )
        )


def assert_tertiary_prio(prognum, priofn, targfn):
    """
    Verifies that a tertiary PROGNUM input prioritiess file is correctly formatted.

    Args:
        prognum: the tertiary PROGNUM (int)
        priofn: full path to the tertiary PROGNUM input priorities file (string)
        targfn: full path to the tertiary PROGNUM input targets file (string)

    Notes:
        Will raise an error if the file is not correctly formatted.
    """
    # AR read the files
    prio = read_priofn(priofn)
    targ, _ = read_targfn(targfn)

    # AR check that all TERTIARY_TARGET values in targ exist in prio
    targ_tertiary_targets = np.unique(targ["TERTIARY_TARGET"])
    prio_tertiary_targets = np.unique(prio["TERTIARY_TARGET"])
    sel = np.array(
        [
            tertiary_target not in prio_tertiary_targets
            for tertiary_target in np.unique(targ["TERTIARY_TARGET"])
        ]
    )
    if sel.sum() > 0:
        msg = "TERTIARY_TARGET={} not present in prio".format(
            ",".join(targ_tertiary_targets[sel])
        )
        log.error(msg)
        raise IOError(msg)

    # AR check:
    # AR - NUMOBS_DONE_MIN=0, NUMOBS_DONE_MAX=99
    # AR - NUMOBS_DONE_MIN, NUMOBS_DONE_MAX consistency
    for tertiary_target in np.unique(prio["TERTIARY_TARGET"]):
        # AR check NUMOBS_DONE_MIN=0, NUMOBS_DONE_MAX=99
        for numobs_key, numobs_val in zip(
            ["NUMOBS_DONE_MIN", "NUMOBS_DONE_MAX"], [0, 99]
        ):
            sel = (prio["TERTIARY_TARGET"] == tertiary_target) & (
                prio["NUMOBS_DONE_MIN"] == 0
            )
            if sel.sum() == 0:
                msg = "NUMOBS_DONE_MIN=0 case not set for TERTIARY_TARGET={}".format(
                    tertiary_target
                )
                log.error(msg)
                raise IOError(msg)
            if sel.sum() > 1:
                msg = "NUMOBS_DONE_MIN=0 case appearing {} times for TERTIARY_TARGET={}".format(
                    sel.sum(), tertiary_target
                )
                log.error(msg)
                raise IOError(msg)
        # AR NUMOBS_DONE_MIN, NUMOBS_DONE_MAX consistency
        sel = prio["TERTIARY_TARGET"] == tertiary_target
        ok = True
        if np.unique(prio["NUMOBS_DONE_MIN"][sel]).size != sel.sum():
            ok = False
        if np.unique(prio["NUMOBS_DONE_MAX"][sel]).size != sel.sum():
            ok = False
        for numobs_done_min, numobs_done_max in zip(
            prio["NUMOBS_DONE_MIN"][sel], prio["NUMOBS_DONE_MAX"][sel]
        ):
            if numobs_done_min > numobs_done_max:
                ok = False
        if not ok:
            msg = "NUMOBS_DONE_MIN, NUMOBS_DONE_MAX ill-designed for TERTIARY_TARGET={}".format(
                tertiary_target
            )
            log.error(msg)
            raise IOError(msg)


def get_numobs_priority(too, prio, prognum, previous_tileids=None, fadir=None):
    """
    Obtain the NUMOBS, NUMOBS_MORE and PRIORITY columns of a tertiary ToO catalog.

    Args:
        too: Table() structured array
        prio: Table() structured array (output from read_priofn())
        prognum: the tertiary PROGNUM (int)
        previous_tileids (optional, defaults to None): comma-separated list of already designed TILEIDs to consider to set NUMOBS and NUMOBS_MORE (string)
        fadir (optional, defaults to None): folder with the fiberassign files for the previous_tileids (string)

    Returns:
        numobss: NUMOBS values (np.array())
        numobs_mores: NUMOBS_MORE values (np.array())
        priorities: PRIORITY values (np.array())

    Notes:
        If fadir is provided the code will first look in fadir, then will look in $DESI_TARGET/fiberassign/tiles/trunk.
        If fadir is not provided, the code will look in $DESI_TARGET/fiberassign/tiles/trunk.
    """
    # AR initialize NUMOBS, NUMOBS_MORE, PRIORITY
    numobss = np.zeros(len(too), dtype=int)
    numobs_mores = too["NUMOBS_INIT"].copy()
    priorities = too["PRIORITY_INIT"].copy()

    # AR previous assignment?
    # AR we proceed as follows:
    # AR - each time a target is assigned, we do NUMOBS_MORE-=1 and NUMOBS+=1
    # AR - when NUMOBS_MORE=0, we change PRIORITY to PRIORITY_DONE
    if previous_tileids is None:
        log.info(
            "previous_tileids is None -> setting NUMOBS=0, NUMOBS_MORE=NUMOBS_INIT, PRIORITY=PRIORITY_INIT"
        )
    else:
        for tileid in previous_tileids.split(","):
            tileid = int(tileid)
            isfn = False
            # AR first try fadir, if provided
            # AR (case where tiles are designed on the same day)
            if fadir is not None:
                fn = os.path.join(
                    fadir,
                    "{}".format("{:06d}".format(tileid)[:3]),
                    "fiberassign-{:06d}.fits.gz".format(tileid),
                )
                if os.path.isfile(fn):
                    isfn = True
                else:
                    log.inf("no {}".format(fn))
            # AR then try the DESI_TARGET/fiberassign if no file in fadir (or fadir not provided)
            if not isfn:
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
                _, prev_prognum, prev_release, _, _, _ = decode_targetid(
                    prev_d["TARGETID"]
                )
                sel = np.in1d(prev_release, release) & np.in1d(prev_prognum, prognum)
                prev_d = prev_d[sel]
                sel = np.in1d(too["TARGETID"], prev_d["TARGETID"])
                # AR update columns
                numobss[sel] += 1
                numobs_mores[sel] -= 1
                numobs_mores = np.clip(numobs_mores, 0, None)
                log.info(
                    "{} targets were already assigned on TILEID={}".format(
                        sel.sum(),
                        tileid,
                    )
                )
            else:
                log.warning(
                    "no fiberassign-{:06d}.fits.gz file found: skipping this tile".format(
                        tileid
                    )
                )

        # AR now update PRIORITY
        for i in range(len(prio)):
            sel = (
                ((too["TERTIARY_TARGET"] == prio["TERTIARY_TARGET"][i]) > 0)
                & (numobss >= prio["NUMOBS_DONE_MIN"][i])
                & (numobss <= prio["NUMOBS_DONE_MAX"][i])
            )
            priorities[sel] = prio["PRIORITY"][i]
            log.info(
                "setting PRIORITY={} for {} targets with TERTIARY_TARGET = {} and {} <= NUMOBS <= {}".format(
                    prio["PRIORITY"][i],
                    sel.sum(),
                    prio["TERTIARY_TARGET"][i],
                    prio["NUMOBS_DONE_MIN"][i],
                    prio["NUMOBS_DONE_MAX"][i],
                )
            )
    return numobss, numobs_mores, priorities


def create_tertiary_too(args):
    """
    Create the ToO-args.prognum-args.tileid.ecsv file.

    Args:
        args: the argument output by get_parse() in fba_tertiary_too
    """
    # AR output file
    toofn = get_toofn(args.prognum, args.tileid, targdir=args.targdir)
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

    # AR targets and priorities files
    targfn = get_targfn(args.prognum, targdir=args.targdir)
    priofn = get_priofn(args.prognum, targdir=args.targdir)

    # AR check input target catalog and priority file
    # AR (not optimal in term of processing time,
    # AR    as we read targfn three times here,
    # AR    in assert_tertiary_targ(),
    # AR    in assert_tertiary_prio(),
    # AR    and in the code after;
    # AR    but that insures that the checks are done
    # AR    using the exact same reading function
    # AR    than what is done in the code;
    # AR    but processing time here is not a limiting factor
    # AR )
    assert_tertiary_targ(args.prognum, targfn)
    assert_tertiary_prio(args.prognum, priofn, targfn)

    # AR read targets and priorities files
    targ, targhdr = read_targfn(targfn)
    prio = read_priofn(priofn)

    # AR scnd_mask_name
    if args.scnd_mask_name is None:
        args.scnd_mask_name = "{}_TOO_{}P".format(
            targhdr["OBSCONDS"], default["TOO_PRIO"]
        )
        log.info(
            "setting args.scnd_mask_name = '{}_TOO_{}P'".format(
                targhdr["OBSCONDS"], default["TOO_PRIO"]
            )
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
            sel.sum(),
            len(targ),
            args.tilera,
            args.tiledec,
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
    for key in ["TARGETID", "RA", "DEC", "PMRA", "PMDEC", "REF_EPOCH"] + ["SCND_ORDER"]:
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
    # AR temporary adding TERTIARY_TARGET
    too["TERTIARY_TARGET"] = targ["TERTIARY_TARGET"]

    # AR SUBPRIORITY
    # AR creating it if not present
    if "SUBPRIORITY" in targ.dtype.names:
        log.info("SUBPRIORITY column present in {} -> propagating it".format(targfn))
        too["SUBPRIORITY"] = targ["SUBPRIORITY"]
    else:
        log.info("SUBPRIORITY column not present in {} -> creating it".format(targfn))
        too["SUBPRIORITY"] = np.random.uniform(size=ntarg)

    # AR NUMOBS_INIT, PRIORITY_INIT
    too["NUMOBS_INIT"], too["PRIORITY_INIT"] = 0, 0
    for tertiary_target in np.unique(prio["TERTIARY_TARGET"]):
        # AR PRIORITY_INIT
        sel = (prio["TERTIARY_TARGET"] == tertiary_target) & (
            prio["NUMOBS_DONE_MIN"] == 0
        )
        priority_init = prio["PRIORITY"][sel][0]
        # AR NUMOBS_INIT
        sel = (prio["TERTIARY_TARGET"] == tertiary_target) & (
            prio["NUMOBS_DONE_MAX"] == 99
        )
        numobs_init = prio["NUMOBS_DONE_MIN"][sel][0]
        #
        sel = too["TERTIARY_TARGET"] == tertiary_target
        too["PRIORITY_INIT"][sel] = priority_init
        too["NUMOBS_INIT"][sel] = numobs_init
        log.info(
            "setting PRIORITY_INIT={} and NUMOBS_INIT={} for {} targets".format(
                priority_init, numobs_init, sel.sum()
            )
        )

    # AR NUMOBS, NUMOBS_MORE, PRIORITY
    numobss, numobs_mores, priorities = get_numobs_priority(
        too,
        prio,
        args.prognum,
        previous_tileids=args.previous_tileids,
        fadir=args.fadir,
    )
    too["NUMOBS"] = numobss
    too["NUMOBS_MORE"] = numobs_mores
    too["PRIORITY"] = priorities

    # AR remove TERTIARY_TARGET
    too.remove_column("TERTIARY_TARGET")

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
    too.meta["PRIOFN"] = priofn
    for hdrkey in req_hdrkeys:
        too.meta[hdrkey] = targhdr[hdrkey]

    # AR write
    too.write(toofn)
