#!/usr/bin/env python

import sys
import os
from glob import glob
import tempfile
import yaml
from pkg_resources import resource_filename
from datetime import datetime
import numpy as np
import fitsio
import healpy as hp
from astropy.io import fits
from astropy.table import Table, vstack, hstack
from astropy import units
from astropy.coordinates import SkyCoord
from desimodel.footprint import is_point_in_desi, tiles2pix
from desitarget.geomask import match, pixarea2nside
from desitarget.io import read_targets_in_tiles
from desitarget.targetmask import obsconditions
from desitarget.gaiamatch import gaia_psflike
from fiberassign.assign import merged_fiberassign_swap
from desitarget.targets import main_cmx_or_sv
from fiberassign.utils import Logger
from desiutil.dust import SFDMap
from desiutil.redirect import stdouterr_redirected


log = Logger.get()


# AR fiberassign extensions
# AR POTENTIAL_ASSIGNMENTS only have TARGETID,FIBER,LOCATION, so nothing to patch
all_exts = ["FIBERASSIGN", "SKY_MONITOR", "GFA_TARGETS", "TARGETS"]

# AR possible sources of desitarget catalogs
# AR all_names refer to the header keywords in the 0-extension of the fiberassign-TILEID.fits.gz files
all_names = ["targ", "sky", "gfa", "scnd", "too"]


def get_tileid_from_fafn(fafn):
    """
    Returns the tileid for a single fiberassign file name.

    Args:
        fafn: fiberassign file name (fiberassign-ABCDE.fits.gz) (str)

    Returns:
        tileid: tileid (int)
    """
    return int(os.path.basename(fafn)[12:18])


def get_fafns_to_check(fa_srcdir, skip_subdirs=None, only_subdirs=None, only_tileids=None, debug=False):
    """
    List of fiberassign files to check for patching diagnosis.

    Args:
        fa_srcdir: source folder for input fiberassign files (str)
        skip_subdirs (optional, defaults to None): comma-separated list of subdirs (e.g., '000,001') to skip (str)
        only_subdirs (optional, defaults to None): comma-separated list of subdirs (e.g., '000,001') to restrict to (str)
        only_tileids (optional, defaults to None): comma-separated list of tileids (e.g. '82405,82406') to restrict to (str)
        debug (optional, defaults to False): if set, picks one fiberassign file per FAFLAVOR (bool)

    Returns:
        fafns: list of fiberassign files to check (np.array())
        tileids: corresponding tileids (np.array())
        subdirs: corresponding subdirs (np.array())

    Notes:
        The function will look for {fa_srcdir}/???/fiberassign-??????.fits.gz files.
    """
    log.info("fa_srcdir = {}".format(fa_srcdir))
    fafns = np.sort(
        glob(
            os.path.join(
                fa_srcdir,
                "???",
                "fiberassign-??????.fits*",
            )
        )
    )
    tileids = np.array([get_tileid_from_fafn(fafn) for fafn in fafns])
    subdirs = np.array(["{:06d}".format(tileid)[:3] for tileid in tileids])
    log.info("start with {} fiberassign-TILEID.fits* files".format(fafns.size))

    # AR 50000 <= TILEID < 80000
    reject = (tileids >= 50000) & (tileids < 80000)
    log.info("reject {} tiles with 50000 <= TILEID < 80000".format(reject.sum()))
    fafns, tileids, subdirs = fafns[~reject], tileids[~reject], subdirs[~reject]

    # AR 82248 <= TILEID < 82258
    reject = (tileids >= 82248) & (tileids < 82258)
    log.info(
        "reject {} tiles with 82248 <= TILEID < 82258 (cmxposmapping)".format(
            reject.sum()
        )
    )
    fafns, tileids, subdirs = fafns[~reject], tileids[~reject], subdirs[~reject]
    # AR
    log.info("keep {} tiles".format(fafns.size))

    # AR skip subdirs?
    if skip_subdirs is not None:
        reject = np.zeros(len(fafns), dtype=bool)
        for subdir in skip_subdirs.split(","):
            reject_i = subdirs == subdir
            reject |= reject_i
            log.info(
                "{}\tremove {} fiberassign files, as per skip_subdirs".format(
                    subdir, reject_i.sum()
                )
            )
        fafns, tileids, subdirs = fafns[~reject], tileids[~reject], subdirs[~reject]

    # AR only subdirs?
    if only_subdirs is not None:
        sel = np.zeros(len(fafns), dtype=bool)
        for subdir in only_subdirs.split(","):
            sel_i = subdirs == subdir
            sel |= sel_i
            log.info(
                "{}\trestrict to {} fiberassign files, as per only_subdirs".format(
                    subdir, sel_i.sum()
                )
            )
        fafns, tileids, subdirs = fafns[sel], tileids[sel], subdirs[sel]

    # AR only tileids?
    if only_tileids is not None:
        sel = np.in1d(tileids, [int(tileid) for tileid in only_tileids.split(",")])
        log.info(
            "{}\trestrict to {} fiberassign files, as per only_tileids".format(
                only_tileids, sel.sum()
            )
        )
        fafns, tileids, subdirs = fafns[sel], tileids[sel], subdirs[sel]

    # AR debug: pick one tileid per faflavor
    if debug:
        d = Table.read(
            os.path.join(
                os.getenv("DESI_ROOT"), "spectro", "redux", "daily", "tiles-daily.csv"
            )
        )
        d = d[~np.in1d(d["FAFLAVOR"], ["cmxposmapping", "unknown"])]
        faflavors, ii = np.unique(d["FAFLAVOR"], return_index=True)
        sel = np.in1d(tileids, d["TILEID"][ii])
        if sel.sum() != ii.size:
            log.warning(
                "debug: only {} / {} FAFLAVORs picked".format(sel.sum(), ii.size),
            )
        fafns, tileids, subdirs = fafns[sel], tileids[sel], subdirs[sel]
        log.warning("debug: keep {} tiles".format(fafns.size))

    # AR what is left
    for subdir in np.unique(subdirs):
        log.info(
            "{}\tworking with {} fiberassign files".format(
                subdir, (subdirs == subdir).sum()
            )
        )

    return fafns, tileids, subdirs


# AR name: targ, sky, gfa, scnd, too
def get_static_desitarget_fns(fafn, name):
    """
    Obtain the full paths of the input static (=not MTL except for ToOs) desitarget catalogs
        which have been used to design the tile.

    Args:
        fafn: full path to the fiberassign file name (fiberassign-ABCDE.fits.gz) (str)
        name: "targ", "scnd", "too", "sky", or "gfa" (str)

    Returns:
        fns: list of static desitarget catalogs (list of strs)

    Notes:
        The function fixes all the "historical" inconsistencies in the paths.
        Note that the return "catalogs" paths can be folders.
    """
    tileid = get_tileid_from_fafn(fafn)

    # AR valid name?
    if name not in all_names:
        msg = "{:06d}\tNAME={} not allowed".format(tileid, name)
        log.error(msg)
        raise ValueError(msg)

    # AR get relevant paths
    fahdr = fits.getheader(fafn, 0)
    hdrkeys = [name.upper()]

    # AR check for possible additional secondary folders
    # AR and also for targ, as it happens for 80901 and 80865
    # AR    (but that could be the only cases)
    for i in range(2, 100):
        if (name == "targ") & ("TARG{}".format(i) in fahdr):
            hdrkeys.append("TARG{}".format(i))
        if (name == "scnd") & ("SCND{}".format(i) in fahdr):
            hdrkeys.append("SCND{}".format(i))
    if (name == "sky") & ("SKYSUPP" in fahdr):
        hdrkeys.append("SKYSUPP")
    log.info(
        "{:06d}\t{}\tcheck {} hdrkeys".format(
            tileid, name.upper().ljust(11), ",".join(hdrkeys)
        )
    )

    # AR fns can be empty (e.g. TOO for early tiles)
    fns = [fahdr[hdrkey] for hdrkey in hdrkeys if hdrkey in fahdr]
    # AR special case for 80615...
    if (tileid == 80615) & (name == "targ"):
        fns.append(
            os.path.join(
                os.getenv("DESI_TARGET"),
                "catalogs",
                "gaiadr2",
                "0.47.0",
                "targets",
                "cmx",
                "resolve",
                "supp",
            )
        )

    # AR fixing paths...
    for i in range(len(fns)):
        fn = fns[i]
        fn = fn.replace(
            "DESIROOT/target/catalogs/mtl/1.1.1/mtl/main/ToO/ToO.ecsv",
            "{}/mtl/main/ToO/ToO.ecsv".format(os.getenv("DESI_SURVEYOPS")),
        )
        fn = fn.replace(
            "DESIROOT/survey/ops/staging/mtl/main/ToO/ToO.ecsv",
            "{}/mtl/main/ToO/ToO.ecsv".format(os.getenv("DESI_SURVEYOPS")),
        )
        # AR case where the path is generically defined, but not used (no mtl for sv1)
        fn = fn.replace(
            "DESIROOT/survey/ops/surveyops/trunk/mtl/sv1/ToO/ToO.ecsv",
            "-",
        )
        # AR case where the path is generically defined, but not used (no mtl for sv1)
        fn = fn.replace(
            "DESIROOT/survey/ops/surveyops/trunk/mtl/sv2/ToO/ToO.ecsv",
            "-",
        )
        fn = fn.replace(
            "DESIROOT/target/catalogs/dr9/0.58.0/targets/main/secondary/dark/maintargets-dark-secondary.fits",
            "DESIROOT/target/catalogs/dr9/0.58.0/targets/main/secondary/dark/targets-dark-secondary.fits",
        )
        # AR case where the path is generically defined, but not used (81096-81099; no secondaries for sv2)
        fn = fn.replace(
            "DESIROOT/target/catalogs/dr9/0.53.0/targets/sv2/secondary/dark/sv2targets-dark-secondary.fits",
            "-",
        )
        if name == "scnd":
            txt = os.path.sep.join(fn.split(os.path.sep)[-3:])
            if txt in ["sv1/secondary/bright", "sv1/secondary/dark"]:
                prog = fn.split(os.path.sep)[-1]
                fn = os.path.join(fn, "sv1targets-{}-secondary.fits".format(prog))
        fn = fn.replace("DESIROOT", os.getenv("DESI_ROOT"))
        fn = fn.replace("/data/target", os.getenv("DESI_TARGET"))
        fn = fn.replace(
            "/data/afternoon_planning/surveyops/trunk", os.getenv("DESI_SURVEYOPS")
        )
        fn = fn.replace(
            "/global/cscratch1/sd/adamyers/gaiadr2/1.3.0.dev5218/targets/main/resolve/backup",
            "{}/catalogs/gaiadr2/2.2.0/targets/main/resolve/backup".format(
                os.getenv("DESI_TARGET")
            ),
        )
        if fn != fns[i]:
            log.info(
                "{:06d}\t{}\treplacing {} by {}".format(
                    tileid, name.upper().ljust(11), fns[i], fn
                ),
            )
        fns[i] = fn
    for fn in fns:
        log.info("{:06d}\t{}\t{}".format(tileid, name.upper().ljust(11), fn))
        if fn != "-":
            if (fn[-5:] == ".fits") | (fn[-5:] == ".ecsv"):
                if not os.path.isfile(fn):
                    msg = "{:06d}\t{}\tmissing {}".format(
                        tileid, name.upper().ljust(11), fn
                    )
                    log.error(msg)
                    raise IOError(msg)
            else:
                if not os.path.isdir(fn):
                    msg = "{:06d}\t{}\tmissing {}".format(
                        tileid, name.upper().ljust(11), fn
                    )
                    log.error(msg)
                    raise IOError(msg)
    return fns


_scnd_radec_cache = {}


def read_scnd_radec(fn, columns=["RA", "DEC"]):
    """
    Reads the RA,DEC columns of a secondary *.fits catalog.

    Args:
        fn: full path to the secondary catalog (str)

    Returns:
        The secondary catalog with RA,DEC columns (Table() array)

    Notes:
        To speed up, the catalog is cached.
    """
    global _scnd_radec_cache
    try:
        log.info("using cached {}".format(fn))
        return _scnd_radec_cache[fn]
    except KeyError:
        pass
    log.info("reading {}".format(fn))
    assert fn[-5:] == ".fits"
    _scnd_radec_cache[fn] = Table(fitsio.read(fn, columns=["RA", "DEC"]))
    return _scnd_radec_cache[fn]


def get_tiles(fafn):
    """
    Generate of "fake" tiles array for a given fiberassign file.

    Args:
        fafn: path to the fiberassign file (fiberassign-TILEID.fits.gz) (str)

    Returns:
        Table() array with TILEID,RA,DEC,IN_DESI,OBSCONDITIONS (Table() array)

    Notes:
        This is to be used in e.g. desitarget.io.read_targets_in_tiles()
    """
    fahdr = fits.getheader(fafn, 0)
    tiles = Table()
    tiles["TILEID"] = [fahdr["TILEID"]]
    tiles["RA"] = [fahdr["TILERA"]]
    tiles["DEC"] = [fahdr["TILEDEC"]]
    tiles["IN_DESI"] = 1  # AR forcing 1, irrespective of location
    tiles["OBSCONDITIONS"] = obsconditions.mask("DARK|GRAY|BRIGHT|BACKUP")
    return tiles


def get_desitarget_fafn_rows(fafn, ras, decs):
    """
    Returns the rows of (ras,decs) falling inside the fafn tile.

    Args:
        fafn: path to the fiberassign file (fiberassign-TILEID.fits.gz) (str)
        ras: R.A. values (np.array())
        decs: Dec. values (np.array())

    Returns:
        rows: indexes of ras,decs (list of ints)

    """
    tiles = get_tiles(fafn)
    nside, nest = pixarea2nside(7.0), True
    pixlist = tiles2pix(nside, tiles=tiles)
    pixs = hp.ang2pix(nside, np.radians((90.0 - decs)), np.radians(ras), nest=nest)
    ii = np.where(np.in1d(pixs, pixlist))[0]
    jj = is_point_in_desi(tiles, ras[ii], decs[ii])
    rows = ii[jj]
    return rows


def read_fafn_static_desitarget_name(fafn, name):
    """
    Reads a desitarget static catalog for a given fiberassign file.

    Args:
        fafn: path to the fiberassign file (fiberassign-TILEID.fits.gz) (str)
        name: "targ", "scnd", "too", "sky", or "gfa" (str)

    Returns:
        d: targets catalog covered by the tile (Table() array)

    Notes:
        If several "inputs" (e.g. SCND and SCND2 in fafn header), the
            catalogs are concatenated into one.
    """
    tileid = get_tileid_from_fafn(fafn)
    fns = get_static_desitarget_fns(fafn, name)
    fns = [fn for fn in fns if fn != "-"]
    if len(fns) == 0:
        d = None
    else:
        # AR fake tiles file
        tiles = get_tiles(fafn)
        # AR
        ds = []
        for fn in fns:
            log.info(
                "{:06d}\t{}\treading {}".format(tileid, name.upper().ljust(11), fn)
            )
            # AR fn is a catalog file
            if fn[-5:] == ".fits":
                assert name == "scnd"
                radecs = read_scnd_radec(fn)
                rows = get_desitarget_fafn_rows(fafn, radecs["RA"], radecs["DEC"])
                d = Table(fitsio.read(fn, rows=rows))
            elif fn[-5:] == ".ecsv":
                d = Table.read(fn)
                rows = get_desitarget_fafn_rows(fafn, d["RA"], d["DEC"])
                d = d[rows]
            # AR fn is a folder
            else:

                d = Table(read_targets_in_tiles(fn, tiles=tiles, quick=True))
            ds.append(d)
        d = vstack(ds)
    return d


def read_fafn_static_desitarget_names(fafn, names=all_names):
    """
    Reads all the desitarget static catalogs for a given fiberassign file.

    Args:
        fafn: path to the fiberassign file (fiberassign-TILEID.fits.gz) (str)
        names (optional, defaults to ["targ", "scnd", "too", "sky", or "gfa"]:
            list of fafn header keywords, restricted to "targ", "scnd", "too", "sky", or "gfa"
            (list of strs)

    Returns:
        ds: dictionary of targets catalog covered by the tile, with one key
            for each element of names, containing a Table() array (dict)
    """

    tileid = get_tileid_from_fafn(fafn)
    # AR fiberassign header
    fahdr = fits.getheader(fafn, 0)

    # AR desitarget input files
    ds = {}
    for name in names:
        if name in fahdr:
            d = read_fafn_static_desitarget_name(fafn, name)
            if d is not None:
                ds[name] = d
        else:
            log.warning(
                "{:06d}\t{}\tnot in header".format(tileid, name.upper().ljust(11)),
            )
    return ds


def get_dith_infos(ext, name, fafn, fa_name):
    """
    Obtain per-row information if dithering is applied or not for a
        fiberassign file.

    Args:
        ext: fiberassign extension (just used for logging) (str)
        name: "targ", "scnd", "too", "sky", or "gfa" (str)
        fafn: path to the fiberassign file (fiberassign-TILEID.fits.gz) (str)
        fa_name: Table() array from the ext extension of fiberassign (Table() array)

    Returns:
        isdith: dithering is applied? (array of bools)
    """

    isdith = np.zeros(len(fa_name), dtype=bool)

    # AR the dithering is for name = "targ" only
    # AR identifying the dithered rows
    fahdr = fits.getheader(fafn, 0)
    if (name == "targ") & (fahdr["FAFLAVOR"] in ["dithprec", "dithlost", "dithfocus"]):
        h = fits.open(fafn)
        if "EXTRA" in [h[i].header["EXTNAME"] for i in range(1, len(h))]:
            fa_undith = h["EXTRA"].data
            sel = np.in1d(fa_name["TARGETID"], fa_undith["TARGETID"])
            isdith[sel] = True
    log.info(
        "{:06d}\t{}\t{}\tdith: found {} rows".format(
            fahdr["TILEID"], ext.ljust(11), name.upper(), isdith.sum()
        )
    )

    return isdith


def get_pmcorr_infos(ext, name, fahdr, fa_name, d_name):
    """
    Obtain per-row information if proper-motion correction is applied or not for a
        fiberassign file.

    Args:
        ext: fiberassign extension (just used for logging) (str)
        name: "targ", "scnd", "too", "sky", or "gfa" (str)
        fahdr: fiberassign 0-extension header (header array)
        fa_name: Table() array from the ext extension of fiberassign (Table() array)
        d_name: Table() from the static desitarget catalog, row-matched to fa_name
            (Table() array)

    Returns:
        ispmcorr: proper-motion is applied? (array of bools)
        iscmx33std: 80615-specific information, is it a standard star (array of bools)

    Notes:
        iscmx33std will be all False, except for TILEID=80615.
    """
    # AR for cmx,sv1 we did correct for proper-motion,
    # AR    except for:
    # AR    - fba_cmx early dith* tiles (targ and gfa at different times)
    # AR    - cmxm33 (80615) for standards which are not science targets
    # AR we identify here objects passing the criterion
    # AR remarks:
    # AR - hard-coding gaia/dr2 here (i.e. not using gaia/edr3)
    # AR - we use the desitarget values, as fa values can be bugged
    # AR - in fba_sv1, for "scnd", we pm-corrected for all rows with REF_EPOCH>0
    # AR - we updated REF_EPOCH for all rows (except for cmxm33 stds)
    pmcorr = False
    ispmcorr = np.zeros(len(d_name), dtype=bool)
    iscmx33std = np.zeros(len(d_name), dtype=bool)
    if (fahdr["FA_SURV"] in ["cmx", "sv1"]) & (name != "sky"):
        if fahdr["FAFLAVOR"] in ["dithprec", "dithlost", "dithfocus"]:
            if (name == "gfa") & (fahdr["TILEID"] > 80171):
                pmcorr = True
            if (name == "targ") & (fahdr["TILEID"] > 80604):
                pmcorr = True
        else:
            pmcorr = True

    if pmcorr:
        if name == "scnd":
            ispmcorr = d_name["REF_EPOCH"] > 0
        else:
            gaia_aens, gaia_gs = (
                d_name["GAIA_ASTROMETRIC_EXCESS_NOISE"],
                d_name["GAIA_PHOT_G_MEAN_MAG"],
            )
            ispmcorr = gaia_psflike(gaia_aens, gaia_gs, dr="dr2")
            # AR exluding cmxm33 standards which are not science targets
            if (name == "targ") & (fahdr["FAFLAVOR"] == "cmxm33"):
                from desitarget.cmx.cmx_targetmask import cmx_mask

                std_mask = cmx_mask.mask("SV0_WD|STD_BRIGHT")
                sci_mask = cmx_mask.mask(
                    "SV0_WD|M33_H2PN|M33_GC|M33_QSO|M33_M33cen|M33_M33out|SV0_QSO|SV0_LRG|SV0_ELG"
                )
                iscmx33std = (d_name["CMX_TARGET"] & std_mask) > 0
                iscmx33std &= (d_name["CMX_TARGET"] & sci_mask) == 0
                ispmcorr &= ~iscmx33std

        # AR verify that in that case, all REF_EPOCH values have been set to the same value
        # AR (exclude the cmxm33 standards)
        if "REF_EPOCH" in fa_name.dtype.names:
            n = np.unique(fa_name["REF_EPOCH"][~iscmx33std]).size
            if n != 1:
                msg = "{:06d}\t{}\t{}\tfound {} unique values of REF_EPOCH; expected one".format(
                    fahdr["TILEID"], ext.ljust(11), name.upper(), n
                )
                log.error(msg)
                raise ValueError(msg)
    log.info(
        "{:06d}\t{}\t{}\tpmcorr: found {}/{} rows".format(
            fahdr["TILEID"], ext.ljust(11), name.upper(), ispmcorr.sum(), ispmcorr.size
        )
    )
    if "REF_EPOCH" in d_name.dtype.names:
        log.info(
            "{:06d}\t{}\t{}\tpmcorr: REF_EPOCH updated for {}/{} rows".format(
                fahdr["TILEID"],
                ext.ljust(11),
                name.upper(),
                (~iscmx33std).sum(),
                ispmcorr.size,
            )
        )
    return ispmcorr, iscmx33std


def get_fbaonly_keys():
    """
    Returns the columns defined by fiberassign which are not present in the desitarget catalogs.

    Args:
        None

    Returns:
        keys: dictionary, with one key per fiberassign extension, containing the list of columns (dict)

    """
    keys = {
        "FIBERASSIGN": [
            "FIBER",
            "LOCATION",
            "FIBERSTATUS",
            "LAMBDA_REF",
            "PETAL_LOC",
            "DEVICE_LOC",
            "DEVICE_TYPE",
            "FA_TARGET",
            "FA_TYPE",
            "FIBERASSIGN_X",
            "FIBERASSIGN_Y",
            "PLATE_RA",
            "PLATE_DEC",
        ]
        + ["OBJTYPE"],
        "SKY_MONITOR": [
            "FIBER",
            "LOCATION",
            "FA_TARGET",
            "FA_TYPE",
            "FIBERASSIGN_X",
            "FIBERASSIGN_Y",
            "FIBERSTATUS",
            "PETAL_LOC",
            "DEVICE_LOC",
            "PRIORITY",
        ],
        "GFA_TARGETS": ["ETC_FLAG", "FOCUS_FLAG", "GFA_LOC", "GUIDE_FLAG"],
        "TARGETS": ["FA_TARGET", "FA_TYPE", "PLATE_RA", "PLATE_DEC"],
    }
    return keys


# AR we rename some columns to match the original name
# AR see fiberassign.assign merged_fiberassign_swap
# AR APFLUX_... disappeared at some point in the fiberassign files...
def fa_key_backswap(fa, ext):
    """
    Renaming some columns to match the original name in the static desitarget catalogs.

    Args:
        fa: Table() array for a fiberassign extension data (Table() array)
        ext: fiberassign extension from which comes fa (str)

    Returns:
        fa: same as input, but with modified column names.
    """
    if ext in ["FIBERASSIGN", "SKY_MONITOR", "GFA_TARGETS"]:
        fa["TARGET_RA"].name, fa["TARGET_DEC"].name = "RA", "DEC"
        if ext == "GFA_TARGETS":
            fa["TARGET_RA_IVAR"].name, fa["TARGET_DEC_IVAR"].name = (
                "RA_IVAR",
                "DEC_IVAR",
            )
    if ext in ["TARGETS"]:
        for key in [
            "APFLUX_G",
            "APFLUX_R",
            "APFLUX_Z",
            "APFLUX_IVAR_G",
            "APFLUX_IVAR_R",
            "APFLUX_IVAR_Z",
        ]:
            if key in fa.dtype.names:
                fa[key].name = key.replace("APFLUX", "FIBERFLUX")
    return fa


def get_names_to_check(ext):
    """
    For a given fiberassign extension, returns the list of names to check.

    Args:
        ext: fiberassign extension (str)

    Returns:
        ext_names: list of names to check (list of strs)

    Notes:
        Names are picked form "targ", "scnd", "too", "sky", "gfa".
    """
    if ext in ["FIBERASSIGN", "TARGETS"]:
        ext_names = ["targ", "sky", "scnd", "too"]
    elif ext in ["SKY_MONITOR"]:
        ext_names = ["sky"]
    elif ext in ["GFA_TARGETS"]:
        ext_names = ["gfa"]
    else:
        raise ValueError("unexpected ext={}; exiting".format(ext))
    return ext_names


# AR black keys
def get_black_keys(name, fafn):
    """
    Returns a list of keys to ignore for fiberassign patching diagnosis.

    Args:
        name: "targ", "scnd", "too", "sky", or "gfa" (str)
        fafn: path to the fiberassign file (fiberassign-TILEID.fits.gz) (str)

    Returns:
        black_keys: list of columns to ignore (list of strs)

    Notes:
        This function is mostly handling the chaotic fiberassign
            evolution in the first months of DESI.
        The PRIORITY,PRIORITY_INIT,NUMOBS_INIT columns are always
            black-listed, as they are fa-sensitive columns.
    """
    # AR check all keys
    # AR but:
    # AR - PRIORITY,PRIORITY_INIT,NUMOBS_INIT: fa-sensitive keys
    # AR        (and often re-set by fiberassign, e.g. cmx/sv1, sv3/main mtl)
    # AR - OBSCONDITIONS: re-set by fiberassign
    # AR - FLUX_R against "gfa" (only for GFA_TARGETS), as fiberassign modifies it
    # AR - EBV against "sky", "scnd" and "too", as those do not have EBV
    # AR        (fiberassign fills the EBV value since fiberassign/5.4.0)
    # AR - SUBPRIORITY, when it was overwritten by fiberassign:
    # AR    - FA_VER < 5.0.0 (fixed in fba_launch in fiberassign/5.0.0)
    # AR    - FA_SURV != "main" (the fix did not apply for fba_cmx/fba_sv1
    # AR    - TILEIDS=82401-82409, designed with adamyers/gaiadr2/1.3.0.dev5218, different SUBPRIORITY...
    # AR OBSCONDITIONS:
    #       - reset in fba_cmx/fba_sv1 for targ,scnd
    #       - (can be) different in mtl vs. static desitarget...
    fahdr = fits.getheader(fafn, 0)
    black_keys = []
    black_keys += [
        "PRIORITY",
        "PRIORITY_INIT",
        "NUMOBS_INIT",
    ]
    if name in ["sky", "scnd", "too"]:
        black_keys += ["EBV"]
    if name == "gfa":
        black_keys += ["FLUX_R"]
    if (
        (fahdr["FA_VER"] < "5.0.0")
        | (fahdr["FA_SURV"] != "main")
        | (
            fahdr["TILEID"]
            in [82401, 82402, 82403, 82404, 82405, 82406, 82407, 82408, 82409]
        )
    ):
        black_keys += ["SUBPRIORITY"]
    if name in ["targ", "scnd"]:
        if (
            (fahdr["FA_SURV"] in ["cmx", "sv1"])
            | ((name == "targ") & ("MTL" in fahdr))
            | ((name == "scnd") & ("SCNDMTL" in fahdr))
        ):
            black_keys += ["OBSCONDITIONS"]
    black_keys = np.unique(black_keys)
    return black_keys


def arrays_equal(key, tileid, a, b):
    """
    Function to test for equality between arrays -- both being non-finite (eg NaN) counts as equal.

    Args:
        a: 1D np.array()
        b: 1D np.array()
    """
    bothnan = False
    try:
        bothnan = np.logical_not(np.isfinite(a)) * np.logical_not(np.isfinite(b))
    except:
        pass
    # AR PARALLAX,PMRA,PMDEC
    # AR - for some 'specialm31', 'sv1dc3r2' tiles, format change in desitarget (dtype= ">f8" or ">f4")
    # AR - else discrepancy in formatting between fiberassign and desitarget
    # AR        (makes difference <= 5e9; we pick 1e8 for simplicity)
    if key in ["PARALLAX", "PMRA", "PMDEC"]:
        if tileid in [
            80971,
            80972,
            80973,
            80974,
            80975,
            80976,
            82634,
            82635,
        ]:
            eq = np.abs(a - b) < 1e-6
        else:
            eq = np.abs(a - b) < 1e-8
    else:
        eq = a == b
    return np.logical_or(eq, bothnan)


# AR shall not use Table.read() ! (it messes up things... with "masking" some columns)
def diagnose_values_fafn(fafn):
    """
    Execute diagnostic of discrepant values between the fiberassign file
        and the input static desitarget catalogs.

    Args:
        fafn: path to the fiberassign file (fiberassign-TILEID.fits.gz) (str)

    Returns:
        d: a Table() array with the discrepant values (Table() array)

    Notes:
        The comparison ignores some expected/understood discrepancies.
        The column content of the output is:
            "FAFN": fiberassign file name
            "TILEID": tileid
            "FAFLAVOR": faflavor
            "FA_VER": fiberassign version
            "PMTIME": time when the tile was designed
            "LASTNIGHT": last night of observation for that tile (-99 if not observed)
            "EXTENSION": fiberassign extension
            "TARGETID": targetid
            "ORIGFN": path to the static desitarget catalog
            "KEY": column name
            "ORIGVAL": value in the desitarget catalog, converted to string ("none" if the column is not present)
            "FAVAL": value in the fiberassign file
    """

    # AR tileid
    tileid = get_tileid_from_fafn(fafn)
    log.info("{:06d}\t{}\t{}".format(tileid, "FAFN".ljust(11), fafn))

    # AR
    fbaonly_keys = get_fbaonly_keys()

    #
    myd = {
        key: []
        for key in ["EXTENSION", "TARGETID", "ORIGFN", "KEY", "ORIGVAL", "FAVAL"]
    }

    # AR fiberassign header
    fahdr = fits.getheader(fafn, 0)

    # AR desitarget input files
    ds = read_fafn_static_desitarget_names(fafn, names=all_names)

    # AR check differences for each ext
    for ext in all_exts:
        fa = Table(fitsio.read(fafn, ext=ext))
        # AR we restrict to TARGETID > 0
        # AR that removes:
        # AR - unassigned / stuck sky fibers
        # AR - TARGETID=-1 in GFA_TARGETS
        fa = fa[fa["TARGETID"] > 0]
        log.info(
            "{:06d}\t{}\t{} rows with TARGETID>0".format(tileid, ext.ljust(11), len(fa))
        )

        # AR we rename some columns to match the original name
        # AR see fiberassign.assign merged_fiberassign_swap
        # AR APFLUX_... disappeared at some point in the fiberassign files...
        fa = fa_key_backswap(fa, ext)

        # AR verify that TARGETIDs are unique
        assert np.unique(fa["TARGETID"]).size == len(fa)

        # AR names to check
        ext_names = get_names_to_check(ext)

        # AR now cutting on existing files
        ext_names = [name for name in ext_names if name in ds]

        # AR to record what rows are matched
        ii_match = []

        for name in ext_names:

            # AR verify that TARGETIDs are unique
            # AR exception for 80865, where three TARG catalogs are used
            try:
                assert np.unique(ds[name]["TARGETID"]).size == len(ds[name]["TARGETID"])
            except AssertionError:
                _, ii = np.unique(ds[name]["TARGETID"], return_index=True)
                ds[name] = ds[name][ii]
                log.warning(
                    "{:06d}\t{}\t{}\tsome duplicates TARGETID; picking a unique list..".format(
                        tileid, ext.ljust(11), name.upper()
                    ),
                )

            # AR cut catalogs on matched objects
            ii, iid = match(fa["TARGETID"], ds[name]["TARGETID"])
            log.info(
                "{:06d}\t{}\t{}\t{}/{} matched rows".format(
                    tileid, ext.ljust(11), name.upper(), ii.size, len(fa)
                )
            )
            ii_match += ii.tolist()

            if ii.size > 0:

                fa_name, d_name = fa[ii], ds[name][iid]
                keys = np.array(fa_name.dtype.names)
                extname_diffs = {key: -99 for key in keys}

                # AR proper-motion correction infos
                ispmcorr, iscmx33std = get_pmcorr_infos(
                    ext, name, fahdr, fa_name, d_name
                )

                # AR dithering infos
                isdith = get_dith_infos(ext, name, fafn, fa_name)

                # AR check all keys, except black_keys
                black_keys = get_black_keys(name, fafn)
                check_keys = [
                    key
                    for key in keys
                    if key not in black_keys and key not in fbaonly_keys[ext]
                ]
                log.info(
                    "{:06d}\t{}\t{}\tblack_keys={}".format(
                        tileid, ext.ljust(11), name.upper(), ",".join(black_keys)
                    )
                )
                log.info(
                    "{:06d}\t{}\t{}\tcheck_keys={}".format(
                        tileid, ext.ljust(11), name.upper(), ",".join(check_keys)
                    )
                )

                for key in check_keys:
                    # AR common key?
                    if key in d_name.dtype.names:
                        ignore = np.zeros(len(fa_name), dtype=bool)
                        # AR PMRA,PMDEC:
                        # - ignore if nan in desitarget and 0 in fa
                        if key in ["PMRA", "PMDEC"]:
                            ignore = (~np.isfinite(d_name[key])) & (fa_name[key] == 0.0)
                        # AR RA,DEC:
                        # AR - check only for rows with:
                        # AR    - no proper-motion correction
                        # AR    - no dithering
                        if key in ["RA", "DEC"]:
                            ignore = (ispmcorr) | (isdith)
                        # AR REF_EPOCH:
                        # AR - if some pmcorr, REF_EPOCH updated for all rows (except for iscmx33std)
                        # AR - ignore cases where fiberassign changed 0 to 2015.5
                        # AR - also ignore few sv3 tiles, as $DESI_SURVEYOPS/mtl/sv3/ToO/ToO.ecsv
                        # AR        has been overwritten (orignal rows had REF_EPOCH=0,
                        # AR        changed to REF_EPOCH=2015.5 in fiberassign; but then the ToO.ecsv
                        # AR        got overwritten, with updated rows having REF_EPOCH=2000)
                        # AR - same story for 80980-80981 and $DESI_SURVEYOPS/mtl/main/ToO/ToO.ecsv
                        if key == "REF_EPOCH":
                            if ispmcorr.sum() > 0:
                                ignore = (~iscmx33std).copy()
                            else:
                                # ignore = ispmcorr.copy()
                                ignore = (d_name[key] == 0) & (fa_name[key] == 2015.5)
                                if (name == "too") & (
                                    tileid
                                    in [466, 518, 519, 532, 533, 534, 80980, 80981]
                                ):
                                    ignore |= (d_name[key] == 2000) & (
                                        fa_name[key] == 2015.5
                                    )
                        # AR
                        iseq = arrays_equal(key, tileid, d_name[key], fa_name[key])
                        sel = (~iseq) & (~ignore)
                    # AR if not, then it should be zero (we already excluded fba-only keys)
                    else:
                        if isinstance(fa_name[key][0], str):
                            sel = fa_name[key] != ""
                        else:
                            sel = fa_name[key] != 0
                    n = sel.sum()
                    extname_diffs[key] = n
                    # log.info("{:06d}\t{}\t{}\tkey = {} -> {} discrepancies".format(tileid, ext.ljust(11), name, key, n))
                    if n > 0:
                        myd["EXTENSION"] += [ext for j in range(n)]
                        myd["TARGETID"] += fa_name["TARGETID"][sel].tolist()
                        tmpfn = ",".join(get_static_desitarget_fns(fafn, name))
                        myd["ORIGFN"] += [tmpfn for j in range(n)]
                        myd["KEY"] += [key for j in range(n)]
                        if key in d_name.dtype.names:
                            myd["ORIGVAL"] += d_name[key][sel].astype(str).tolist()
                        else:
                            myd["ORIGVAL"] += ["none" for j in range(n)]
                        myd["FAVAL"] += fa_name[key][sel].astype(str).tolist()
                vals = np.unique([extname_diffs[key] for key in keys])
                vals = [val for val in vals if val not in [-99, 0]]
                for val in vals:
                    txt = "ndiff={}".format(str(val).ljust(4))
                    keystxt = ",".join(
                        np.sort([key for key in keys if extname_diffs[key] == val])
                    )
                    log.info(
                        "{:06d}\t{}\t{}\t{}\t{}".format(
                            tileid, ext.ljust(11), name.upper(), txt, keystxt
                        )
                    )

        # AR there should not be duplicates in ii_match
        n_match = len(ii_match)
        if np.unique(ii_match).size != n_match:
            msg = "{:06d}\t{}\tmatches have duplicated {} TARGETIDs".format(
                tileid, name, n_match - np.unique(ii_match).size
            )
            log.error(msg)
            raise ValueError(msg)

        # AR verify that all rows have been matched
        log.info(
            "{:06d}\t{}\tall_names\t{}/{} matched rows".format(
                tileid, ext.ljust(11), n_match, len(fa)
            )
        )
        if n_match != len(fa):
            msg = "{}:06d\t{}\tn_match={} != len(fa)={}".format(
                tileid, ext.ljust(11), n_match, len(fa)
            )
            # AR mismatch happens for 82405,82406; we do not track it down, just ignore..
            if tileid in [
                82401,
                82402,
                82403,
                82404,
                82405,
                82406,
                82407,
                82408,
                82409,
            ]:
                log.warning(msg)
                log.warning(
                    "{:06d}\t{}\tignoring mismatch, as tile designed with /global/cscratch1/sd/adamyers/gaiadr2/1.3.0.dev5218".format(
                        tileid, ext.ljust(11)
                    ),
                )
            else:
                log.error(msg)
                raise ValueError(msg)

    # AR store in a table
    d = Table()
    n = len(myd["TARGETID"])
    log.info("{:06d}\t{}\tall_names\tndiff={}".format(tileid, "all_exts".ljust(11), n))
    if n > 0:
        d["FAFN"] = [fafn for j in range(n)]
        for key in ["TILEID", "FAFLAVOR", "FA_VER", "PMTIME"]:
            if key not in fahdr:
                d[key] = "-"  # AR PMTIME
            else:
                d[key] = [fahdr[key] for j in range(n)]
        tiles = Table.read(
            os.path.join(
                os.getenv("DESI_ROOT"), "spectro", "redux", "daily", "tiles-daily.csv"
            )
        )
        ii = np.where(tiles["TILEID"] == tileid)[0]
        keys = ["LASTNIGHT"]
        if ii.size == 0:
            for key in keys:
                d[key] = -99
        else:
            for key in keys:
                d[key] = [tiles[key][ii[0]] for j in range(n)]
        for key in myd:
            d[key] = myd[key]

    return d


def diagnose_columns_fafn(fafn, ref_tileid=3001):
    """
    Execute diagnostic of the column content of the fiberassign file, with
        comparing with a fiducial Main tile.

    Args:
        fafn: path to the fiberassign file (fiberassign-TILEID.fits.gz) (str)
        ref_tileid (optional, defaults to 3001): Main tileid used as a reference (int)

    Returns:
        d: a Table() array with the column diagnosis.

    Notes:
        The column content of the output is:
            "TILEID": tileid
            "FAFLAVOR": faflavor
            "FA_VER": fiberassign version
            "EXTENSION": fiberassign extension
            "MISS": comma-separated list of missing columns
            "EXTRA": comma-separated list of extra columns
    """

    # AR survey-specific keys
    main_dt_keys = ["DESI_TARGET", "BGS_TARGET", "MWS_TARGET", "SCND_TARGET"]
    main2surv_dtkeys = {
        "cmx": {
            key: "CMX_TARGET" if key == "DESI_TARGET" else None for key in main_dt_keys
        },
        "main": {key: key for key in main_dt_keys},
    }
    for surv in ["sv1", "sv2", "sv3"]:
        main2surv_dtkeys[surv] = {
            key: "{}_{}".format(surv.upper(), key) for key in main_dt_keys
        }

    # AR survey for fafn
    h = fits.open(fafn)
    _, _, surv = main_cmx_or_sv(h["FIBERASSIGN"].data)
    if surv not in main2surv_dtkeys:
        msg = "{}\tunexpected surv={}".format(fafn, surv)
        log.error(msg)
        raise ValueError(msg)

    # AR reference set of columns
    ref_tileidpad = "{:06d}".format(ref_tileid)
    ref_fafn = os.path.join(
        os.getenv("DESI_TARGET"),
        "fiberassign",
        "tiles",
        "trunk",
        ref_tileidpad[:3],
        "fiberassign-{}.fits.gz".format(ref_tileidpad),
    )
    ref_fahdr = fits.getheader(ref_fafn, 0)
    if ref_fahdr["FA_SURV"] != "main":
        msg = "ref_tileid={} (FA_SURV={}) is not a Main tile".format(
            ref_tileid, ref_fahdr["FA_SURV"]
        )
        log.error(msg)
        raise ValueError(msg)

    ref_h = fits.open(ref_fafn)
    keys = {}
    for ext in all_exts:
        keys[ext] = []
        for key in ref_h[ext].columns.names:
            if key in main_dt_keys:
                if main2surv_dtkeys[surv][key] is not None:
                    keys[ext].append(main2surv_dtkeys[surv][key])
                # AR DESI_TARGET,BGS_TARGET,MWS_TARGET are always present in FIBERASSIGN and TARGETS because of skies...
                if (
                    (ext in ["FIBERASSIGN", "TARGETS"])
                    & (surv != "main")
                    & (key in ["DESI_TARGET", "BGS_TARGET", "MWS_TARGET"])
                ):
                    keys[ext].append(key)
            else:
                keys[ext].append(key)

    #
    d = Table()
    for key in ["TILEID", "FA_VER", "FAFLAVOR"]:
        d[key] = [h[0].header[key] for ext in all_exts]
    d["EXTENSION"] = all_exts
    d["MISS"] = np.array(
        [
            ",".join([key for key in keys[ext] if key not in h[ext].columns.names])
            for ext in all_exts
        ]
    )
    d["EXTRA"] = np.array(
        [
            ",".join([key for key in h[ext].columns.names if key not in keys[ext]])
            for ext in all_exts
        ]
    )
    # AR replace "" by "-"
    d["MISS"][d["MISS"] == ""] = "-"
    d["EXTRA"][d["EXTRA"] == ""] = "-"

    return d


def get_patching_params(fn="patching_202210.yaml"):
    """
    Obtain the patching_root, fixcols, addcols, and populate_ebv patching parameters.

    Args:
        fn (optional, defaults to "patching_202210.yaml"): path to a .yaml file with patching_root, fixcols, addcols, populate_ebv (str)

    Returns:
        params: dictionary with the patching_root, fixcols, addcols, and populate_ebv patching parameters.

    Notes:
        fn: the code will first look for:
            - first look for: resource_filename("fiberassign", os.path.join("data", fn))
            - then look for: fn
        See fiberassign/data/patching_202210.yaml for the formatting.
    """
    myfn = resource_filename("fiberassign", os.path.join("data", fn))
    if not os.path.isfile(myfn):
        log.warning("no {}".format(myfn))
        myfn = fn
        if not os.path.isfile(myfn):
            log.warning("no {}".format(myfn))
            msg = "did not find {}".format(fn)
            log.error(msg)
            raise ValueError(msg)
    log.info("reading {}".format(myfn))
    with open(myfn, "r") as f:
        params = yaml.safe_load(f)
    f.close()
    return params


def patch(in_fafn, out_fafn, params_fn):
    """
    Repair data corruption in a fiberassign file.

    Args:
        in_fafn: full path to the fiberassign file to be repaired (str)
        out_fafn: full path where the repaired fiberassign file will be written, if any (str)
        params_fn: path to a .yaml file with fixcols, addcols, populate_ebv (str)

    Notes:
        Originally developed by D. Lang.
        Targets are matched on TARGETID, and values updated as necessary.
        If any rows are changed, then an updated file is written out, preserving all other HDUs.
        params_fn:
            - the code will first look for:
                - first look for: resource_filename("fiberassign", os.path.join("data", fn))
                - then look for: fn
            - see fiberassign/data/patching_202210.yaml for the formatting.
    """

    # AR safe
    assert(out_fafn != in_fafn)

    tileid = get_tileid_from_fafn(in_fafn)

    # AR default fixcols, addcols, populate_ebv
    params = get_patching_params(fn=params_fn)
    for key in ["patching_root", "fixcols", "addcols", "populate_ebv"]:
        log.info("{:06d}\tPatching with: {} = {}".format(tileid, key, params[key]))

    log.info("{:06d}\tin_fafn={}".format(tileid, in_fafn))
    log.info("{:06d}\tout_fafn={}".format(tileid, out_fafn))

    # AR read in_fafn
    F = fitsio.FITS(in_fafn)

    # AR read all the desitarget static catalogs
    ds = read_fafn_static_desitarget_names(in_fafn)

    patched_tables = {}
    added_tables = {}

    myd = {}
    for key in ["EXTENSION", "NAME", "KEY", "TARGETID", "OLDVAL", "NEWVAL"]:
        myd[key] = []

    # AR extensions to parse (in the "fiducial order")
    tmp_exts = list(params["fixcols"].keys())
    tmp_exts += list(params["addcols"].keys())
    tmp_exts += list(params["populate_ebv"].keys())
    exts = [ext for ext in all_exts if ext in tmp_exts]

    # AR loop on extensions
    for ext in exts:

        names = get_names_to_check(ext)

        ok_cols, patched_cols, added_cols = params["fixcols"][ext].copy(), [], []

        ff = F[ext]
        tab = ff.read()

        # AR adding columns?
        if ext in params["addcols"]:
            addtab = Table()
            cols, dtypes = [val[0] for val in params["addcols"][ext]], [val[1] for val in params["addcols"][ext]]
            for col, dtype in zip(cols, dtypes):
                if col in tab.dtype.names:
                    log.info("{:06d}\t{}\tColumn {}: already present".format(tileid, ext.ljust(11), col)) 
                else:
                    addtab[col] = np.zeros(len(tab), dtype=dtype)
                    added_cols.append(col)
                    log.info("{:06d}\t{}\tColumn {}: adding the column".format(tileid, ext.ljust(11), col)) 

        # AR we restrict the matching to positive TARGETIDs
        targetids = tab["TARGETID"]
        unq_targetids = np.unique(targetids[targetids > 0])
        log.info("{:06d}\t{}\tThere are {} TARGETID>0".format(tileid, ext.ljust(11), unq_targetids.size))

        # AR to record what rows are matched
        ii_match = []

        for name in names:

            if name not in ds:
                log.info("{:06d}\t{}\t{}\tno catalog".format(tileid, ext.ljust(11), name))
                continue
            # AR desitarget data
            T = ds[name]
            if len(T) == 0:
                log.info("{:06d}\t{}\t{}\tZero rows in catalog".format(tileid, ext.ljust(11), name))
                continue
            tids = T["TARGETID"]
            I = np.flatnonzero([t in unq_targetids for t in tids])
            log.info("{:06d}\t{}\t{}\tFound {}/{} matching TARGETIDs".format(tileid, ext.ljust(11), name, len(I), len(T)))
            if len(I) == 0:
                continue
            T = T[I]
            ii_match += np.where(np.in1d(unq_targetids, T["TARGETID"]))[0].tolist()

            # AR "targetids" are the ones from the FA file
            # AR "tids" are the ones from the desitarget files (always > 0)
            tid_map = dict([(tid, i) for i, tid in enumerate(T["TARGETID"])])
            I = np.array([tid_map.get(t, -1) for t in targetids])
            J = np.flatnonzero(I >= 0)
            I = I[J]
            # log.info("Checking {} matched TARGETIDs".format(len(I)))

            # AR loop on columns (the ones to be fixed or, the added ones)
            ext_cols = params["fixcols"][ext]
            if ext in params["addcols"]:
                ext_cols += [val[0] for val in params["addcols"][ext]]
            ext_cols = np.unique(ext_cols)
            for col in ext_cols:
                isaddcol = False
                if col not in tab.dtype.names:
                    if col not in addtab.dtype.names:
                        continue
                    else:
                        old = addtab[col][J]
                        isaddcol = True
                else:
                    old = tab[col][J]

                # AR case where the column is not present in desitarget
                if not col in T.dtype.names:
                    new = np.zeros_like(old)
                else:
                    new = T[col][I]

                # AR check "equality"
                eq = arrays_equal(col, tileid, old, new)
                if not np.all(eq):
                    if isaddcol:
                        log.info("{:06d}\t{}\t{}\tColumn {}: filling {} row(s)".format(tileid, ext.ljust(11), name, col, (~eq).sum()))
                    else:
                        if col not in patched_cols:
                            ok_cols.remove(col)
                            patched_cols.append(col)
                        log.info("{:06d}\t{}\t{}\tColumn {}: patching {} row(s)".format(tileid, ext.ljust(11), name, col, (~eq).sum()))
                    diff = np.flatnonzero(np.logical_not(eq))
                    myd["EXTENSION"] += [ext for x in range(diff.size)]
                    myd["NAME"] += [name for x in range(diff.size)]
                    myd["KEY"] += [col for x in range(diff.size)]
                    myd["TARGETID"] += tab["TARGETID"][J[diff]].tolist()
                    myd["NEWVAL"] += new[diff].astype(str).tolist()
                    if isaddcol:
                        myd["OLDVAL"] += addtab[col][J[diff]].astype(str).tolist()
                        addtab[col][J[diff]] = new[diff]
                    else:
                        myd["OLDVAL"] += tab[col][J[diff]].astype(str).tolist()
                        tab[col][J[diff]] = new[diff]

            # AR populate EBV for FIBERASSIGN?
            # AR note that J is already cut on TARGETID>0, so we are safe
            if ext in params["populate_ebv"]:
                if params["populate_ebv"][ext]:
                    col = "EBV"
                    old = tab[col][J]
                    diff = np.flatnonzero(old == 0)
                    if diff.size > 0:
                        if col not in patched_cols:
                            if col in ok_cols:
                                ok_cols.remove(col)
                            patched_cols.append(col)
                        ras, decs = tab["TARGET_RA"][J][diff], tab["TARGET_DEC"][J][diff]
                        assert((~np.isfinite(ras)).sum() == 0)
                        assert((~np.isfinite(decs)).sum() == 0)
                        cs = SkyCoord(ra=ras * units.deg, dec=decs * units.deg, frame="icrs")
                        ebvs = SFDMap(scaling=1).ebv(cs)
                        log.info("{:06d}\t{}\t{}\tColumn {}: populating {} row(s)".format(tileid, ext.ljust(11), name, col, diff.size))
                        myd["EXTENSION"] += [ext for x in range(diff.size)]
                        myd["NAME"] += [name for x in range(diff.size)]
                        myd["KEY"] += [col for x in range(diff.size)]
                        myd["TARGETID"] += tab["TARGETID"][J[diff]].tolist()
                        myd["OLDVAL"] += tab[col][J[diff]].astype(str).tolist()
                        myd["NEWVAL"] += ebvs.astype(str).tolist()
                        tab[col][J[diff]] = ebvs

        # AR store in patched/added_tables the fixed/added table
        patched_tables[ext] = tab
        if len(added_cols) > 0:
            added_tables[ext] = addtab

        # AR there should not be duplicates in ii_match
        n_match = len(ii_match)
        if np.unique(ii_match).size != n_match:
            msg = "{:06d}\t{}\tmatches have duplicated {} TARGETIDs".format(
                tileid, ext.ljust(11), n_match - np.unique(ii_match).size
            )
            log.error(msg)
            raise ValueError(msg)


        # AR verify that all rows have been matched
        log.info(
            "{:06d}\t{}\tall_names\t{}/{} matched rows".format(
                tileid, ext.ljust(11), n_match, len(unq_targetids)
            )
        )
        if n_match != len(unq_targetids):
            msg = "{:06d}\t{}\tn_match={} != len(fa)={}".format(
                tileid, ext.ljust(11), n_match, len(unq_targetids)
            )

        # AR columns summary
        log.info("{:06d}\t{}\tColumns ok:\t{}".format(tileid, ext.ljust(11), ",".join(ok_cols)))
        log.info("{:06d}\t{}\tColumns patched:\t{}".format(tileid, ext.ljust(11), ",".join(patched_cols)))
        log.info("{:06d}\t{}\tColumns added:\t{}".format(tileid, ext.ljust(11), ",".join(added_cols)))


    if len(myd["TARGETID"]) > 0:
        # DL make sure output directory exists
        outdir = os.path.dirname(out_fafn)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # DL write out output file, leaving other HDUs unchanged.
        f, tempout = tempfile.mkstemp(dir=outdir, suffix=".fits")
        os.close(f)
        Fout = fitsio.FITS(tempout, "rw", clobber=True)
        # DL/AR fitsio will add its own headers about the FITS format
        # DL/AR     so we do not propagate those
        # AR    note that some tiles (mostly dithers) have twice these comments
        # AR        so we only "remove them once"...
        # AR    also, we remove specifically those comments, as other meaningful
        # AR    information can be stored in the COMMENT keywords (see e.g. TILEID=80611)
        comment2trims = [
            "  FITS (Flexible Image Transport System) format is defined in 'Astronomy",
            "  and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H",
        ]
        for ext in F:
            extname = ext.get_extname()
            hdr = ext.read_header()
            data = ext.read()
            if extname == "PRIMARY":
                newhdr = fitsio.FITSHDR()
                istrim = {comment2trim : False for comment2trim in comment2trims}
                for r in hdr.records():
                    if (r["name"] == "COMMENT") & (r["comment"] in comment2trims):
                        if ~istrim[r["comment"]]:
                            istrim[r["comment"]] = True
                            continue
                    newhdr.add_record(r)
                hdr = newhdr
            if extname in patched_tables:
                # Swap in our updated FIBERASSIGN table!
                if extname in added_tables:
                    data = np.array(hstack([Table(patched_tables[extname]), added_tables[extname]]))
                else:
                    data = patched_tables[extname]
            Fout.write(data, header=hdr, extname=extname)
        Fout.close()
        os.rename(tempout, out_fafn)
        log.info("Wrote {}".format(out_fafn))
    else:
        log.info("No changes for tile {} file {}".format(tileid, in_fafn))

    # AR store in a table all the changes
    d = Table()
    n = len(myd["TARGETID"])
    log.info("{:06d}\t{}\tall_names\tndiff={}".format(tileid, "all_exts".ljust(11), n))
    if n > 0:
        fahdr = F["PRIMARY"].read_header()
        d["OLDFAFN"] = [in_fafn for j in range(n)]
        for key in ["TILEID", "FAFLAVOR", "FA_VER", "PMTIME"]:
            if key not in fahdr:
                d[key] = "-"  # AR PMTIME
            else:
                d[key] = [fahdr[key] for j in range(n)]
        tiles = Table.read(
            os.path.join(
                os.getenv("DESI_ROOT"), "spectro", "redux", "daily", "tiles-daily.csv"
            )
        )
        ii = np.where(tiles["TILEID"] == tileid)[0]
        keys = ["LASTNIGHT"]
        if ii.size == 0:
            for key in keys:
                d[key] = -99
        else:
            for key in keys:
                d[key] = [tiles[key][ii[0]] for j in range(n)]
        for key in myd:
            d[key] = myd[key]
        d.write(out_fafn.replace(".fits.gz", "-{}.ecsv".format(params["patching_root"])))
