#!/usr/bin/env python

import matplotlib

matplotlib.use("Agg")
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from matplotlib import gridspec
import healpy as hp
from glob import glob
from astropy import units
from astropy.coordinates import SkyCoord
from desitarget.sv1.sv1_targetmask import desi_mask
from desitarget.cmx.cmx_targetmask import cmx_mask
import fitsio
from pathlib import Path
import astropy.table
import scipy.ndimage
from astropy.table import Table
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
from argparse import ArgumentParser

# AR reading arguments
parser = ArgumentParser()
parser.add_argument(
    "--outdir",
    help="output directory",
    type=str,
    default=None,
    required=True,
    metavar="OUTDIR",
)
parser.add_argument(
    "--tiles",
    help="merge the TILEID-tiles.fits into one file? (y/n)",
    type=str,
    default="y",
    metavar="TILES",
)
parser.add_argument(
    "--exposures",
    help="do fits with per-exposure stats? (y/n)",
    type=str,
    default="y",
    metavar="EXPOSURES",
)
parser.add_argument(
    "--wiki", help="do wiki tables? (y/n)", type=str, default="y", metavar="WIKI"
)
parser.add_argument(
    "--skymap", help="do skymap png? (y/n)", type=str, default="y", metavar="SKYMAP"
)
parser.add_argument(
    "--depth", help="do r_depth png? (y/n)", type=str, default="y", metavar="DEPTH"
)
parser.add_argument(
    "--update",
    help="start from existing args.outroot+'sv1-exposures.fits' (y/n)",
    type=str,
    default="y",
    metavar="UPDATE",
)
args = parser.parse_args()
for kwargs in args._get_kwargs():
    print(kwargs)


# AR : all exposure depths routines are copied from DK, with minor modifications

# AR safe
if args.outdir[-1] != "/":
    args.outdir += "/"

# AR folders / files
sv1dir = "/global/cfs/cdirs/desi/users/raichoor/fiberassign-sv1/"
dailydir = "/global/cfs/cdirs/desi/spectro/redux/daily/"
pixwfn = "/global/cfs/cdirs/desi/target/catalogs/dr9/0.47.0/pixweight/sv1/resolve/dark/sv1pixweight-dark.fits"
desfn = os.path.join(sv1dir, "misc", "des_footprint.txt")
gfafn = np.sort(
    glob(
        "/global/cfs/cdirs/desi/users/ameisner/GFA/conditions/offline_all_guide_ccds_SV1-thru_20??????.fits"
    )
)[-1]

# AR output products
outfns = {}
outfns["tiles"] = os.path.join(args.outdir, "sv1-tiles.fits")
outfns["exposures"] = os.path.join(args.outdir, "sv1-exposures.fits")
outfns["wiki"] = os.path.join(args.outdir, "sv1-tables.wiki")
outfns["skymap"] = os.path.join(args.outdir, "sv1-skymap.png")
outfns["depth"] = {}
for flavshort in ["QSO+LRG", "ELG", "BGS+MWS"]:
    outfns["depth"][flavshort] = os.path.join(
        args.outdir, "sv1-depth-{}.png".format(flavshort.lower().replace("+", ""))
    )

# AR sv1 first night (for exposures search)
# AR if args.update="y":
# AR - assumes args.outdir+"sv1-exposures.fits" exists
# AR - recompute the last night of args.outdir+"sv1-exposures.fits", plus following nights
if args.update == "y":
    if not os.path.isfile(outfns["exposures"]):
        print("{} is missing; exiting".format(outfns["exposures"]))
        sys.exit()
    else:
        d = fits.open(outfns["exposures"])[1].data
        firstnight = d["NIGHT"].max().astype(str)
else:
    firstnight = "20201214"
print("firstnight = {}".format(firstnight))

# AR for the exposure and the wiki cases
targets = ["TGT", "SKY", "STD", "WD", "LRG", "ELG", "QSO", "BGS", "MWS"]
cmx_msks = [
    "TGT",
    "SKY",
    "STD",
    "SV0_WD",
    "SV0_LRG",
    "SV0_ELG",
    "SV0_QSO",
    "SV0_BGS",
    "SV0_MWS",
]
sv1_msks = ["TGT", "SKY", "STD", "STD_WD", "LRG", "ELG", "QSO", "BGS_ANY", "MWS_ANY"]
std_cmx_msks = ["SV0_WD", "STD_BRIGHT"]
std_sv1_msks = ["STD_WD", "STD_BRIGHT", "STD_FAINT"]

# AR short names for faflavor
flavdict = {
    "QSO+LRG": {"FAFLAVORS": ["cmxlrgqso", "sv1lrgqso"], "COLOR": "r"},
    "ELG": {"FAFLAVORS": ["cmxelg", "sv1elg"], "COLOR": "b"},
    "BGS+MWS": {"FAFLAVORS": ["cmxbgsmws", "sv1bgsmws"], "COLOR": "g"},
    "M33+Dark": {"FAFLAVORS": ["cmxm33"], "COLOR": "magenta"},
}
# AR field names - hand-written...
fielddict = {
    "XMM-LSS": [80605, 80606],
    "Lynx": [80607, 80608, 80613],
    "COSMOS": [80609, 80610],
    "Triangulum": [80611],
    "Eridanus": [80612],
    "Sextans": [80614],
    "Triangulum-CMX": [80615],
    "Pegasus1": [80616],
    "Pegasus2": [80617],
    "NGC2419": [80618],
    "UMajor": [80619, 80620, 80621],
    "LeoMinor": [80622, 80623],
}


# AR/DK settings for exposure depths
wmin, wmax, wdelta = 3600, 9824, 0.8
fullwave = np.round(np.arange(wmin, wmax + wdelta, wdelta), 1)
cslice = {"b": slice(0, 2751), "r": slice(2700, 5026), "z": slice(4900, 7781)}

# AR/DK exposure depths utilities
def load_spec_thru(path=os.getenv("DESIMODEL") + "/data/throughput/"):
    thru = {}
    for camera, color in zip("brz", "brk"):
        data = fitsio.read(f"{path}/thru-{camera}.fits", "THROUGHPUT")
        thru[camera] = np.interp(
            fullwave[cslice[camera]], data["wavelength"], data["throughput"]
        )
    return thru


# AR/DK exposure depths utilities
def load_spec(path):
    spec = {}
    with fitsio.FITS(str(path)) as hdus:
        for camera in "brz":
            spec[camera] = hdus[camera].read()
    return spec


# AR/DK settings for exposure depths
wmin, wmax, wdelta = 3600, 9824, 0.8
fullwave = np.round(np.arange(wmin, wmax + wdelta, wdelta), 1)
cslice = {"b": slice(0, 2751), "r": slice(2700, 5026), "z": slice(4900, 7781)}
spec_thru = load_spec_thru()
det_eso = load_spec(os.path.join(sv1dir, "misc", "dark_eso.fits"))
det_desimodel = load_spec(os.path.join(sv1dir, "misc", "dark_desimodel.fits"))
_sky_cache = {}
_depths = {}


# AR/DK exposure depths utilities
class Spectrum(object):
    def __init__(self, stype, flux=None, ivar=None, mask=None):
        assert stype == "full" or stype in cslice, "invalid stype"
        self.stype = stype
        self.wave = fullwave[cslice[stype]] if stype in cslice else fullwave
        if flux is None and ivar is None:
            self._flux = np.zeros(len(self.wave))
            self.ivar = np.zeros(len(self.wave))
        elif flux is not None and ivar is not None:
            self._flux = np.asarray(flux)
            self.ivar = np.asarray(ivar)
            assert (
                self.ivar.shape == self._flux.shape
            ), "flux and ivar have different shapes."
        else:
            raise ValueError("flux and ivar must both be specified.")
        if mask is None:
            self.mask = np.zeros_like(self._flux, bool)
        else:
            self.mask = np.asarray(mask)
            assert (
                self.mask.shape == self._flux.shape
            ), "flux and mask have different shapes."

    def copy(self):
        return Spectrum(
            self.stype, self.flux.copy(), self.ivar.copy(), self.mask.copy()
        )

    def __itruediv__(self, factor):
        np.divide(self.flux, factor, out=self._flux, where=factor != 0)
        self.ivar *= factor ** 2
        return self

    def __truediv__(self, factor):
        result = self.copy()
        result /= factor
        return result

    @property
    def flux(self):
        return self._flux


# AR/DK exposure depths utilities
class CoAdd(Spectrum):
    def __init__(self, stype):
        super(CoAdd, self).__init__(stype)
        self._weighted_flux_sum = np.zeros(len(self.wave))
        self._finalized = False

    def __iadd__(self, other):
        if other.stype == self.stype:
            self_slice = slice(None, None)
        elif self.stype == "full":
            self_slice = cslice[other.stype]
        else:
            raise ValueError(f'Cannot add "{other.stype}" to "{self.stype}".')
        self._weighted_flux_sum[self_slice] += other.ivar * other.flux
        self.ivar[self_slice] += other.ivar
        self._finalized = False
        return self

    @property
    def flux(self):
        if not self._finalized:
            np.divide(
                self._weighted_flux_sum, self.ivar, out=self._flux, where=self.ivar > 0
            )
            self._finalized = True
        return self._flux


# AR/DK exposure depths utilities
# AR/DK Estimate the average sky for a single exposure in phot/sec detected in each camera and incident on M1
def get_sky(night, expid, specs=range(10), use_cache=True, fill_cache=True):
    """
    Estimate the sky spectrum for one exposure in units of phot/sec per wavelength bin.
    Returns a tuple (flux_inc, ivar_inc, flux_det, ivar_det) where "det" is detected phot/sec
    in each camera with flat-field corrections applied, and "inc" corrects for the average
    spectrograph throughput in each camera, then coadds over cameras.
    """
    night = str(night)
    expid = str(expid).zfill(8)
    # print("running get_sky for {}-{}".format(night,expid))
    if use_cache and (night, expid) in _sky_cache:
        return _sky_cache[(night, expid)]
    incident = CoAdd("full")
    detected = {}
    # Loop over cameras.
    for camera in "brz":
        detected[camera] = CoAdd(camera)
        # Loop over spectrographs.
        for spec in specs:
            # Read the flat-fielded (constant) sky model in this spectrograph.
            skypath = os.path.join(
                dailydir,
                "exposures",
                night,
                expid,
                "sky-{}{}-{}.fits".format(camera, spec, expid),
            )
            if not os.path.isfile(skypath):
                print(f"Skipping non-existent {camera}{spec}.")
                continue
            with fitsio.FITS(str(skypath)) as hdus:
                exptime = hdus[0].read_header()["EXPTIME"]
                flux = hdus["SKY"].read()
                ivar = hdus["IVAR"].read()
                mask = hdus["MASK"].read()
                # Verify that we have the expected wavelengths.
                assert np.allclose(detected[camera].wave, hdus["WAVELENGTH"].read())
                # There are actually small variations in flux!
                # TODO: figure out where these variations come from.
                # For now, take the median over fibers.
                detected[camera] += Spectrum(
                    camera, np.median(flux, axis=0), np.median(ivar, axis=0)
                )
        # Scale to the exposure time.
        detected[camera] /= exptime
        # Correct for throughput and accumulate over cameras.
        incident += detected[camera] / spec_thru[camera]
    if fill_cache:
        _sky_cache[(night, expid)] = (incident, detected)
    return incident, detected


# AR/DK exposure depths utilities
def determine_tile_depth2(
    tileid,
    night,
    expid,
    exptime,
    transparency,
    fiber_fracflux,
    darkref=det_eso,
    ffracref=0.56,
    smoothing=125,
):
    tileid = str(tileid)
    night = str(night)
    expid = str(expid)
    exptimes = np.empty((6))
    exptimes[0] = exptime
    exptimes[1] = exptimes[0] * transparency ** 2
    exptimes[2] = exptimes[1] * (fiber_fracflux / ffracref) ** 2
    inc, det = get_sky(night, expid)
    for j, camera in enumerate("brz"):
        wave = det[camera].wave
        smoothref = scipy.ndimage.gaussian_filter1d(darkref[camera], smoothing)
        smooth = scipy.ndimage.gaussian_filter1d(det[camera].flux, smoothing)
        mean_ratio = np.sum(smooth) / np.sum(smoothref)
        exptimes[3 + j] = exptimes[2] / mean_ratio
    _depths[(tileid, night)] = exptimes
    return (
        tileid,
        night,
        expid,
        np.round(exptimes[3], 1),
        np.round(exptimes[4], 1),
        np.round(exptimes[5], 1),
    )


# AR mollweide plot setting
# AR http://balbuceosastropy.blogspot.com/2013/09/the-mollweide-projection.html
def set_mwd(ax, org=0):
    # org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder(tick_labels + 360 + org, 360)
    ax.set_xticklabels(tick_labels)  # we add the scale on the x axis
    ax.set_xlabel("R.A [deg]")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec. [deg]")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)
    return True


# AR mollweide coordinates conversion
def get_radec_mw(ra, dec, org):
    ra = np.remainder(ra + 360 - org, 360)  # shift ra values
    ra[ra > 180] -= 360  # scale conversion to [-180, 180]
    ra = -ra  # reverse the scale: East to the left
    return np.radians(ra), np.radians(dec)


# AR per-tile information
tiles = {}
tiles["FN"] = np.sort(glob(sv1dir + "202?????/fiberassign-??????.fits.gz"))
nt = len(tiles["FN"])
# AR initialising
for key in ["TILEID", "TILERA", "TILEDEC", "FAFLAVOR", "TARGETS", "COLOR", "FIELD"]:
    if key in ["TILEID"]:
        dtype = int
    elif key in ["FAFLAVOR", "TARGETS", "COLOR", "FIELD"]:
        dtype = object
    else:
        dtype = float
    tiles[key] = np.zeros(nt, dtype=dtype)
for target in targets:
    tiles[target] = np.zeros(nt, dtype=int)
# AR populating
for i in range(nt):
    # AR general
    hdr = fits.getheader(tiles["FN"][i])
    for key in ["TILEID", "TILERA", "TILEDEC", "FAFLAVOR"]:
        tiles[key][i] = hdr[key]
    # AR number of targets per tracer
    d = fits.open(tiles["FN"][i])[1].data
    if tiles["FAFLAVOR"][i] == "cmxm33":
        mask, key, msks, std_msks = cmx_mask, "cmx_target", cmx_msks, std_cmx_msks
    else:
        mask, key, msks, std_msks = desi_mask, "sv1_desi_target", sv1_msks, std_sv1_msks
    for target, msk in zip(targets, msks):
        if target in ["TGT", "SKY"]:
            tiles[target][i] = (d["objtype"] == target).sum()
        elif target == "STD":
            keep = np.zeros(len(d), dtype=bool)
            for std_msk in std_msks:
                keep |= (d[key] & mask[std_msk]) > 0
            tiles[target][i] = keep.sum()
        else:
            tiles[target][i] = ((d[key] & mask[msk]) > 0).sum()
# AR rounding coordinates to get a unique ra,dec for close tiles
prec = 1.0
tiles["radec"] = np.array(
    [
        "{:.1f},{:.1f}".format(prec * np.round(ra / prec), prec * np.round(dec / prec))
        for ra, dec in zip(tiles["TILERA"], tiles["TILEDEC"])
    ]
)
ref_faflavors = [
    ["cmxm33"],
    ["cmxlrgqso", "sv1lrgqso"],
    ["cmxelg", "sv1elg"],
    ["sv1bgsmws"],
]
ref_cols = ["magenta", "r", "b", "g"]
# AR extra infos
for key in list(flavdict.keys()):
    keep = np.in1d(tiles["FAFLAVOR"], flavdict[key]["FAFLAVORS"])
    tiles["TARGETS"][keep] = key
    tiles["COLOR"][keep] = flavdict[key]["COLOR"]
for key in list(fielddict.keys()):
    tiles["FIELD"][np.in1d(tiles["TILEID"], fielddict[key])] = key
# AR sorting by increasing TILEID
ii = np.argsort(tiles["TILEID"])
for key in list(tiles.keys()):
    tiles[key] = tiles[key][ii]


# AR gathering all TILEID-tiles.fits files in one
if args.tiles == "y":
    # AR each tile appears only once
    fns = [
        glob(os.path.join(sv1dir, "202?????", "{:06}-tiles.fits".format(tileid)))[0]
        for tileid in tiles["TILEID"]
    ]
    h = fits.open(fns[0])
    keys, fmts = h[1].columns.names, h[1].columns.formats
    t = {}
    for key in keys:
        t[key] = []
    for fn in fns:
        d = fits.open(fn)[1].data
        for key in keys:
            t[key] += [d[key]]
    # AR building/writing fits
    cols = []
    for key, fmt in zip(keys, fmts):
        cols += [fits.Column(name=key, format=fmt, array=t[key])]
    h = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
    h.writeto(outfns["tiles"], overwrite=True)


# AR per-exposure information (various from header, gfas + depths)
if args.exposures == "y":
    # AR listing existing nights and expids for the considered tiles
    # AR based on the presence of the sframe-??-EXPID.fits files
    nights = [
        int(fn.split("/")[-1])
        for fn in np.sort(glob(os.path.join(dailydir, "exposures", "202?????")))
        if int(fn.split("/")[-1]) >= int(firstnight)
    ]
    # AR GFA file
    gfa = fits.open(gfafn)[1].data
    keep = [program[:2] == "SV" for program in gfa["PROGRAM"]]
    gfa = gfa[keep]
    gfa_eci = np.array(
        ["{}-{}".format(e, c) for e, c in zip(gfa["EXPID"], gfa["CUBE_INDEX"])]
    )
    # AR quantities we store
    exposures = {}
    hdrkeys = ["NIGHT", "EXPID", "TILEID", "TILERA", "TILEDEC", "EXPTIME", "MJDOBS"]
    ownkeys = [
        "FIELD",
        "TARGETS",
        "EBV",
        "NGFA",
        "B_DEPTH",
        "R_DEPTH",
        "Z_DEPTH",
    ] + targets
    gfakeys = [
        "AIRMASS",
        "MOON_SEP_DEG",
        "TRANSPARENCY",
        "FWHM_ASEC",
        "SKY_MAG_AB",
        "FIBER_FRACFLUX",
    ]
    for key in hdrkeys + ownkeys:
        exposures[key] = []
    for key in gfakeys:
        exposures["{}_MIN".format(key)] = []
        exposures["{}_MED".format(key)] = []
        exposures["{}_MAX".format(key)] = []
    # AR looping on nights
    for night in nights:
        # AR first listing all exposures
        expids = np.unique(
            [
                fn.split("/")[-1]
                for fn in np.sort(
                    glob(
                        os.path.join(
                            dailydir, "exposures", "{}".format(night), "????????"
                        )
                    )
                )
            ]
        )
        # AR looping on all exposures
        for i in range(len(expids)):
            fns = glob(
                os.path.join(
                    dailydir,
                    "exposures",
                    "{}".format(night),
                    expids[i],
                    "sframe-??-{}.fits".format(expids[i]),
                )
            )
            if len(fns) > 0:
                hdr = fits.getheader(fns[0], 0)
                if hdr["TILEID"] in tiles["TILEID"]:
                    print(night, expids[i], hdr["TILEID"])
                    # AR header informations
                    for key in hdrkeys:
                        if key == "MJDOBS":
                            exposures[key] += [hdr["MJD-OBS"]]
                        else:
                            exposures[key] += [hdr[key]]
                    # AR field
                    it = np.where(tiles["TILEID"] == hdr["TILEID"])[0][0]
                    exposures["FIELD"] += [tiles["FIELD"][it]]
                    # AR targets
                    exposures["TARGETS"] += [tiles["TARGETS"][it]]
                    # AR ebv
                    exposures["EBV"] += [
                        float(
                            "{:.2f}".format(
                                np.median(
                                    fitsio.read(tiles["FN"][it], columns=["EBV"])["EBV"]
                                )
                            )
                        )
                    ]
                    # AR number of targets per tracer
                    for target in targets:
                        exposures[target] += [tiles[target][it]]
                    # AR GFA information
                    keep = gfa["SPECTRO_EXPID"] == hdr["EXPID"]
                    exposures["NGFA"] += [keep.sum()]
                    if keep.sum() > 0:
                        gfa_i = gfa[keep]
                        eci_i = gfa_eci[keep]
                        for key in gfakeys:
                            # AR first binning by expid-cube_index (see Aaron s email 09mar2020)
                            x = []
                            for eci in np.unique(eci_i):
                                tmp = eci_i == eci
                                x += [np.nanmedian(gfa_i[key][tmp])]
                            # AR then taking min,med,max
                            exposures["{}_MIN".format(key)] += [np.nanmin(x)]
                            exposures["{}_MED".format(key)] += [np.nanmedian(x)]
                            exposures["{}_MAX".format(key)] += [np.nanmax(x)]
                        # AR/DK exposure depths (needs gfa information)
                        for band, ib in zip(["B", "R", "Z"], [3, 4, 5]):
                            exposures["{}_DEPTH".format(band)] += [
                                determine_tile_depth2(
                                    exposures["TILEID"][-1],
                                    exposures["NIGHT"][-1],
                                    exposures["EXPID"][-1],
                                    exposures["EXPTIME"][-1],
                                    exposures["TRANSPARENCY_MED"][-1],
                                    exposures["FIBER_FRACFLUX_MED"][-1],
                                )[ib]
                            ]
                    else:
                        for key in gfakeys:
                            exposures["{}_MIN".format(key)] += [-99]
                            exposures["{}_MED".format(key)] += [-99]
                            exposures["{}_MAX".format(key)] += [-99]
                        for band in ["B", "R", "Z"]:
                            exposures["{}_DEPTH".format(band)] += [-99]

    # AR if update=y, pre-append the results from previous nights
    if args.update == "y":
        d = fits.open(outfns["exposures"])[1].data
        keep = d["NIGHT"] < int(firstnight)
        d = d[keep]
        for key in hdrkeys + ownkeys:
            exposures[key] = d[key].tolist() + exposures[key]
        for key in gfakeys:
            for suffix in ["_MIN", "_MED", "_MAX"]:
                exposures["{}{}".format(key, suffix)] = (
                    d["{}{}".format(key, suffix)].tolist()
                    + exposures["{}{}".format(key, suffix)]
                )
    # AR building/writing fits
    cols = []
    for key in hdrkeys + ownkeys:
        if key in ["NIGHT", "EXPID", "TILEID", "NGFA"] + targets:
            fmt = "K"
        elif key in ["FIELD", "TARGETS"]:
            fmt = "{}A".format(np.max([len(x) for x in exposures[key]]))
        else:
            fmt = "E"
        cols += [fits.Column(name=key, format=fmt, array=exposures[key])]
    for key in gfakeys:
        for suffix in ["_MIN", "_MED", "_MAX"]:
            cols += [
                fits.Column(
                    name="{}{}".format(key, suffix),
                    format="E",
                    array=exposures["{}{}".format(key, suffix)],
                )
            ]
    h = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
    h.writeto(outfns["exposures"], overwrite=True)


# AR tables to be copied-pasted in the wiki
if args.wiki == "y":
    f = open(outfns["wiki"], "w")

    # AR tile design (fiberassign, log, QA plot, viewer, split per tracer)
    d = fits.open(outfns["exposures"])[1].data
    _, ii = np.unique(d["TILEID"], return_index=True)
    ii = ii[d["TILEID"][ii].argsort()]
    d = d[ii]
    f.write("=================== TILE DESIGN ==================\n")
    f.write("\n")
    fields = [
        "TILEID",
        "Name",
        "Targets",
        "RA",
        "Dec",
        "Fits",
        "QA plot",
        "Log",
        "Viewer",
    ] + targets
    f.write("||= **{}** =||\n".format(" =||=".join(fields)))
    for i in range(len(d)):
        tmparr = ["{:06}".format(d["TILEID"][i])]
        tmparr += [d["FIELD"][i]]
        tmparr += [d["TARGETS"][i]]
        tmparr += ["{:.3f}".format(d["TILERA"][i])]
        tmparr += ["{:.3f}".format(d["TILEDEC"][i])]
        tmparr += [
            "[https://desi.lbl.gov/svn/data/tiles/trunk/{}/fiberassign-{:06}.fits.gz fiberassign-{:06}.fits.gz]".format(
                str(d["TILEID"][i]).zfill(6)[:3], d["TILEID"][i], d["TILEID"][i]
            )
        ]
        tmparr += [
            "[https://desi.lbl.gov/svn/data/tiles/trunk/{}/fiberassign-{:06}.png fiberassign-{:06}.png]".format(
                str(d["TILEID"][i]).zfill(6)[:3], d["TILEID"][i], d["TILEID"][i]
            )
        ]
        tmparr += [
            "[https://desi.lbl.gov/svn/data/tiles/trunk/{}/{:06}.log {:06}.log]".format(
                str(d["TILEID"][i]).zfill(6)[:3], d["TILEID"][i], d["TILEID"][i]
            )
        ]
        tmparr += [
            "[https://www.legacysurvey.org/viewer-dev/?ra={:.3f}&dec={:.3f}&layer=ls-dr9&zoom=9 Viewer]".format(
                d["TILERA"][i], d["TILEDEC"][i]
            )
        ]
        tmparr += ["{}".format(d[target][i]) for target in targets]
        f.write("||{} ||\n".format(" ||".join(tmparr)))
    f.write("\n")
    f.write("\n")
    f.write("\n")

    # AR observed exposures
    d = fits.open(outfns["exposures"])[1].data
    tileids = np.sort(np.unique(d["TILEID"]))
    f.write("=================== NB OF EXPOSURES  =============\n")
    f.write("\n")
    fields = ["TILEID", "Name", "Targets"] + ["Total nb exp.", "Nb. exp. per night"]
    f.write("||= **{}** =||\n".format(" =||=".join(fields)))
    for tileid in tileids:
        ii = d["TILEID"] == tileid
        di = d[ii]
        di = di[di["EXPID"].argsort()]
        tmparr = ["{:06}".format(di["TILEID"][0])]
        tmparr += [di["FIELD"][0]]
        tmparr += [di["TARGETS"][0]]
        # total nb exp
        tmparr += ["{}".format(len(di))]
        # 1st night
        nights = np.unique(di["NIGHT"])
        j = 0
        jj = di["NIGHT"] == nights[j]
        texps = di["EXPTIME"][jj].astype(int)
        ts, cs = np.unique(texps, return_counts=True)
        tmparr += [
            "{}:{}".format(
                nights[j], ",".join(["{}x{}s".format(c, t) for c, t in zip(cs, ts)])
            )
        ]
        f.write("||{} ||\n".format(" ||".join(tmparr)))
        # next nights, if any
        if len(nights) > 1:
            for j in range(1, len(nights)):
                tmparr = ["" for k in range(len(fields) - 1)]
                jj = di["NIGHT"] == nights[j]
                texps = di["EXPTIME"][jj].astype(int)
                ts, cs = np.unique(texps, return_counts=True)
                tmparr += [
                    "{}:{}".format(
                        nights[j],
                        ",".join(["{}x{}s".format(c, t) for c, t in zip(cs, ts)]),
                    )
                ]
                f.write("||{} ||\n".format(" ||".join(tmparr)))
    f.write("\n")
    f.write("\n")
    f.write("\n")

    # AR observing conditions
    d = fits.open(outfns["exposures"])[1].data
    tileids = np.sort(np.unique(d["TILEID"]))
    f.write("=================== OBSERVING CONDITIONS  =============\n")
    f.write("\n")
    keys = [
        "AIRMASS",
        "MOON_SEP_DEG",
        "TRANSPARENCY",
        "FWHM_ASEC",
        "SKY_MAG_AB",
        "FIBER_FRACFLUX",
    ]
    fields = (
        ["TILEID", "NIGHT", "EXPID", "EXPTIME", "EBV"]
        + keys
        + ["B_DEPTH", "R_DEPTH", "Z_DEPTH"]
    )
    for tileid in tileids:
        f.write("||= **{}** =||\n".format(" =||=".join(fields)))
        ii = np.where(d["TILEID"] == tileid)[0]
        ii = ii[d["EXPID"][ii].argsort()]
        di = d[ii]
        for j in range(len(di)):
            tmparr = ["{:06}".format(tileid)]
            tmparr += ["{}".format(di["NIGHT"][j])]
            tmparr += ["{}".format(di["EXPID"][j])]
            tmparr += ["{:.0f}".format(di["EXPTIME"][j])]
            tmparr += ["{:.2f}".format(di["EBV"][j])]
            tmparr += ["{:.2f}".format(di[key + "_MED"][j]) for key in keys]
            tmparr += [
                "{:.0f}s".format(di[band + "_DEPTH"][j]) for band in ["B", "R", "Z"]
            ]
            f.write("||{} ||\n".format(" ||".join(tmparr)))
    f.close()


# AR sky map
if args.skymap == "y":
    # dr9
    h = fits.open(pixwfn)
    nside, nest = h[1].header["HPXNSIDE"], h[1].header["HPXNEST"]
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix, dtype=int), nest=nest)
    hpdict = {}
    hpdict["ra"], hpdict["dec"] = 180.0 / np.pi * phi, 90.0 - 180.0 / np.pi * theta
    for key in [
        "fracarea",
        "stardens",
        "ebv",
        "psfsize_g",
        "psfsize_r",
        "psfsize_z",
        "galdepth_g",
        "galdepth_r",
        "galdepth_z",
        "psfdepth_w1",
        "psfdepth_w2",
    ]:
        if key == "stardens":
            hpdict[key] = np.log10(h[1].data[key])
            hpdict[key + "lab"] = "log10(stardens)"
        elif (key[:8] == "galdepth") | (key[:8] == "psfdepth"):
            hpdict[key] = 22.5 - 2.5 * np.log10(5.0 / np.sqrt(h[1].data[key]))
            hpdict[key + "lab"] = r"5$\sigma$ " + key
        else:
            hpdict[key] = h[1].data[key]
            hpdict[key + "lab"] = key
        if key[:7] == "psfsize":
            hpdict[key][hpdict[key] == 0] = np.nan
    ## north/south
    c = SkyCoord(
        hpdict["ra"] * units.degree, hpdict["dec"] * units.degree, frame="icrs"
    )
    hpdict["north"] = (
        (hpdict["fracarea"] > 0) & (hpdict["dec"] > 32.375) & (c.galactic.b.value > 0)
    )
    hpdict["south"] = (hpdict["fracarea"] > 0) & (~hpdict["north"])

    # plotting skymap
    projection = "mollweide"
    org = 120  # centre ra for mollweide plots
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111, projection=projection)
    _ = set_mwd(ax, org)
    # dr9
    ramw, decmw = get_radec_mw(hpdict["ra"], hpdict["dec"], org)
    for s, a in zip([hpdict["north"], hpdict["south"]], [0.05, 0.6]):
        ax.scatter(ramw[s], decmw[s], s=1, c="0.8", zorder=0, alpha=a, rasterized=True)
    # des
    ra, dec = np.loadtxt(desfn, unpack=True)
    ramw, decmw = get_radec_mw(ra, dec, org)
    ax.plot(ramw, decmw, c="k", lw=0.5)
    #
    ramw, decmw = get_radec_mw(tiles["TILERA"], tiles["TILEDEC"], org)
    for radec in np.unique(tiles["radec"]):
        ii = np.where(tiles["radec"] == radec)[0]
        n = len(ii)
        dy = 5
        y = np.array([float(radec.split(",")[1]) + (n - 1) / 2.0 * dy])
        for i in ii:
            ax.scatter(ramw[i], decmw[i], c=tiles["COLOR"][i], marker="X", s=50)
            if tiles["TILEID"][i] in [80611, 80614, 80617]:
                dx, ha = +3, "right"
            else:
                dx, ha = -3, "left"
            x = np.array([float(radec.split(",")[0]) + dx])
            mwx, mwy = get_radec_mw(x, y, org)
            ax.text(
                mwx,
                mwy,
                "{}".format(tiles["TILEID"][i]),
                color=tiles["COLOR"][i],
                ha=ha,
                va="center",
            )
            y -= dy

    for faflavors, col in zip(ref_faflavors, ref_cols):
        faflavors = np.unique(
            [faflavor.replace("cmx", "").replace("sv1", "") for faflavor in faflavors]
        )
        ax.scatter(100, 100, marker="X", s=50, c=col, label=",".join(faflavors))
    ax.legend(loc=2)
    plt.savefig(outfns["skymap"], bbox_inches="tight")
    plt.close()


# AR r_depth plot
if args.depth == "y":
    d = fits.open(outfns["exposures"])[1].data
    xlim = (0, 30)
    xs = 0.5 + np.arange(xlim[0], xlim[1])
    # AR color-coding by night
    ref_cols = ["r", "g", "b"]
    nights = np.unique(d["NIGHT"])
    cols = np.zeros(len(d), dtype=object)
    for i in range(len(nights)):
        keep = d["NIGHT"] == nights[i]
        cols[keep] = ref_cols[i % len(ref_cols)]
    key = "R_DEPTH"
    for flavshort, ymax in zip(["QSO+LRG", "ELG", "BGS+MWS"], [2000, 2000, 1000]):
        keep = d["TARGETS"] == flavshort
        tileids = np.unique(d["TILEID"][keep])
        nightmin, nightmax = d["NIGHT"][keep].min(), d["NIGHT"][keep].max()
        fig = plt.figure(figsize=(25, 1 * len(tileids)))
        gs = gridspec.GridSpec(len(tileids), 1, hspace=0)
        for i in range(len(tileids)):
            ax = plt.subplot(gs[i])
            ax.text(
                0.98,
                0.80,
                tileids[i],
                color="k",
                fontweight="bold",
                ha="right",
                transform=ax.transAxes,
            )
            ax.text(
                0.98,
                0.55,
                d["FIELD"][d["TILEID"] == tileids[i]][0],
                color="k",
                fontweight="bold",
                ha="right",
                transform=ax.transAxes,
            )
            jj = np.where(d["TILEID"] == tileids[i])[0]
            ax.plot(xs[: len(jj)], d[key][jj], color="k", lw=1)
            x = 0
            for j in jj:
                ax.text(
                    0.5 + x,
                    0.75 * ymax,
                    "{}\n{}\n{:.0f}s".format(d["NIGHT"][j], d["EXPID"][j], d[key][j]),
                    ha="center",
                    va="center",
                    fontsize=7,
                    color=cols[j],
                )
                ax.scatter(0.5 + x, d[key][j], c=cols[j], marker="o", s=5)
                ax.plot(
                    [x, x + 1], d["EXPTIME"][j] + np.zeros(2), c="k", ls="--", lw=0.5
                )
                x += 1
            ax.grid(True)
            ax.set_axisbelow(True)
            ax.set_xlim(0, 30)
            if i == int(len(tileids) / 2):
                ax.set_ylabel("{} [s]".format(key))
            for x in range(xlim[0], xlim[1]):
                ax.axvline(x, c="k", lw=0.1)
            if i == 0:
                ax.set_title(
                    "SV1 {} ({} exposures from {} tiles between {} and {})".format(
                        flavshort, keep.sum(), len(tileids), nightmin, nightmax
                    )
                )
            if i == len(tileids) - 1:
                ax.set_xlabel("Exposure #")
            else:
                ax.set_xticks([])
            ax.set_ylim(0.01, ymax - 0.01)
            if flavshort == "BGS+MWS":
                ax.yaxis.set_major_locator(MultipleLocator(250))
            else:
                ax.yaxis.set_major_locator(MultipleLocator(500))
        plt.savefig(outfns["depth"][flavshort], bbox_inches="tight")
        plt.close()
