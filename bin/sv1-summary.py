#!/usr/bin/env python

import matplotlib

matplotlib.use("Agg")
import sys
import os
from glob import glob
import numpy as np
import fitsio
import astropy.io.fits as fits
from astropy.time import Time
from astropy.table import Table
from astropy import units, constants
from astropy.coordinates import SkyCoord
import healpy as hp
from desitarget.sv1.sv1_targetmask import desi_mask
from desitarget.cmx.cmx_targetmask import cmx_mask
import fitsio
import scipy.ndimage
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
import textwrap
from speclite import filters
import yaml
from desispec.io import read_frame, fiberflat
from desispec.fiberflat import apply_fiberflat
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
    "--plot",
    help="do plots (skymap, observing conditions, r_depth)? (y/n)",
    type=str,
    default="y",
    metavar="PLOT",
)
parser.add_argument(
    "--html",
    help="create html pages, one per tile (y/n)",
    type=str,
    default="y",
    metavar="HTML",
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
sv1dir = os.getenv("DESI_ROOT")+"/users/raichoor/fiberassign-sv1/"
dailydir = os.getenv("DESI_ROOT")+"/spectro/redux/daily/"
pixwfn = os.getenv("DESI_TARGET")+"/catalogs/dr9/0.47.0/pixweight/sv1/resolve/dark/sv1pixweight-dark.fits"
desfn = os.path.join(sv1dir, "misc", "des_footprint.txt")
gfafn = np.sort(
    glob(
        os.getenv("DESI_ROOT")+"/users/ameisner/GFA/conditions/offline_all_guide_ccds_SV1-thru_20??????.fits"
    )
)[-1]


# AR for the fiberflat routine
os.environ["DESI_LOGLEVEL"] = "warning"


# AR output products
outfns = {}
outfns["tiles"] = os.path.join(args.outdir, "sv1-tiles.fits")
outfns["exposures"] = os.path.join(args.outdir, "sv1-exposures.fits")
outfns["wiki"] = os.path.join(args.outdir, "sv1-tables.wiki")
if os.path.isdir(os.path.join(args.outdir, "sv1-plots")) == False:
    os.mkdir(os.path.join(args.outdir, "sv1-plots"))
outfns["skymap"] = os.path.join(args.outdir, "sv1-plots", "sv1-skymap.png")
months = ["202012", "202101", "202102", "202103", "202104"]  # AR for splitting plots
outfns["obscond"] = {}
for month in months:
    outfns["obscond"][month] = os.path.join(
        args.outdir, "sv1-plots", "sv1-obscond-{}.png".format(month)
    )
outfns["depth"] = {}
for flavshort in ["QSO+LRG", "ELG", "BGS+MWS", "QSO+ELG"]:
    outfns["depth"][flavshort] = os.path.join(
        args.outdir,
        "sv1-plots",
        "sv1-depth-{}.png".format(flavshort.lower().replace("+", "")),
    )
outfns["html"] = os.path.join(args.outdir, "sv1-per-tile")
if (args.html == "y") & (os.path.isdir(outfns["html"]) == False):
    os.mkdir(outfns["html"])

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
    "QSO+ELG": {"FAFLAVORS": ["sv1elgqso"], "COLOR": "c"},
    "BGS+MWS": {"FAFLAVORS": ["cmxbgsmws", "sv1bgsmws"], "COLOR": "g"},
    "M33+Dark": {"FAFLAVORS": ["cmxm33"], "COLOR": "magenta"},
    "M31": {"FAFLAVORS": ["sv1m31"], "COLOR": "y"},
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
    # BGS+MWS 20210101 , Dark 20210105
    "Dust cloud": [80624],
    "Far north": [80650, 80655, 80693, 80694],
    "Far south": [80630, 80638, 80671, 80701, 80672, 80702],
    "GAMA G02": [80633, 80635],
    "GAMA G12": [80661, 80662, 80663, 80705, 80706],
    "G.Plane b=+15": [80642],
    "G.Plane b=+16": [80641, 80675, 80676],
    "G.Plane b=+18": [80644, 80683, 80684],
    "G.Plane b=-24": [80640, 80673, 80674],
    "G.Plane b=-27": [80639],
    "HSC": [80632, 80634],
    "HSC/S82": [80625, 80637, 80669, 80670],
    "MWS south edge of NGC": [80647],
    "MWS GD1_LOW_2 blk": [80648],
    "MWS MONOC._LOW blk": [80645],
    "MWS Sag. stream blk hp": [80628, 80629],
    "MWS TRIAND_STR._1 blk": [80626],
    "MWS TRIAND_STR._2 blk": [80627],
    "Moon avoidance": [80631, 80636, 80667, 80668],
    "Sag. stream": [80665, 80666, 80709, 80710],
    "MWS GD1-C-1 blk": [80658],
    "MWS North G.Pole blk": [80664],
    "MWS ORPHAN-A-1 blk": [80657],
    "Overlap": [
        80646,
        80649,
        80651,
        80654,
        80659,
        80689,
        80691,
        80697,
        80699,
        80703,
        80690,
        80692,
        80698,
        80700,
        80704,
    ],
    "Overlap and b=+23": [80643],
    "MWS GD-B-1 str. blk hp": [80652],
    "MWS GD-B-2 str. blk hp": [80653],
    "MWS GD-B-3 str. blk hp": [80656],
    "MWS GD-C-2 str. blk hp": [80660],
    "Beehive cluster": [80687, 80688],
    "Gal. extinction": [80695, 80696],
    "Monoc. stream": [80685, 80686],
    "N. Gal. Pole": [80707, 80708],
    "2-pass g-depth": [80679, 80680],
    "2-pass r-depth + G.Plane b=+17": [80677, 80678],
    "DEEP2 EGS": [80711, 80712],
    "Overlap+Monoc. stream": [80681, 80682],
    "M31cen": [80713, 80714, 80715, 80716],
}


# AR/DK DESI spectra wavelengths
wmin, wmax, wdelta = 3600, 9824, 0.8
fullwave = np.round(np.arange(wmin, wmax + wdelta, wdelta), 1)
cslice = {"b": slice(0, 2751), "r": slice(2700, 5026), "z": slice(4900, 7781)}


# AR DESI telescope geometric area (cm2) and fiber area (arcsec2)
# AR for computing SKY_RMAG_AB
fn = os.path.join(os.getenv("DESIMODEL"), "data", "desi.yaml")
f = open(fn, "r")
desi = yaml.safe_load(f)
f.close()
telap_cm2 = desi["area"]["geometric_area"] * 1e4  # AR telescope geometric area in cm2
fiber_area_arcsec2 = (
    np.pi * (desi["fibers"]["diameter_arcsec"] / 2.0) ** 2
)  # fiber area in arcsec2


# AR/DK exposure depths utilities
def load_spec_thru(path=os.getenv("DESIMODEL") + "/data/throughput/"):
    thru = {}
    for camera in ["b", "r", "z"]:
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
spec_thru = load_spec_thru()
det_eso = load_spec(os.path.join(sv1dir, "misc", "dark_eso.fits"))
det_desimodel = load_spec(os.path.join(sv1dir, "misc", "dark_desimodel.fits"))
_sky_cache = {}


# AR/ES ebv and airmass coefficient in depth_ebvair
# AR/ES ebv coeffs : taking SDSS grz from Schlafly & Finkbeinger (2011)
# AR/ES airmass coeffs : [decam-chatter 15497]
depth_coeffs = {"EBV": {"B":3.303, "R":2.285, "Z":1.263}, "AIRMASS": {"B":0.195, "R":0.096, "Z":0.055}}


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
    # print("running get_sky for {}-{}".format(night,expid))
    if use_cache and (night, expid) in _sky_cache:
        return _sky_cache[(night, expid)]
    incident = CoAdd("full")
    detected = {}
    # Loop over cameras.
    for camera in ["b", "r", "z"]:
        detected[camera] = CoAdd(camera)
        # Loop over spectrographs.
        for spec in specs:
            # Read the flat-fielded (constant) sky model in this spectrograph.
            skypath = os.path.join(
                dailydir,
                "exposures",
                "{}".format(night),
                "{:08}".format(expid),
                "sky-{}{}-{:08}.fits".format(camera, spec, expid),
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
# AR transpfrac = transparency * fiber_fracflux
def determine_tile_depth2(
    night, expid, exptime, transpfrac, darkref=det_eso, ffracref=0.56, smoothing=125,
):
    inc, det = get_sky(night, expid)
    depths = {}
    for camera in ["b", "r", "z"]:
        wave = det[camera].wave
        smoothref = scipy.ndimage.gaussian_filter1d(darkref[camera], smoothing)
        smooth = scipy.ndimage.gaussian_filter1d(det[camera].flux, smoothing)
        mean_ratio = np.sum(smooth) / np.sum(smoothref)
        depths[camera] = np.round(
            exptime * (transpfrac / (1.0 * ffracref)) ** 2 / mean_ratio, 1
        )
    return depths


# AR r-band sky mag / arcsec2 from sky-....fits files
def get_sky_rmag_ab(night, expid, exptime, ftype, redux="daily"):
    # AR ftype = "data" or "model"
    # AR redux = "daily" or "blanc"
    # AR if ftype = "data" : read the sky fibers from frame*fits + apply flat-field
    # AR if ftype = "model": read the sky model from sky*.fits for the first fiber of each petal (see DJS email from 29Dec2020)
    # AR those are in electron / angstrom
    # AR to integrate over the decam-r-band, we need cameras b and r
    if ftype not in ["data", "model"]:
        sys.exit("ftype should be 'data' or 'model'")
    sky = np.zeros(len(fullwave))
    reduxdir = dailydir.replace("daily", redux)
    # AR see AK email [desi-data 5218]
    if redux == "blanc":
        specs = ["0", "1", "3", "4", "5", "7", "8", "9"]
    else:
        specs = np.arange(10, dtype=int).astype(str)
    for camera in ["b", "r"]:
        norm_cam = np.zeros(len(fullwave[cslice[camera]]))
        sky_cam = np.zeros(len(fullwave[cslice[camera]]))
        for spec in specs:
            # AR data
            if ftype == "data":
                frfn = os.path.join(
                    reduxdir,
                    "exposures",
                    "{}".format(night),
                    expid,
                    "frame-{}{}-{}.fits".format(camera, spec, expid),
                )
                flfn = os.path.join(
                    reduxdir,
                    "calibnight",
                    "{}".format(night),
                    "fiberflatnight-{}{}-{}.fits".format(camera, spec, night),
                )
                if not os.path.isfile(frfn) or not os.path.isfile(flfn):
                    print("Skipping non-existent {}, {}".format(frfn, flfn))
                else:
                    fr = read_frame(frfn, skip_resolution=True)
                    fl = fiberflat.read_fiberflat(flfn)
                    apply_fiberflat(fr, fl)
                    # AR cutting on sky fibers with at least one valid pixel
                    ii = (fr.fibermap["OBJTYPE"] == "SKY") & (fr.ivar.sum(axis=1) > 0)
                    # AR frame*fits are in e- / angstrom ; adding the N sky fibers
                    # sky_cam += fr.flux[ii, :].sum(axis=0)
                    # nspec += ii.sum()
                    sky_cam += (fr.flux[ii, :] * fr.ivar[ii, :]).sum(axis=0)
                    norm_cam += fr.ivar[ii, :].sum(axis=0)
            # AR model
            if ftype == "model":
                fn = os.path.join(
                    reduxdir,
                    "exposures",
                    "{}".format(night),
                    expid,
                    "sky-{}{}-{}.fits".format(camera, spec, expid),
                )
                if not os.path.isfile(fn):
                    print("Skipping non-existent {}".format(fn))
                else:
                    fd = fitsio.FITS(fn)
                    assert np.allclose(
                        fullwave[cslice[camera]], fd["WAVELENGTH"].read()
                    )
                    fd = fitsio.FITS(fn)
                    # AR sky*fits are in e- / angstrom
                    # AR handling some cases where SKY=IVAR=0
                    if fd["IVAR"][0, :][0].max() > 0:
                        sky_cam += fd["SKY"][0, :][0]  # AR reading the first fiber only
                        # nspec += 1
                        norm_cam += np.ones(len(fullwave[cslice[camera]]))
                    else:
                        print(
                            "{}-{}-{}{}: no spectra for {}".format(
                                night, expid, camera, spec, ftype
                            )
                        )
                    fd.close()
        # AR sky model flux in incident photon / angstrom / s
        # if nspec > 0:
        keep = norm_cam > 0
        if keep.sum() > 0:
            sky[cslice[camera]][keep] = (
                sky_cam[keep] / norm_cam[keep] / exptime / spec_thru[camera][keep]
            )
        else:
            print("{}-{}-{}: no spectra for {}".format(night, expid, camera, ftype))
    # AR sky model flux in erg / angstrom / s (using the photon energy in erg)
    e_phot_erg = (
        constants.h.to(units.erg * units.s)
        * constants.c.to(units.angstrom / units.s)
        / (fullwave * units.angstrom)
    )
    sky *= e_phot_erg.value
    # AR sky model flux in erg / angstrom / s / cm**2 / arcsec**2
    sky /= telap_cm2 * fiber_area_arcsec2
    # AR integrate over the DECam r-band
    filts = filters.load_filters("decam2014-r")
    # AR zero-padding spectrum so that it covers the DECam r-band range
    sky_pad, fullwave_pad = filts.pad_spectrum(sky, fullwave, method="zero")
    return filts.get_ab_magnitudes(
        sky_pad * units.erg / (units.cm ** 2 * units.s * units.angstrom),
        fullwave_pad * units.angstrom,
    ).as_array()[0][0]


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


# AR/ADM from desitarget/QA.py
def _javastring():
    """Return a string that embeds a date in a webpage
    """

    js = textwrap.dedent(
        """
    <SCRIPT LANGUAGE="JavaScript">
    var months = new Array(13);
    months[1] = "January";
    months[2] = "February";
    months[3] = "March";
    months[4] = "April";
    months[5] = "May";
    months[6] = "June";
    months[7] = "July";
    months[8] = "August";
    months[9] = "September";
    months[10] = "October";
    months[11] = "November";
    months[12] = "December";
    var dateObj = new Date(document.lastModified)
    var lmonth = months[dateObj.getMonth() + 1]
    var date = dateObj.getDate()
    var fyear = dateObj.getYear()
    if (fyear < 2000)
    fyear = fyear + 1900
    if (date == 1 || date == 21 || date == 31)
    document.write(" " + lmonth + " " + date + "st, " + fyear)
    else if (date == 2 || date == 22)
    document.write(" " + lmonth + " " + date + "nd, " + fyear)
    else if (date == 3 || date == 23)
    document.write(" " + lmonth + " " + date + "rd, " + fyear)
    else
    document.write(" " + lmonth + " " + date + "th, " + fyear)
    </SCRIPT>
    """
    )
    return js


# AR html tile design
def write_html_tiledesign(html, tiles, ii, nexps, style, h2title, main=True):
    fields = [
        "TILEID",
        "NEXP",
        "Name",
        "Targets",
        "RA",
        "Dec",
        "Fiber Assign Fits",
        "Fiber Assign QA plot",
        "Fiber Assign Log",
        "Viewer",
    ] + targets
    html.write(
        "<h2><a id='tile-nexp-design' href='#tile-nexp-design' > {}</a>\n".format(
            h2title
        )
    )
    html.write(
        "<a href='#top' style='position: absolute; right: 0;'>Top of the page</a></h2>\n"
    )
    if main:
        html.write(
            "<p style='{}'>Click on the TILEID to access its detailed SV1 html page.</p>".format(
                style
            )
        )
    html.write(
        "<p style='{}'>The exact number of tracers with FIBERSTATUS == 0 fluctuates from exposure to exposure, but is overall 80%.</p>".format(
            style
        )
    )
    html.write("<table>\n")
    html.write("<style>\n")
    html.write("th, td {border:1px solid black; font-size: 0.95em}\n")
    html.write("tr:nth-child(even) {background-color: " + bkgcol + ";}\n")
    html.write("</style>\n")
    count = 0
    for i, nexp in zip(ii, nexps):
        if count % 20 == 0:
            html.write("<tr>\n")
            html.write(" ".join(["<th> {} </th>".format(x) for x in fields]) + "\n")
            html.write("</tr>\n")
        count += 1
        fafits = "https://desi.lbl.gov/svn/data/tiles/trunk/{}/fiberassign-{:06}.fits.gz".format(
            str(tiles["TILEID"][i]).zfill(6)[:3], tiles["TILEID"][i]
        )
        if main:
            tmparr = [
                "<a href='{:06}.html'>{:06}</a>".format(
                    tiles["TILEID"][i], tiles["TILEID"][i]
                )
            ]
        else:
            tmparr = ["{:06}".format(tiles["TILEID"][i])]
        tmparr += ["{}".format(nexp)]
        tmparr += [tiles["FIELD"][i]]
        tmparr += [tiles["TARGETS"][i]]
        tmparr += ["{:.3f}".format(tiles["TILERA"][i])]
        tmparr += ["{:.3f}".format(tiles["TILEDEC"][i])]
        tmparr += [
            "<a href='{}' target='external'> fiberassign-{:06}.fits.gz </a>".format(
                fafits, tiles["TILEID"][i]
            )
        ]
        tmparr += [
            "<a href='{}' target='external'> fiberassign-{:06}.png".format(
                fafits.replace(".fits.gz", ".png"), tiles["TILEID"][i]
            )
        ]
        tmparr += [
            "<a href='{}' target='external'> {:06}.log".format(
                os.path.join(
                    "/".join(fafits.split("/")[:-1]),
                    "{:06}.log".format(tiles["TILEID"][i]),
                ),
                tiles["TILEID"][i],
            )
        ]
        tmparr += [
            "<a href='https://www.legacysurvey.org/viewer-dev/?ra={:.3f}&dec={:.3f}&layer=ls-dr9&zoom=8&tile={}' target='external'> Viewer".format(
                tiles["TILERA"][i], tiles["TILEDEC"][i], tiles["TILEID"][i]
            )
        ]
        tmparr += ["{}".format(tiles[target][i]) for target in targets]
        html.write("<tr>\n")
        html.write(" ".join(["<td> {} </td>".format(x) for x in tmparr]) + "\n")
        html.write("</tr>" + "\n")
    html.write("</table>\n")
    html.write("\n")
    return True


# AR html per-exposure table
def write_html_perexp(html, d, style, h2title):
    html.write(
        "<h2><a id='per-exposure-properties' href='#per-exposure-properties' > {} </a> ({} exposure(s) over {} night(s))".format(
            h2title, len(d), len(np.unique(d["NIGHT"]))
        )
    )
    html.write(
        "<a href='#top' style='position: absolute; right: 0;'>Top of the page</a></h2>\n"
    )
    html.write(
        "<p style='{}'>Click on the TILEID / NIGHT / EXPID / NIGHTWATCH to access its folder with processed files.</p>".format(
            style
        )
    )
    html.write("<p style='{}'>GFA fits file: {}</p>".format(style, gfafn))
    html.write(
        "<p style='{}'>FIBER_FRACFLUX: fraction of light in a fiber-sized aperture given the PSF shape, assuming that the PSF is perfectly aligned with the fiber (i.e. does not capture any astrometry/positioning errors).</p>".format(
            style
        )
    )
    html.write(
        "<p style='{}'>R_DEPTH: EXPTIME x (TRANSPARENCY x FIBER_FRACFLUX / (1.0 x 0.56))^2 x FIDSKY_DARK/EXPSKY (does not correct for EBV, AIRMASS).</p>".format(
            style
        )
    )
    html.write(
        "<p style='{}'>R_DEPTH_EBVAIR: R_DEPTH x 10^(-2 x 0.4 x 2.285 x EBV) x 10^(-2 x 0.4 x 0.096 x AIRMASS).</p>".format(
            style
        )
    )
    # ADM write out a list of the target categories.
    keys = [
        "AIRMASS",
        "MOON_SEP_DEG",
        "TRANSPARENCY",
        "FWHM_ASEC",
        "SKY_MAG_AB",
        "FIBER_FRACFLUX",
    ]
    fields = (
        ["TILEID", "NIGHT", "EXPID", "NIGHTWATCH", "EXPTIME", "EBV"]
        + keys
        #+ ["B_DEPTH", "R_DEPTH", "Z_DEPTH"]
        + ["R_DEPTH", "R_DEPTH_EBVAIR"]
    )
    html.write("<table>\n")
    night = ""
    for i in range(len(d))[::-1]:
        night_prev = night
        night = d["NIGHT"][i]
        if night != night_prev:
            html.write("<tr>\n")
            html.write(" ".join(["<th> {} </th>".format(x) for x in fields]) + "\n")
            html.write("</tr>\n")
        html.write("<tr>")
        tmparr = [
            "<a href='{}' target='external'> {}".format(
                os.path.join(
                    dailydir.replace(
                        os.getenv("DESI_ROOT"), "https://data.desi.lbl.gov/desi"
                    ),
                    "tiles",
                    "{}".format(d["TILEID"][i]),
                ),
                d["TILEID"][i],
            )
        ]
        tmparr += [
            "<a href='{}' target='external'> {}".format(
                os.path.join(
                    dailydir.replace(
                        os.getenv("DESI_ROOT"), "https://data.desi.lbl.gov/desi"
                    ),
                    "exposures",
                    "{}".format(d["NIGHT"][i]),
                ),
                d["NIGHT"][i],
            )
        ]
        tmparr += [
            "<a href='{}' target='external'> {}".format(
                os.path.join(
                    dailydir.replace(
                        os.getenv("DESI_ROOT"), "https://data.desi.lbl.gov/desi"
                    ),
                    "exposures",
                    "{}".format(d["NIGHT"][i]),
                    "{:08}".format(d["EXPID"][i]),
                ),
                d["EXPID"][i],
            )
        ]
        tmparr += [
            "<a href='{}' target='external'> {}".format(
                os.path.join(
                    "https://nightwatch.desi.lbl.gov/" "{}".format(d["NIGHT"][i]),
                    "{:08}".format(d["EXPID"][i]),
                    "qa-summary-{:08}.html".format(d["EXPID"][i]),
                ),
                "Nightwatch",
            )
        ]
        tmparr += ["{:.0f}".format(d["EXPTIME"][i])]
        tmparr += ["{:.2f}".format(d["EBV"][i])]
        tmparr += ["{:.2f}".format(d["GFA_" + key + "_MED"][i]) for key in keys]
        #tmparr += ["{:.0f}s".format(d[band + "_DEPTH"][i]) for band in ["B", "R", "Z"]]
        tmparr += ["{:.0f}s".format(d["R_DEPTH"][i])]   
        tmparr += ["{:.0f}s".format(d["R_DEPTH_EBVAIR"][i])]
        html.write(" ".join(["<td> {} </td>".format(x) for x in tmparr]) + "\n")
        html.write("</tr>\n")
    html.write("</table>\n")
    html.write("\n")
    return True


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
    # AR see Aaron s email from 01Jan2021: using ext=2, which already contains
    # AR the median over CUBE_INDEX
    # AR now using ext=3, which computes the median on a "cleaner" sample
    # AR cutting on CONTRAST and N_SOURCES_FOR_PSF
    # AR see Aaron message https://desisurvey.slack.com/archives/C01HNN87Y7J/p1610865136043700?thread_ts=1610839476.023800&cid=C01HNN87Y7J
    gfa = fits.open(gfafn)[3].data
    # AR 20210109 : two exposures with "CMX LRG+QSO" instead of "SV1 LRG+QSO" (71594, 71595)
    # AR 20210110 : need to adapt...
    keep = np.array([program[:2] == "SV" for program in gfa["PROGRAM"]])
    keep |= np.array([program[:3] == "sv1" for program in gfa["PROGRAM"]])
    keep |= np.in1d(gfa["EXPID"], [71594, 71595])
    keep |= gfa["PROGRAM"] == "M31"
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
        # "SPECDATA_SKY_RMAG_AB",
        "SPECMODEL_SKY_RMAG_AB",
        #"BLANC_SPECMODEL_SKY_RMAG_AB",
        "NGFA",
        "B_DEPTH",
        "R_DEPTH",
        "Z_DEPTH",
        "B_DEPTH_EBVAIR",
        "R_DEPTH_EBVAIR",
        "Z_DEPTH_EBVAIR",
    ] + targets
    gfakeys = [
        "AIRMASS",
        "MOON_ILLUMINATION",
        "MOON_ZD_DEG",
        "MOON_SEP_DEG",
        "TRANSPARENCY",
        "FWHM_ASEC",
        "SKY_MAG_AB",
        "FIBER_FRACFLUX",
        "FIBER_FRACFLUX_ELG",
        "TRANSPFRAC",
    ]
    gfalabs = ["MIN", "MEAN", "MED", "MAX"]
    gfafuncs = [np.nanmin, np.nanmean, np.nanmedian, np.nanmax]
    allkeys = hdrkeys + ownkeys
    for key in gfakeys:
        allkeys += ["GFA_{}_{}".format(key, gfalab) for gfalab in gfalabs]
    for key in allkeys:
        exposures[key] = []
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
        # AR special dealing with 20210114, where positioners were frozen after the first
        # AR exposure (72381); hence we reject all the subsequent ones
        # AR https://desisurvey.slack.com/archives/C6C320XMK/p1610713799187700
        if night == 20210114:
            expids = np.array(["00072381"])
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
                    # AR SKY_RMAG_AB from integrating the sky fibers over the decam r-band
                    # exposures["SPECDATA_SKY_RMAG_AB"] += [
                    #    get_sky_rmag_ab(
                    #        night, expids[i], exposures["EXPTIME"][-1], "data"
                    #    )
                    # ]
                    # AR SKY_RMAG_AB from integrating the sky model over the decam r-band - daily
                    exposures["SPECMODEL_SKY_RMAG_AB"] += [
                        get_sky_rmag_ab(
                            night,
                            expids[i],
                            exposures["EXPTIME"][-1],
                            "model",
                            redux="daily",
                        )
                    ]
                    # AR GFA information
                    keep = gfa["EXPID"] == hdr["EXPID"]
                    exposures["NGFA"] += [keep.sum()]
                    if keep.sum() > 0:
                        for key in gfakeys:
                            # AR already considering median per CUBE_INDEX
                            # AR SKY_MAG_AB: converting to linear flux
                            if key == "SKY_MAG_AB":
                                x = 10.0 ** (-0.4 * (gfa["SKY_MAG_AB"][keep] - 22.5))
                            # AR TRANSPARENCY x FIBER_FRACFLUX
                            elif key == "TRANSPFRAC":
                                x = (
                                    gfa["TRANSPARENCY"][keep]
                                    * gfa["FIBER_FRACFLUX"][keep]
                                )
                            else:
                                x = gfa[key][keep]
                            # AR taking the min/mean/median/max
                            for gfalab, gfafunc in zip(gfalabs, gfafuncs):
                                # AR going back to mag for SKY_MAG_AB, after having done the stats
                                if key == "SKY_MAG_AB":
                                    exposures["GFA_{}_{}".format(key, gfalab)] += [
                                        22.5 - 2.5 * np.log10(gfafunc(x))
                                    ]
                                else:
                                    exposures["GFA_{}_{}".format(key, gfalab)] += [
                                        gfafunc(x)
                                    ]
                        # AR/DK exposure depths (needs gfa information)
                        # AR/DK adding also depth including ebv+airmass
                        depths_i = determine_tile_depth2(
                            exposures["NIGHT"][-1],
                            exposures["EXPID"][-1],
                            exposures["EXPTIME"][-1],
                            exposures["GFA_TRANSPFRAC_MEAN"][-1],
                        )
                        for camera in ["B", "R", "Z"]:
                            exposures["{}_DEPTH".format(camera)] += [
                                depths_i[camera.lower()]
                            ]
                            ebv = exposures["EBV"][-1]
                            fact_ebv = 10. ** (-2 * 0.4 * depth_coeffs["EBV"][camera] * ebv)
                            airmass = exposures["GFA_AIRMASS_MEAN"][-1]
                            fact_air = 10. ** (-2 * 0.4 * depth_coeffs["AIRMASS"][camera] * (airmass - 1.0))
                            exposures["{}_DEPTH_EBVAIR".format(camera)] += [
                                depths_i[camera.lower()] * fact_ebv * fact_air
                            ]
                        # for camera in ["B", "R", "Z"]: exposures["{}_DEPTH".format(camera)] += [-99]
                        # for camera in ["B", "R", "Z"]: exposures["{}_DEPTH_EBVAIR".format(camera)] += [-99]
                    else:
                        for key in gfakeys:
                            for gfalab in gfalabs:
                                exposures["GFA_{}_{}".format(key, gfalab)] += [-99]
                        for band in ["B", "R", "Z"]:
                            exposures["{}_DEPTH".format(band)] += [-99]
                            exposures["{}_DEPTH_EBVAIR".format(band)] += [-99]

    # AR if update=y, pre-append the results from previous nights
    if args.update == "y":
        d = fits.open(outfns["exposures"])[1].data
        keep = d["NIGHT"] < int(firstnight)
        d = d[keep]
        for key in allkeys:
            exposures[key] = d[key].tolist() + exposures[key]
    # AR building/writing fits
    cols = []
    for key in allkeys:
        if key in ["NIGHT", "EXPID", "TILEID", "NGFA"] + targets:
            fmt = "K"
        elif key in ["FIELD", "TARGETS"]:
            fmt = "{}A".format(np.max([len(x) for x in exposures[key]]))
        else:
            fmt = "E"
        cols += [fits.Column(name=key, format=fmt, array=exposures[key])]
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
            "[https://www.legacysurvey.org/viewer-dev/?ra={:.3f}&dec={:.3f}&layer=ls-dr9&zoom=8 Viewer]".format(
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
            tmparr += ["{:.2f}".format(di["GFA_" + key + "_MED"][j]) for key in keys]
            tmparr += [
                "{:.0f}s".format(di[band + "_DEPTH"][j]) for band in ["B", "R", "Z"]
            ]
            f.write("||{} ||\n".format(" ||".join(tmparr)))
    f.close()


# AR plots
if args.plot == "y":

    # AR sky map
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
            hpdict[key] = -99 + 0.0 * h[1].data[key]
            keep = h[1].data[key] > 0
            hpdict[key][keep] = np.log10(h[1].data[key][keep])
            hpdict[key + "lab"] = "log10(stardens)"
        elif (key[:8] == "galdepth") | (key[:8] == "psfdepth"):
            hpdict[key] = -99 + 0.0 * h[1].data[key]
            keep = h[1].data[key] > 0
            hpdict[key][keep] = 22.5 - 2.5 * np.log10(
                5.0 / np.sqrt(h[1].data[key][keep])
            )
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
    # AR
    for key in list(flavdict.keys()):
        ax.scatter(100, 100, marker="X", s=50, c=flavdict[key]["COLOR"], label=key)
    ax.legend(loc=2)
    plt.savefig(outfns["skymap"], bbox_inches="tight")
    plt.close()

    # AR observing conditions (one plot per month)
    keys = [
        "GFA_MOON_ILLUMINATION",
        "GFA_MOON_ZD_DEG",
        "GFA_MOON_SEP_DEG",
        "GFA_AIRMASS",
        "GFA_TRANSPARENCY",
        "GFA_FWHM_ASEC",
        "GFA_SKY_MAG_AB",
        "GFA_FIBER_FRACFLUX",
        "R_DEPTH / EXPTIME",
        "R_DEPTH_EBVAIR / EXPTIME"
    ]
    mlocs = [0.20, 25, 25, 0.20, 0.20, 0.50, 1.0, 0.20, 0.5, 0.5]
    ylims = [
        (0, 1),
        (0, 180),
        (0, 180),
        (0.9, 2),
        (0, 1.1),
        (0, 3),
        (17, 22),
        (0, 1),
        (0, 2.5),
        (0, 2.5)
    ]
    for month in months:
        d = fits.open(outfns["exposures"])[1].data
        keep = np.array(["{}".format(night // 100) for night in d["NIGHT"]]) == month
        if keep.sum() > 0:
            d = d[keep]
            xlim = (
                np.floor(d["MJDOBS"].min() - 1).astype(int),
                np.ceil(d["MJDOBS"].max()).astype(int) + 1,
            )
            cols = np.zeros(len(d), dtype=object)
            for fkey in list(flavdict.keys()):
                cols[d["TARGETS"] == fkey] = flavdict[fkey]["COLOR"]
            fig = plt.figure(figsize=(25, 1 * len(keys)))
            gs = gridspec.GridSpec(len(keys), 1, hspace=0.1)
            for i in range(len(keys)):
                ax = plt.subplot(gs[i])
                if keys[i] == "R_DEPTH / EXPTIME":
                    y = d["R_DEPTH"] / d["EXPTIME"]
                elif keys[i] == "R_DEPTH_EBVAIR / EXPTIME":
                    y = d["R_DEPTH_EBVAIR"] / d["EXPTIME"]
                else:
                    y = d["{}_MED".format(keys[i])]
                ylim = ylims[i]
                y = np.clip(y, ylim[0], ylim[1])
                ax.scatter(d["MJDOBS"], y, c=cols, marker="o", s=5)
                ax.grid(True)
                ax.set_axisbelow(True)
                ax.set_xlim(xlim)
                ax.yaxis.set_major_locator(MultipleLocator(mlocs[i]))
                ax.set_ylim(ylim)
                ax.text(
                    0.01,
                    0.80,
                    keys[i],
                    color="k",
                    fontweight="bold",
                    ha="left",
                    transform=ax.transAxes,
                )
                for x in range(xlim[0], xlim[1]):
                    ax.axvline(x, c="k", lw=0.1)
                if i == 0:
                    ax.set_title(
                        "SV1 observing conditions from {} to {} ({} exposures)".format(
                            d["NIGHT"].min(), d["NIGHT"].max(), len(d),
                        )
                    )
                    _, jj = np.unique(d["NIGHT"], return_index=True)
                    for j in jj:
                        ax.text(
                            np.floor(d["MJDOBS"][j]) + 0.5,
                            0.8 * ax.get_ylim()[1],
                            d["NIGHT"][j],
                            color="k",
                            ha="center",
                        )
                    for fkey in list(flavdict.keys()):
                        ax.scatter(
                            None,
                            None,
                            c=flavdict[fkey]["COLOR"],
                            marker="o",
                            s=5,
                            label=fkey,
                        )
                    ax.legend(loc=4)
                if i == len(keys) - 1:
                    ax.set_xlabel("MJD-OBS")
                    ax.set_xticks(np.linspace(xlim[0], xlim[1], xlim[1] - xlim[0] + 1))
                    ax.ticklabel_format(useOffset=False, style="plain")
                else:
                    ax.set_xticks([])
            plt.savefig(outfns["obscond"][month], bbox_inches="tight")
        plt.close()

    # AR r_depth
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
    for flavshort, ymax in zip(
        ["QSO+LRG", "QSO+ELG", "ELG", "BGS+MWS"], [2000, 2000, 2000, 1000]
    ):
        keep = d["TARGETS"] == flavshort
        if keep.sum() > 0:
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
                        "{}\n{}\n{:.0f}s".format(
                            d["NIGHT"][j], d["EXPID"][j], d[key][j]
                        ),
                        ha="center",
                        va="center",
                        fontsize=7,
                        color=cols[j],
                    )
                    ax.scatter(0.5 + x, d[key][j], c=cols[j], marker="o", s=5)
                    ax.plot(
                        [x, x + 1],
                        d["EXPTIME"][j] + np.zeros(2),
                        c="k",
                        ls="--",
                        lw=0.5,
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


# AR html pages
# initially copied from desitarget/QA.py from ADM!
if args.html == "y":
    d = fits.open(outfns["exposures"])[1].data
    tileids, ii = np.unique(tiles["TILEID"], return_index=True)
    nexps = np.array([(d["TILEID"] == tileid).sum() for tileid in tileids])

    # AR html settings
    bkgcol = "#f0f0f5"  # "#e6fff2"
    style = "line-height:25%; font-size:1vw"

    # ADM set up the html file and write preamble to it.
    htmlfile = os.path.join(outfns["html"], "index.html")

    # ADM grab the magic string that writes the last-updated date to a webpage.
    js = _javastring()

    # ADM html preamble.
    htmlmain = open(htmlfile, "w")
    htmlmain.write("<html><body>\n")
    htmlmain.write("<h1>SV1 Overview Page</h1>\n")
    htmlmain.write("\n")

    # AR page menu
    htmlmain.write("<nav>\n")
    htmlmain.write("\t<ul>\n")
    htmlmain.write("\t\t<li><a href='#fitsfiles' >Fits files</a></li>\n")
    htmlmain.write("\t\t<li><a href='#skymap' >Tiles sky map</a></li>\n")
    htmlmain.write(
        "\t\t<li><a href='#tile-nexp-design' >Tiles: NEXP and design</a></li>\n"
    )
    htmlmain.write(
        "\t\t<li><a href='#per-exposure-properties' > Per-exposure properties</a> ({} exposure(s) over {} night(s))</li>\n".format(
            len(d), len(np.unique(d["NIGHT"]))
        )
    )
    for month in months:
        if os.path.isfile(outfns["obscond"][month]):
            htmlmain.write(
                "\t\t<li><a href='#obsconds-{}' > Observing conditions for {}</a></li>\n".format(
                    month, month
                )
            )
    for flavshort in ["QSO+LRG", "ELG", "BGS+MWS", "QSO+ELG"]:
        htmlmain.write(
            "\t\t<li><a href='#depths-{}' > Per-tile exposure depths: {}</a></li>\n".format(
                flavshort, flavshort
            )
        )
    htmlmain.write("\t</ul>\n")
    htmlmain.write("</nav>\n")
    htmlmain.write("\n")

    # AR Fits files
    htmlmain.write("<h2><a id='fitsfiles' href='#fitsfiles' > Fits files</a>\n")
    htmlmain.write(
        "<a href='#top' style='position: absolute; right: 0;'>Top of the page</a></h2>\n"
    )
    htmlmain.write(
        "<p style='font-size:1vw'><a href='../sv1-exposures.fits' target='external'> sv1-exposures.fits </a> : main file with all exposures. Column content: {} ;  + GFA quantities (MIN, MEAN, MED, MAX for each quantity): {}.</p>\n".format(
            ", ".join([key for key in d.dtype.names if key[:4] != "GFA_"]),
            ", ".join(
                [
                    key.replace("GFA_", "").replace("_MIN", "")
                    for key in d.dtype.names
                    if key[:4] == "GFA_" and key[-3:] == "MIN"
                ]
            ),
        )
    )
    htmlmain.write(
        "<p style='font-size:1vw'><a href='../sv1-tiles.fits' target='external'> sv1-tiles.fits </a> : minimal file with the following tile informations: {}.</p>\n".format(
            ", ".join(fits.open(outfns["tiles"])[1].columns.names)
        )
    )

    # AR Sky map
    htmlmain.write("<h2><a id='skymap' href='#skymap' > Tiles: sky map</a>\n")
    htmlmain.write(
        "<a href='#top' style='position: absolute; right: 0;'>Top of the page</a></h2>\n"
    )
    tmppng = outfns["skymap"].replace(args.outdir, "../")
    htmlmain.write("<tr>\n")
    htmlmain.write(
        "<td align=center><a href='{}'><img SRC='{}' width=80% height=auto></a></td>\n".format(
            tmppng, tmppng
        )
    )
    htmlmain.write("</tr>\n")
    htmlmain.write("\n")

    # AR nexp + tile design
    _ = write_html_tiledesign(
        htmlmain, tiles, ii, nexps, style, "Tiles: NEXP and design", main=True
    )

    # AR Exposure properties
    _ = write_html_perexp(htmlmain, d, style, "Per-exposure properties")

    # AR Observing conditions
    for month in months:
        if os.path.isfile(outfns["obscond"][month]):
            htmlmain.write(
                "<h2><a id='obsconds-{}' href='#obsconds-{}' > Observing conditions for {}</a>\n".format(
                    month, month, month
                )
            )
            htmlmain.write(
                "<a href='#top' style='position: absolute; right: 0;'>Top of the page</a></h2>\n"
            )
            tmppng = outfns["obscond"][month].replace(args.outdir, "../")
            htmlmain.write("<tr>\n")
            htmlmain.write(
                "<td align=center><a href='{}'><img SRC='{}' width=100% height=auto></a></td>\n".format(
                    tmppng, tmppng
                )
            )
            htmlmain.write("</tr>\n")
    htmlmain.write("\n")

    # AR Depths per flavshort
    for flavshort in ["QSO+LRG", "QSO+ELG", "ELG", "BGS+MWS"]:
        htmlmain.write(
            "<h2><a id='depths-{}' href='#depths-{}' > Per-tile exposure depths: {}\n".format(
                flavshort, flavshort, flavshort
            )
        )
        htmlmain.write(
            "<a href='#top' style='position: absolute; right: 0;'>Top of the page</a></h2>\n"
        )
        tmppng = outfns["depth"][flavshort].replace(args.outdir, "../")
        htmlmain.write("<tr>\n")
        htmlmain.write(
            "<td align=center><a href='{}'><img SRC='{}' width=100% height=auto></a></td>\n".format(
                tmppng, tmppng
            )
        )
        htmlmain.write("</tr>\n")
    htmlmain.write("\n")

    # ADM for each tileid, make a separate page.
    for i in ii:
        tileid = tiles["TILEID"][i]
        jj = np.where(d["TILEID"] == tileid)[0]
        jj = jj[d["EXPID"][jj].argsort()]
        di = d[jj]

        # ADM call each page by the target class name, out it in the requested directory.
        htmlfile = os.path.join(outfns["html"], "{:06}.html".format(tileid))
        html = open(htmlfile, "w")

        # ADM html preamble.
        html.write("<html><body>\n")
        html.write("<h1>Tile {:06}</h1>\n".format(tileid))

        # AR page menu
        html.write("<nav>\n")
        html.write("\t<ul>\n")
        html.write(
            "\t\t<li><a href='index.html' target='external' > Back to the SV1 overview page </a></li>\n"
        )
        html.write(
            "\t\t<li><a href='#tile-nexp-design' >Tile {}: NEXP and design</a></li>\n".format(
                tileid
            )
        )
        html.write(
            "\t\t<li><a href='#per-exposure-properties' > Tile {}: per-exposure properties</a> ({} exposure(s) over {} night(s))</li>\n".format(
                tileid, len(di), len(np.unique(di["NIGHT"]))
            )
        )
        html.write(
            "\t\t<li><a href='#fa-qa-plot' > Tile {}: Fiber Assign QA plot</a></li>\n".format(
                tileid
            )
        )
        html.write("\t</ul>\n")
        html.write("</nav>\n")
        html.write("\n")

        # AR Tile design
        _ = write_html_tiledesign(
            html,
            tiles,
            [i],
            [len(di)],
            style,
            "Tile {}: NEXP and design".format(tileid),
            main=False,
        )

        # AR Exposure properties
        _ = write_html_perexp(
            html, di, style, "Tile {}: per-exposure properties".format(tileid)
        )

        # AR QA plot
        qapng = "https://desi.lbl.gov/svn/data/tiles/trunk/{}/fiberassign-{:06}.png".format(
            str(tiles["TILEID"][i]).zfill(6)[:3], tiles["TILEID"][i]
        )
        html.write(
            "<h2><a id='fa-qa-plot' href='#fa-qa-plot' > Tile {}: Fiber Assign QA plot</a>\n".format(
                tiles["TILEID"][i]
            )
        )
        html.write(
            "<a href='#top' style='position: absolute; right: 0;'>Top of the page</a></h2>\n"
        )
        html.write(
            "<td align=center><a href='{}'><img SRC='{}' width=100% height=auto></a></td>\n".format(
                qapng, qapng
            )
        )
        html.write("</tr>\n")

        # ADM html postamble
        html.write(
            "<p style='font-size:1vw; text-align:right'><i>Last updated: {}</p></i>\n".format(
                js
            )
        )
        html.write("</html></body>\n")
        html.close()

    # ADM html postamble for main page.
    htmlmain.write(
        "<p style='font-size:1vw; text-align:right'><i>Last updated: {}</p></i>\n".format(
            js
        )
    )
    htmlmain.write("</html></body>\n")
    htmlmain.close()

    # ADM make sure all of the relevant directories and plots can be read by a web-browser.
    cmd = "chmod 644 {}/*".format(outfns["html"])
    ok = os.system(cmd)
    cmd = "chmod 775 {}".format(outfns["html"])
    ok = os.system(cmd)
