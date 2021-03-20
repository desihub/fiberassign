#!/usr/bin/env python

import matplotlib

matplotlib.use("Agg")
import sys
import os
from glob import glob
import numpy as np
import fitsio
import astropy.io.fits as fits
from astropy.table import Table, vstack
from astropy import units, constants
from astropy.coordinates import SkyCoord
from astropy.time import Time
import pytz
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
import datetime
import desisurvey
from desisurvey.scripts.collect_etc import get_conditions
import multiprocessing
from desitarget.internal import sharedmem

# AR https://github.com/desihub/desicmx/blob/master/analysis/gfa/DESI-Online-Database-Tutorial.ipynb
try:
    import desietcimg
except ImportError:
    sys.exit(
        "desietcimg not installed\nto install: {} -m pip install --user git+https://github.com/dkirkby/desietcimg\nexiting".format(
            sys.executable
        )
    )
from desietcimg.db import DB
from desispec.io.fluxcalibration import read_average_flux_calibration
from argparse import ArgumentParser


## AR nb of processors
# nproc = multiprocessing.cpu_count() // 2


# AR reading arguments
parser = ArgumentParser()
parser.add_argument(
    "--outdir",
    help="output directory (fiducial location: $DESI_ROOT/survey/observations/SV1/)",
    type=str,
    default=None,
    required=True,
)
parser.add_argument(
    "--tiles",
    help="merge the TILEID-tiles.fits into one file? (y/n)",
    type=str,
    default="y",
)
parser.add_argument(
    "--exposures", help="do fits with per-exposure stats? (y/n)", type=str, default="y",
)
parser.add_argument(
    "--plot",
    help="do plots (skymap, observing conditions, r_depth)? (y/n)",
    type=str,
    default="y",
)
parser.add_argument(
    "--html", help="create html pages, one per tile (y/n)", type=str, default="y",
)
parser.add_argument(
    "--update",
    help="start from existing args.outroot+'sv1-exposures.fits' (y/n)",
    type=str,
    default="y",
)
parser.add_argument(
    "--numproc",
    help="number of concurrent processes to use (default=1)",
    type=int,
    default=1,
)
parser.add_argument(
    "--refspecprod",
    help="specprod to use (falls back to daily if exposure not in specprod; (default=cascades)",
    type=str,
    default="cascades",
)
args = parser.parse_args()
for kwargs in args._get_kwargs():
    print(kwargs)


# AR : all exposure depths routines are copied from DK, with minor modifications

# AR safe
if args.outdir[-1] != "/":
    args.outdir += "/"

# AR folders / files
rawdir = os.path.join(os.getenv("DESI_ROOT"), "spectro", "data")
surveydir = os.path.join(os.getenv("DESI_ROOT"), "survey")
nightoutdir = os.path.join(args.outdir, "sv1-nights")
if not os.path.isdir(nightoutdir):
    os.mkdir(nightoutdir)
ephemdir = os.path.join(
    surveydir, "observations", "misc", "ephem-bgd-jan2021"
)  # jan2021: using dark_max_sun_altitude=-18, bright_max_sun_altitude=-12
# ephemdir = os.path.join(surveydir, "observations", "misc", "ephem-bgd-master")
pixwfn = (
    os.getenv("DESI_TARGET")
    + "/catalogs/dr9/0.47.0/pixweight/sv1/resolve/dark/sv1pixweight-dark.fits"
)
desfn = os.path.join(surveydir, "observations", "misc", "des_footprint.txt")
dbyamlfn = os.path.join(surveydir, "observations", "misc", "db.yaml")
thrufn = {
    camera: "{}/thru-fluxcalib/v2/cascades-transpfrac-airmass-noseeing/fluxcalib-cascades-{}-0.58transpfrac0.78.fits".format(
        os.getenv("CSCRATCH"), camera
    )
    for camera in ["b", "r", "z"]
}

# AR SKY MONITOR measurements
db = DB(config_name=dbyamlfn)
tmpout = db.query(
    "select time_recorded,skycam0,skycam1,average from sky_skylevel", maxrows=1e10
)
skymon = {}
skymon["MJD"] = np.array(
    [
        tmpout["time_recorded"][i].to_julian_date() - 2400000.5
        for i in range(len(tmpout))
    ]
)
for key in ["SKYCAM0", "SKYCAM1", "AVERAGE"]:
    skymon[key] = tmpout[key.lower()].values

# AR GFA
# AR Feb. 01, 2021: now using the *matched_coadd* file [desi-survey 1302]
# AR Feb. 09, 2021: also using the *acq* file to provide measurements
# AR                for exposures not present in *matched_coadd*
# AR                https://desisurvey.slack.com/archives/C01HNN87Y7J/p1612806269343500
gfafn = {}
gfafn["acq"] = np.sort(
    glob(surveydir + "/GFA/offline_acq_ccds_SV1-thru_20??????.fits")
)[-1]
gfafn["matched_coadd"] = np.sort(
    glob(surveydir + "/GFA/offline_matched_coadd_ccds_SV1-thru_20??????.fits")
)[-1]
print("gfafn = {}".format(gfafn))
gfa = {key: fits.open(gfafn[key])[3].data for key in list(gfafn.keys())}

# AR ephemerides:
# AR using {ephemdir}/config.yaml to define bright/gray/dark
# AR ephem-bgd-jan2021/config.yaml should be similar to NTS as of Jan of 2021
os.environ["DESISURVEY_OUTPUT"] = ephemdir
# AR mimicking desisurvey.ephem.get_ephem()
# AR keeping the same start/stop dates as in desisurvey, as those
# AR      are kind of hard-coded there
# AR      and are used in get_obsconditions()
START_DATE, STOP_DATE = datetime.date(2019, 1, 1), datetime.date(2025, 12, 31)
start_iso, stop_iso = START_DATE.isoformat(), STOP_DATE.isoformat()
ephem_config = desisurvey.config.Configuration(
    file_name=os.path.join(ephemdir, "config.yaml")
)
ephem_fn = ephem_config.get_path("ephem_{}_{}.fits".format(start_iso, stop_iso))
desisurvey.utils.freeze_iers()
# AR should not happen, unless new settings are entered
if not os.path.isfile(ephem_fn):
    print("Generating {}".format(ephem_fn))
    _ephem = desisurvey.ephem.Ephemerides(START_DATE, STOP_DATE)
    _ephem._table.write(ephem_fn)
print("Reading {}".format(ephem_fn))
ephem = desisurvey.ephem.Ephemerides(START_DATE, STOP_DATE, restore=ephem_fn)


# AR for the fiberflat routine
os.environ["DESI_LOGLEVEL"] = "warning"


# AR output products
outfns = {}
outfns["tiles"] = os.path.join(args.outdir, "sv1-tiles.fits")
outfns["exposures"] = os.path.join(args.outdir, "sv1-exposures.fits")
if os.path.isdir(os.path.join(args.outdir, "sv1-plots")) == False:
    os.mkdir(os.path.join(args.outdir, "sv1-plots"))
outfns["skymap"] = os.path.join(args.outdir, "sv1-plots", "sv1-skymap.png")
outfns["onsky"] = os.path.join(args.outdir, "sv1-plots", "sv1-onsky.png")
for flavshort in ["QSO+LRG", "ELG", "BGS+MWS", "QSO+ELG", "SSV"]:
    fsh = flavshort.lower().replace("+", "")
    outfns["efftime-{}".format(fsh)] = os.path.join(
        args.outdir, "sv1-plots", "sv1-efftime-{}.png".format(fsh)
    )
months = ["202012", "202101", "202102", "202103", "202104"]  # AR for splitting plots
outfns["obscond"] = {}
for month in months:
    outfns["obscond"][month] = os.path.join(
        args.outdir, "sv1-plots", "sv1-obscond-{}.png".format(month)
    )
outfns["obscond"]["cumul"] = os.path.join(
    args.outdir, "sv1-plots", "sv1-obscond-cumul.png"
)
outfns["depth"] = {}
for flavshort in ["QSO+LRG", "ELG", "BGS+MWS", "QSO+ELG", "test-PETAL_LOC_0", "SSV"]:
    outfns["depth"][flavshort] = os.path.join(
        args.outdir,
        "sv1-plots",
        "sv1-depth-{}.png".format(flavshort.lower().replace("+", "")),
    )
outfns["depth"]["cumul"] = os.path.join(
    args.outdir, "sv1-plots", "sv1-depth-cumul.png",
)
outfns["html"] = os.path.join(args.outdir, "sv1-html")
if (args.html == "y") & (os.path.isdir(outfns["html"]) == False):
    os.mkdir(outfns["html"])

# AR sv1 first night (for exposures search)
# AR sv1 ref_specprod last night
firstnight = 20201214
if args.refspecprod == "cascades":
    refspecprod_lastnight = 20210224
print("firstnight = {}".format(firstnight))
print("{}, lastnight={}".format(args.refspecprod, refspecprod_lastnight))

# AR for the exposure and the wiki cases
targnames = ["TGT", "SKY", "STD", "WD", "LRG", "ELG", "QSO", "BGS", "MWS"]
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
    "QSO+LRG": {"FAFLAVORS": ["cmxlrgqso", "sv1lrgqso", "sv1lrgqso2"], "COLOR": "r"},
    "ELG": {"FAFLAVORS": ["cmxelg", "sv1elg"], "COLOR": "b"},
    "QSO+ELG": {"FAFLAVORS": ["sv1elgqso"], "COLOR": "c"},
    "BGS+MWS": {"FAFLAVORS": ["cmxbgsmws", "sv1bgsmws"], "COLOR": "g"},
    "M33+Dark": {"FAFLAVORS": ["cmxm33"], "COLOR": "magenta"},
    "M31": {"FAFLAVORS": ["sv1m31"], "COLOR": "y"},
    "test-PETAL_LOC_0": {
        "FAFLAVORS": ["sv1orion", "sv1rosette", "sv1praesepe", "sv1umaii"],
        "COLOR": "0.5",
    },
    "SSV": {"FAFLAVORS": ["sv1ssv"], "COLOR": "orange"},
    "SCND": {"FAFLAVORS": ["sv1mwclusgaldeep", "sv1scndhetdex", "sv1unwisegreen", "sv1unwisebluebright", "sv1unwisebluefaint", "sv1scndcosmos", "sv1unwise"], "COLOR": "brown"},
    "BACKUP": {"FAFLAVORS": ["sv1backup1"], "COLOR": "y"},
}
# AR field names - hand-written...
fielddict = {
    "XMM-LSS": [80605, 80606],
    "Lynx": [80607, 80608, 80613],
    "COSMOS": [80609, 80610, 80742, 80871, 80872],
    "Triangulum": [80611],
    "Eridanus": [80612],
    "Sextans": [80614, 80737],
    "Triangulum-CMX": [80615],
    "Pegasus1": [80616],
    "Pegasus2": [80617],
    "NGC2419": [80618, 80721],
    "UMajor": [80619, 80620, 80621],
    "LeoMinor": [80622, 80623],
    # BGS+MWS 20210101 , Dark 20210105
    "Dust cloud": [80624],
    "Far north": [80650, 80655, 80693, 80694],
    "Far south": [80630, 80638, 80671, 80701, 80672, 80702],
    "GAMA G02": [80633, 80635],
    "GAMA G09": [80739, 80740, 80741],
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
    "Orion": [80717],
    "Rosette": [80718],
    "Prasepe": [80719, 80723],
    "UMAII": [80720, 80726],
    "GD1_LOW_1": [80722],
    "GD1_LOW_2": [80724],
    "M67": [80725, 80864],
    "GD1-B-1": [80727],
    "ORPHAN-A-2": [80728],
    "GD1-C-1": [80729],
    "GD1-C-3": [80730],
    "NGP": [80731],
    "BOSS7456": [80732],
    "M53+N5053": [80733, 80863],
    "M5": [80734],
    "M13": [80735],
    "M92": [80736],
    "Draco": [80738, 80862],
    "ELAIS_N1": [80865, 80866, 80867],
    "CFHTLS_W3": [80868],
    "HETDEX": [80869, 80870]
}


# AR/DK DESI spectra wavelengths
wmin, wmax, wdelta = 3600, 9824, 0.8
fullwave = np.round(np.arange(wmin, wmax + wdelta, wdelta), 1)
cslice = {"b": slice(0, 2751), "r": slice(2700, 5026), "z": slice(4900, 7781)}
# AR (wmin,wmax) to "stich" all three cameras
wstich = {"b": (wmin, 5780), "r": (5780, 7570), "z": (7570, 9824)}


# AR DESI telescope geometric area (cm2) and fiber area (arcsec2)
# AR for converting the spectroscopic sky to correct units
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
det_eso = load_spec(os.path.join(surveydir, "observations", "misc", "dark_eso.fits"))
det_desimodel = load_spec(
    os.path.join(surveydir, "observations", "misc", "dark_desimodel.fits")
)
_sky_cache = {}


# AR/ES ebv and airmass coefficient in depth_ebvair
# AR/ES ebv coeffs : taking SDSS grz from Schlafly & Finkbeinger (2011)
# AR/ES airmass coeffs : [decam-chatter 15497]
depth_coeffs = {
    "EBV": {"B": 3.303, "R": 2.285, "Z": 1.263},
    "AIRMASS": {"B": 0.195, "R": 0.096, "Z": 0.055},
}


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
def get_sky(night, expid, specprod, specs=range(10)):
    """
    Estimate the sky spectrum for one exposure in units of phot/sec per wavelength bin.
    Returns a tuple (flux_inc, ivar_inc, flux_det, ivar_det) where "det" is detected phot/sec
    in each camera with flat-field corrections applied, and "inc" corrects for the average
    spectrograph throughput in each camera, then coadds over cameras.
    """
    # AR specprod : cascades, daily
    incident = CoAdd("full")
    detected = {}
    # Loop over cameras.
    for camera in ["b", "r", "z"]:
        detected[camera] = CoAdd(camera)
        # Loop over spectrographs.
        for spec in specs:
            # Read the flat-fielded (constant) sky model in this spectrograph.
            skypath = os.path.join(
                os.getenv("DESI_ROOT"),
                "spectro",
                "redux",
                specprod,
                "exposures",
                "{}".format(night),
                "{:08}".format(expid),
                "sky-{}{}-{:08}.fits".format(camera, spec, expid),
            )
            if not os.path.isfile(skypath):
                print("\t\tSkipping non-existent {}.".format(skypath))
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
    return incident, detected


# AR/DK exposure depths utilities
# AR transpfrac = transparency * fiber_fracflux
# AR 2021-03-05: updating ffracref=0.582, as GFA now provides measurement for fiber_diameter=1.52"
def determine_tile_depth2(
    night,
    expid,
    specprod,
    exptime,
    transpfrac,
    darkref=det_eso,
    ffracref=0.582,
    smoothing=125,
):
    inc, det = get_sky(night, expid, specprod)
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


# AR grz-band sky mag / arcsec2 from sky-....fits files
# AR now using work-in-progress throughput
# AR still provides a better agreement with GFAs than previous method
def get_sky_grzmag_ab(night, expid, specprod, exptime, ftype, fiber=0):
    # AR ftype = "data" or "model"
    # AR specprod = cascades, daily
    # AR if ftype = "data" : read the sky fibers from frame*fits + apply flat-field
    # AR if ftype = "model": read the sky model from sky*.fits for the first fiber of each petal (see DJS email from 29Dec2020)
    # AR those are in electron / angstrom
    if ftype not in ["data", "model"]:
        sys.exit("ftype should be 'data' or 'model'")
    sky = np.zeros(len(fullwave))
    # AR looking for a petal with brz sky and ivar>0
    spec, ncam = -1, 0
    while (ncam < 3) & (spec < 9):
        spec += 1
        tmpfns = glob(
            os.path.join(
                os.getenv("DESI_ROOT"),
                "spectro",
                "redux",
                specprod,
                "exposures",
                "{}".format(night),
                "{:08d}".format(expid),
                "sky-?{}-{:08d}.fits".format(spec, expid),
            )
        )
        ncam = 0
        for tmpfn in tmpfns:
            if fits.open(tmpfn)["IVAR"].data.max() > 0:
                ncam += 1
    if ncam < 3:
        print(
            "Skipping {}-{:08d} : no petal with sky for the three cameras".format(
                night, expid
            )
        )
        return -99, -99, -99
    else:
        for camera in ["b", "r", "z"]:
            # AR model
            if ftype == "model":
                fn = os.path.join(
                    os.getenv("DESI_ROOT"),
                    "spectro",
                    "redux",
                    specprod,
                    "exposures",
                    "{}".format(night),
                    "{:08d}".format(expid),
                    "sky-{}{}-{:08d}.fits".format(camera, spec, expid),
                )
                h = fits.open(fn)
                acal = read_average_flux_calibration(thrufn[camera])
                flux = np.interp(
                    fullwave[cslice[camera]], h["WAVELENGTH"].data, h["SKY"].data[fiber]
                )
            sky[cslice[camera]] = flux / exptime / acal.throughput()
    # AR sky model flux in erg / angstrom / s (using the photon energy in erg)
    e_phot_erg = (
        constants.h.to(units.erg * units.s)
        * constants.c.to(units.angstrom / units.s)
        / (fullwave * units.angstrom)
    )
    sky *= e_phot_erg.value
    # AR sky model flux in erg / angstrom / s / cm**2 / arcsec**2
    sky /= telap_cm2 * fiber_area_arcsec2
    # AR integrate over the DECam grz-bands
    # AR using the curves with no atmospheric extinction (email from DK from 09Feb2021)
    filts = filters.load_filters(
        "decamDR1noatm-g", "decamDR1noatm-r", "decamDR1noatm-z"
    )
    # filts = filters.load_filters("decam2014-g", "decam2014-r", "decam2014-z")
    # AR zero-padding spectrum so that it covers the DECam grz passbands
    # AR looping through filters while waiting issue to be solved (https://github.com/desihub/speclite/issues/64)
    sky_pad, fullwave_pad = sky.copy(), fullwave.copy()
    for i in range(len(filts)):
        sky_pad, fullwave_pad = filts[i].pad_spectrum(
            sky_pad, fullwave_pad, method="zero"
        )
    return filts.get_ab_magnitudes(
        sky_pad * units.erg / (units.cm ** 2 * units.s * units.angstrom),
        fullwave_pad * units.angstrom,
    ).as_array()[0]


# AR
# AR https://desi.lbl.gov/trac/wiki/SurveyOps/SurveySpeed
# AR 2021-03-05: GFA values are now computed for a fiber diameter of 1.52"
def get_efftime_speed(
    exptime,
    sky,
    ebv,
    transparency,
    airmass,
    ffrac_psf,
    ffrac_elg,
    ffrac_bgs,
    ffrac_psf_nom,
    ffrac_elg_nom,
    ffrac_bgs_nom,
    kterm,
    ebv_r_coeff=2.165,
):
    # AR : sky <- spectro. sky in nMgy / arcsec**2 (corrected by throughput)
    # AR : transparency, airmass <- GFA
    # AR : ffrac_psf <- GFA fiber_fracflux
    # AR : ffrac_elg <- GFA fiber_fracflux_elg
    # AR : ffrac_bgs <- GFA fiber_fracflux_bgs
    # AR : ffrac_*_nom and kterm <- from GFA [ffrac_psf_nom~0.582, ffrac_elg_nom~0.424, ffrac_bgs_nom~0.195]
    #
    # AR nominal values
    exptime_nom = 1000.0  # AR seconds
    sky_nom = 3.73  # AR nMgy/arcsec**2
    fflux_bright_nom = (
        15.8  # AR nMgy/arcsec**2 (r=19.5 mag for de Vaucouleurs rhalf=1.5" BGS)
    )
    fflux_backup_nom = 27.5  # AR nMgy/arcsec**2 (r=18.9 mag star)
    # AR airmass term
    airfac = 10.0 ** (kterm * (airmass - 1.0) / 2.5)
    # AR ebv term
    ebvfac = 10.0 ** (ebv_r_coeff * ebv / 2.5)
    # AR sky readnoise
    sky_rdn = 0.932  # AR nMgy/arcsec**2
    # AR "limit" fiber flux
    fflux_bright = (
        ffrac_bgs * transparency / airfac * fflux_bright_nom / fiber_area_arcsec2
    )
    fflux_backup = (
        ffrac_psf * transparency / airfac * fflux_backup_nom / fiber_area_arcsec2
    )
    # AR effective sky
    effsky_dark = (sky + sky_rdn * exptime_nom / exptime) / (1.0 + sky_rdn / sky_nom)
    effsky_bright = (sky + sky_rdn * exptime_nom / exptime + fflux_bright) / (
        1.0 + sky_rdn / sky_nom + fflux_bright / sky_nom
    )
    effsky_backup = (sky + sky_rdn * exptime_nom / exptime + fflux_backup) / (
        1.0 + sky_rdn / sky_nom + fflux_backup / sky_nom
    )
    # AR effective exposure time
    efftime_dark = (
        exptime
        * (ffrac_elg / ffrac_elg_nom * transparency / airfac) ** 2
        * (sky_nom / effsky_dark)
        / ebvfac ** 2
    )
    efftime_bright = (
        exptime
        * (ffrac_bgs / ffrac_bgs_nom * transparency / airfac) ** 2
        * (sky_nom / effsky_bright)
        / ebvfac ** 2
    )
    efftime_backup = (
        exptime
        * (ffrac_psf / ffrac_psf_nom * transparency / airfac) ** 2
        * (sky_nom / effsky_backup)
        / ebvfac ** 2
    )
    # AR survey speed
    speed_dark = (ffrac_elg / ffrac_elg_nom * transparency) ** 2 * (sky_nom / sky)
    speed_bright = (ffrac_bgs / ffrac_bgs_nom * transparency) ** 2 * (sky_nom / sky)
    speed_backup = (ffrac_psf / ffrac_psf_nom * transparency) ** 2 * (sky_nom / sky)
    # AR
    return (
        efftime_dark,
        efftime_bright,
        efftime_backup,
        speed_dark,
        speed_bright,
        speed_backup,
    )


# AR get the corresonding night in ephem
# AR copying https://github.com/desihub/desisurvey/blob/f83e098e0b31c0d496a6f240a6aa68131191ec76/py/desisurvey/scripts/collect_etc.py#L229-L258
def get_nightind(mjd):
    # taken in 2019 or later, with known MJD---removes some test exposures
    okmjd = np.isfinite(mjd) & (mjd > 58484)
    nights = ephem.table
    indices = np.repeat(np.arange(len(nights)), 2)
    startstop = np.concatenate(
        [[night["brightdusk"], night["brightdawn"]] for night in nights]
    )
    nightind = np.zeros(len(mjd), dtype="f8")
    nightind[~okmjd] = -1
    nightind[okmjd] = np.interp(mjd[okmjd], startstop, indices)
    # a lot of exposures apparently taken during the day?
    # We should mark these as belonging to the "day" program?
    mday = nightind != nightind.astype("i4")
    nightind = nightind.astype("i4")
    return nightind


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
    ] + targnames
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
    html.write("th, td {border:1px solid black; font-size: 0.90em}\n")
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
            "<a href='https://www.legacysurvey.org/viewer-dev/?&layer=ls-dr9&zoom=8&tile={}' target='external'> Viewer".format(
                tiles["TILEID"][i]
            )
        ]
        tmparr += ["{}".format(tiles[targname][i]) for targname in targnames]
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
    html.write("<p style='{}'>GFA acq fits file: {}</p>".format(style, gfafn["acq"]))
    html.write(
        "<p style='{}'>GFA matched_coadd fits file: {}</p>".format(
            style, gfafn["matched_coadd"]
        )
    )
    html.write(
        "<p style='{}'>FIBER_FRACFLUX: fraction of light in a 1.52 arcsec diameter fiber-sized aperture given the PSF shape, assuming that the PSF is perfectly aligned with the fiber (i.e. does not capture any astrometry/positioning errors).</p>".format(
            style
        )
    )
    html.write(
        "<p style='{}'>AIRMASS, MOON_SEP_DEG, TRANSPARENCY, FWHM_ASEC, FIBER_FRACFLUX: from GFA.</p>".format(
            style
        )
    )
    html.write(
        "<p style='{}'>SPECTRO_SKY: spectroscopic sky, corrected for throughput, convolved with DECam r-band filter.</p>".format(
            style
        )
    )
    html.write(
        "<p style='{}'>EFFTIME_DARK, SPEED_DARK: see https://desi.lbl.gov/trac/wiki/SurveyOps/SurveySpeed.</p>".format(
            style
        )
    )
    # ADM write out a list of the target categories.
    fields = (
        ["TILEID", "TARGETS", "NIGHT", "EXPID", "NIGHTWATCH", "EXPTIME", "EFFTIME_DARK", "SPEED_DARK"]
        + ["EBV", "GFA_AIRMASS", "GFA_MOON_SEP_DEG", "GFA_TRANSPARENCY", "GFA_FWHM_ASEC", "GFA_FIBER_FRACFLUX"]
        + ["SPECTRO_SKY"]
    )
    html.write("<table>\n")
    night = ""
    for i in range(len(d))[::-1]:
        night_prev = night
        night = d["NIGHT"][i]
        specprod = d["SPECPROD"][i]
        # AR night header
        if night != night_prev:
            html.write("<tr>\n")
            html.write(" ".join(["<th> {} </th>".format(x.replace("GFA_", "")) for x in fields]) + "\n")
            html.write("</tr>\n")
        html.write("<tr>")
        # AR redux path
        redux_path = os.path.join(
                    "https://data.desi.lbl.gov/desi",
                    "spectro",
                    "redux",
                    specprod,
        )
        # AR building array
        tmparr = []
        for field in fields:
            if field == "TILEID":
                tmparr += [
                    "<a href='{}' target='external'> {}".format(
                    os.path.join(
                        redux_path,
                        specprod,
                        "tiles",
                        "{}".format(d["TILEID"][i]),
                    ),
                    d["TILEID"][i],
                    )
                ]
            elif field == "NIGHT":
                tmparr += [
                    "<a href='{}' target='external'> {}".format(
                        os.path.join(
                            redux_path,
                            "exposures",
                            "{}".format(d["NIGHT"][i]),
                        ),
                        d["NIGHT"][i],
                    )
                ]
            elif field == "EXPID":
                tmparr += [
                    "<a href='{}' target='external'> {}".format(
                        os.path.join(
                            redux_path,
                            "exposures",
                            "{}".format(d["NIGHT"][i]),
                            "{:08}".format(d["EXPID"][i]),
                        ),
                        d["EXPID"][i],
                    )
                ]
            elif field == "NIGHTWATCH":
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
            elif field == "SPECTRO_SKY":
                tmparr += ["{:.1f}".format(22.5-2.5*np.log10(d["SPECMODEL_SKY_RFLUX"][i]))]
            elif field in ["TARGETS"]:
                tmparr += ["{}".format(d[field][i])]
            elif field in ["SPEED_DARK"]:
                tmparr += ["{:.1f}".format(d[field][i])]
            elif field in ["EBV", "GFA_AIRMASS", "GFA_TRANSPARENCY", "GFA_FWHM_ASEC", "GFA_FIBER_FRACFLUX"]:
                tmparr += ["{:.2f}".format(d[field][i])]
            else:
                tmparr += ["{:.0f}".format(d[field][i])]
        html.write(" ".join(["<td> {} </td>".format(x) for x in tmparr]) + "\n")
        html.write("</tr>\n")
    html.write("</table>\n")
    html.write("\n")
    return True


# AR per-tile information
tiles = {}
tiles["FN"] = np.sort(
    glob(
        os.path.join(
            surveydir, "fiberassign", "SV1", "202?????/fiberassign-??????.fits.gz"
        )
    )
)
nt = len(tiles["FN"])
# AR initialising
for key in [
    "TILEID",
    "TILERA",
    "TILEDEC",
    "FAFLAVOR",
    "TARGETS",
    "COLOR",
    "FIELD",
    "RUNDATE",
    "DESIGNDATE"
]:
    if key in ["TILEID", "DESIGNDATE"]:
        tiles[key] = np.zeros(nt, dtype=int)
    elif key in ["FAFLAVOR", "TARGETS", "COLOR", "FIELD", "RUNDATE"]:
        tiles[key] = np.array(["-" for x in range(nt)], dtype=object)
    else:
        tiles[key] = np.zeros(nt, dtype=float)
for targname in targnames:
    tiles[targname] = np.zeros(nt, dtype=int)
# AR populating
for i in range(nt):
    # AR DESIGNDATE
    tiles["DESIGNDATE"][i] = int(tiles["FN"][i].split("/")[-2])
    # AR general
    hdr = fits.getheader(tiles["FN"][i])
    for key in ["TILEID", "TILERA", "TILEDEC", "FAFLAVOR", "RUNDATE"]:
        tiles[key][i] = hdr[key]
    # AR number of targets per tracer
    d = fits.open(tiles["FN"][i])[1].data
    if tiles["FAFLAVOR"][i] == "cmxm33":
        mask, key, msks, std_msks = cmx_mask, "cmx_target", cmx_msks, std_cmx_msks
    else:
        mask, key, msks, std_msks = desi_mask, "sv1_desi_target", sv1_msks, std_sv1_msks
    for targname, msk in zip(targnames, msks):
        if targname in ["TGT", "SKY"]:
            tiles[targname][i] = (d["objtype"] == targname).sum()
        elif targname == "STD":
            keep = np.zeros(len(d), dtype=bool)
            for std_msk in std_msks:
                keep |= (d[key] & mask[std_msk]) > 0
            tiles[targname][i] = keep.sum()
        else:
            tiles[targname][i] = ((d[key] & mask[msk]) > 0).sum()
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
        glob(
            os.path.join(
                surveydir,
                "fiberassign",
                "SV1",
                "202?????",
                "{:06}-tiles.fits".format(tileid),
            )
        )[0]
        for tileid in tiles["TILEID"]
    ]
    h = fits.open(fns[0])
    keys, fmts = h[1].columns.names, h[1].columns.formats
    ownkeys, ownfmts = ["FAFLAVOR", "TARGETS", "FIELD", "RUNDATE", "DESIGNDATE"], ["-", "-", "-", "-", "K"]
    t = {}
    for key in keys + ownkeys:
        t[key] = []
    for fn in fns:
        d = fits.open(fn)[1].data
        for key in keys:
            t[key] += [d[key]]
        i = np.where(tiles["TILEID"] == d["TILEID"])[0][0]
        for key in ownkeys:
            t[key] += [tiles[key][i]]
    # AR building/writing fits
    cols = []
    for key, fmt in zip(keys + ownkeys, fmts + ownfmts):
        if key in ["FAFLAVOR", "TARGETS", "FIELD", "RUNDATE"]:
            fmt = "{}A".format(np.max([len(x) for x in t[key]]))
        cols += [fits.Column(name=key, format=fmt, array=t[key])]
    h = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
    h.writeto(outfns["tiles"], overwrite=True)


# AR output file name for a night
def get_night_outfn(nightoutdir, night):
    return "{}/sv1-exposures-{}.fits".format(nightoutdir, night)


# AR listing the expids for that night
# AR Feb. 05, 2021 : changing the exposure listing approach
# AR                 hopefully more reliable...
def get_night_expids(night, rawdir, tiles):
    nightrawdir = "{}/{}/".format(rawdir, night)
    # AR planned
    fns = np.sort(
        glob(os.path.join(nightrawdir, "????????", "fiberassign-??????.fits*"))
    )
    plan_expids = np.array([fn.split("/")[-2] for fn in fns])
    plan_tileids = np.array(
        [fn.split("/")[-1].replace("-", ".").split(".")[1] for fn in fns]
    )
    # AR taking unique because multiple files for some expids of 20201214
    # AR because of the presence of some fiberassign-TILEID.fits.orig files
    _, ii = np.unique(plan_expids, return_index=True)
    plan_expids = plan_expids[ii]
    plan_tileids = plan_tileids[ii]
    # AR cutting on our science TILEIDs (e.g. removing dithering)
    keep = np.in1d(plan_tileids.astype(int), tiles["TILEID"])
    plan_expids, plan_tileids = plan_expids[keep], plan_tileids[keep]
    # AR observed
    fns = np.sort(glob(os.path.join(nightrawdir, "????????", "desi-????????.fits.fz")))
    obs_expids = np.array([fn.split("/")[-2] for fn in fns])
    # AR keeping planned + observed
    keep = np.in1d(plan_expids, obs_expids)
    expids = plan_expids[keep].astype(int)
    # AR special dealing with 20210114, where positioners were frozen after the first
    # AR exposure (72381); hence we reject all the subsequent ones
    # AR https://desisurvey.slack.com/archives/C6C320XMK/p1610713799187700
    if night == 20210114:
        expids = np.array([72381])
    print("{} : {} exposures".format(night, len(expids)))
    return expids


# AR processing exposures from a given night
# AR per-exposure information (various from header, skymon, gfas + depths)
def process_night(night, nightoutdir, skymon, gfa, ephem, rawdir, tiles):
    # AR output file name (removing if already existing)
    outfn = get_night_outfn(nightoutdir, night)
    if os.path.isfile(outfn):
        os.remove(outfn)
    # AR using args.refspecprod for night <= refspecprod_lastnight
    # AR       daily otherwise
    if night <= refspecprod_lastnight:
        specprod = args.refspecprod
    else:
        specprod = "daily"
    # AR listing the expids for that night
    expids = get_night_expids(night, rawdir, tiles)
    nexp = len(expids)
    print("processing night={} ({} exposures)".format(night, nexp))
    if nexp > 0:
        # AR quantities we store
        # AR listing all camera+petal, to check existing sky/sframe/cframe files in the reduction
        campet_names = []
        for camera in ["b", "r", "z"]:
            campet_names += [camera + p for p in np.arange(10, dtype=int).astype(str)]
        campet_bits = np.arange(len(campet_names), dtype=int)
        # AR from exposure header
        hdrkeys = ["TILEID", "TILERA", "TILEDEC", "EXPTIME", "MJDOBS"]
        # AR own keys
        ownkeys = [
            "NIGHT",
            "EXPID",
            "FIELD",
            "TARGETS",
            "OBSCONDITIONS",
            "ARIZONA_TIMEOBS",
            "EBV",
            "SPECMODEL_SKY_GFLUX",
            "SPECMODEL_SKY_RFLUX",
            "SPECMODEL_SKY_ZFLUX",
            "GFA_ORIGIN",
            "B_DEPTH",
            "R_DEPTH",
            "Z_DEPTH",
            "B_DEPTH_EBVAIR",
            "R_DEPTH_EBVAIR",
            "Z_DEPTH_EBVAIR",
            "EFFTIME_DARK",
            "EFFTIME_BRIGHT",
            "EFFTIME_BACKUP",
            "SPEED_DARK",
            "SPEED_BRIGHT",
            "SPEED_BACKUP",
            "SPECPROD",
            "SPECPROD_BITPSFFN",
            "SPECPROD_BITFRAMEFN",
            "SPECPROD_BITSKYFN",
            "SPECPROD_BITSFRAMEFN",
            "SPECPROD_BITFLUXCALIBFN",
            "SPECPROD_BITCFRAMEFN",
        ]
        # AR nb of assigned targets
        targkeys = ["N_ASSGN_{}".format(targname) for targname in targnames]
        # AR from SKY MONITOR
        skymonkeys = [
            "NEXP",
            "SKYCAM0_MEAN",
            "SKYCAM0_MEAN_ERR",
            "SKYCAM1_MEAN",
            "SKYCAM1_MEAN_ERR",
            "AVERAGE_MEAN",
            "AVERAGE_MEAN_ERR",
        ]
        # AR from GFA
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
            "FIBER_FRACFLUX_BGS",
            "TRANSPFRAC",
            "MAXCONTRAST",
            "MINCONTRAST",
            "KTERM",
            "FRACFLUX_NOMINAL_POINTSOURCE",
            "FRACFLUX_NOMINAL_ELG",
            "FRACFLUX_NOMINAL_BGS",
            "RADPROF_FWHM_ASEC",
            "FIBERFAC",
            "FIBERFAC_ELG",
        ]
        # AR *1d-quantities* from ephem
        ephkeys = [
            key.upper()
            for key in ephem.table.dtype.names
            if len(ephem.table[key].shape) == 1
        ]
        # AR TSNR2
        tsnr2keys = [
            "TSNR2_{}".format(tracer) for tracer in ["BGS", "LRG", "ELG", "QSO"]
        ]
        # AR all keys
        allkeys = ownkeys
        allkeys += tsnr2keys
        allkeys += hdrkeys
        allkeys += targkeys
        allkeys += ["SKYMON_{}".format(key) for key in skymonkeys]
        allkeys += ["GFA_{}".format(key) for key in gfakeys]
        allkeys += ["EPHEM_{}".format(key) for key in ephkeys]
        allkeys = np.array(allkeys)
        # AR all formats
        allfmts = np.array(["E" for key in allkeys], dtype=object)
        sel = np.in1d(allkeys, ["OBSCONDITIONS", "SKYMON_NEXP"] + targkeys)
        allfmts[sel] = "I"
        sel = np.in1d(
            allkeys,
            [
                "NIGHT",
                "EXPID",
                "TILEID",
                "SPECPROD_BITPSFFN",
                "SPECPROD_BITFRAMEFN",
                "SPECPROD_BITSKYFN",
                "SPECPROD_BITSFRAMEFN",
                "SPECPROD_BITFLUXCALIBFN",
                "SPECPROD_BITCFRAMEFN",
            ],
        )
        allfmts[sel] = "K"
        i = np.where(allkeys == "FIELD")[0][0]
        allfmts[i] = "{}A".format(np.max([len(x) for x in tiles["FIELD"]]))
        i = np.where(allkeys == "TARGETS")[0][0]
        allfmts[i] = "{}A".format(np.max([len(x) for x in tiles["TARGETS"]]))
        i = np.where(allkeys == "GFA_ORIGIN")[0][0]
        allfmts[i] = "{}A".format(np.max([len(x) for x in list(gfa.keys())]))
        i = np.where(allkeys == "ARIZONA_TIMEOBS")[0][0]
        allfmts[i] = "19A"  # AR e.g. 2021-02-18T04:28:06
        i = np.where(allkeys == "SPECPROD")[0][0]
        allfmts[i] = "{}A".format(np.max([len("daily"), len(args.refspecprod)]))
        sel = np.in1d(allkeys, ["EPHEM_{}".format(key) for key in ephkeys])
        allfmts[sel] = "D"
        # AR all units
        allunits = np.array(["" for key in allkeys], dtype=object)
        sel = np.in1d(allkeys, ["OBSCONDITIONS"])
        allunits[sel] = "1=DARK, 2=GRAY, 4=BRIGHT, -1=ELSE"
        sel = np.in1d(
            allkeys,
            [
                "EXPTIME",
                "B_DEPTH",
                "R_DEPTH",
                "Z_DEPTH",
                "B_DEPTH_EBVAIR",
                "R_DEPTH_EBVAIR",
                "Z_DEPTH_EBVAIR",
                "EFFTIME_DARK",
                "EFFTIME_BRIGHT",
                "EFFTIME_BACKUP",
            ],
        )
        allunits[sel] = "seconds"
        sel = np.in1d(
            allkeys,
            ["SPECMODEL_SKY_GFLUX", "SPECMODEL_SKY_RFLUX", "SPECMODEL_SKY_ZFLUX"],
        )
        allunits[sel] = "nMgy / arcsec ** 2"
        sel = np.in1d(
            allkeys, ["SKYMON_{}".format(key) for key in skymonkeys if key != "NEXP"]
        )
        allunits[sel] = "flux units from the database"
        sel = np.in1d(
            allkeys, ["GFA_{}".format(key) for key in gfakeys if "FRAC" in key]
        )
        allunits[sel] = "for an 1.52 arcsecond diameter aperture size"
        # AR dictionary
        exposures = {}
        for key, fmt in zip(allkeys, allfmts):
            if key == "NIGHT":
                exposures[key] = night + np.zeros(nexp, dtype=int)
            elif key == "EXPID":
                exposures[key] = expids.copy()
            elif key == "SPECPROD":
                exposures[key] = np.array([specprod for i in range(nexp)])
            # AR setting depth/efftime/speed to zero by default
            elif key in [
                "B_DEPTH",
                "R_DEPTH",
                "Z_DEPTH",
                "B_DEPTH_EBVAIR",
                "R_DEPTH_EBVAIR",
                "Z_DEPTH_EBVAIR",
                "EFFTIME_DARK",
                "EFFTIME_BRIGHT",
                "EFFTIME_BACKUP",
                "SPEED_DARK",
                "SPEED_BRIGHT",
                "SPEED_BACKUP",
                ]:
                exposures[key] = np.zeros(nexp, dtype=float)
            elif fmt[-1] == "A":
                exposures[key] = np.array(
                    ["-" for i in range(nexp)], dtype="S{}".format(fmt[:-1])
                )
            elif fmt in ["L"]:
                exposures[key] = np.zeros(nexp, dtype=bool)
            elif fmt in ["I", "J", "K"]:
                exposures[key] = -99 + np.zeros(nexp, dtype=int)
            elif fmt in ["E", "D"]:
                exposures[key] = -99 + np.zeros(nexp, dtype=float)
            else:
                sys.exit(
                    "(key,fmt)=({},{}) -> fmt not handled; exiting".format(key, fmt)
                )
        # AR looping on all exposures
        for iexp, expid in enumerate(expids):
            # AR listing existing specprod processed files
            for ftype in ["PSF", "FRAME", "SKY", "SFRAME", "FLUXCALIB", "CFRAME"]:
                value = 0
                for campet_name, campet_bit in zip(campet_names, campet_bits):
                    fn = os.path.join(
                        os.getenv("DESI_ROOT"),
                        "spectro",
                        "redux",
                        specprod,
                        "exposures",
                        "{}".format(night),
                        "{:08d}".format(expid),
                        "{}-{}-{:08d}.fits".format(ftype.lower(), campet_name, expid),
                    )
                    if os.path.isfile(fn):
                        value += 2 ** campet_bit
                exposures["SPECPROD_BIT{}FN".format(ftype)][iexp] = value
            print(
                "\t",
                night,
                expid,
                specprod,
                exposures["SPECPROD_BITSKYFN"][iexp],
                exposures["SPECPROD_BITSFRAMEFN"][iexp],
                exposures["SPECPROD_BITCFRAMEFN"][iexp],
            )
            # AR TSNR2
            # AR TSNR2: first reading all the values
            tsnr2 = {}
            for key in ["BGS", "LRG", "ELG", "QSO"]:
                tsnr2[key] = {
                    campet_name: np.nan + np.zeros(500) for campet_name in campet_names
                }
            for campet_name, campet_bit in zip(campet_names, campet_bits):
                if (exposures["SPECPROD_BITCFRAMEFN"][iexp] & 2 ** campet_bit) > 0:
                    fn = os.path.join(
                        os.getenv("DESI_ROOT"),
                        "spectro",
                        "redux",
                        specprod,
                        "exposures",
                        "{}".format(night),
                        "{:08d}".format(expid),
                        "cframe-{}-{:08d}.fits".format(campet_name, expid),
                    )
                    d = fits.open(fn)["SCORES"].data
                    for key in ["BGS", "LRG", "ELG", "QSO"]:
                        tsnr2[key][campet_name] = d[
                            "TSNR2_{}_{}".format(key, campet_name[0].upper())
                        ]
            # AR TSNR2: for each petal, for each fiber summing the 3 cameras,
            # AR                        then taking the mean over the 500 fibers
            # AR        then taking the mean over the 10 petals
            for key in ["BGS", "LRG", "ELG", "QSO"]:
                tsnr2_allpetals = np.nan + np.zeros(10)
                for i in range(10):
                    tsnr2_allpetals[i] = np.nanmean(
                        tsnr2[key]["b{}".format(i)]
                        + tsnr2[key]["r{}".format(i)]
                        + tsnr2[key]["z{}".format(i)]
                    )
                if np.isfinite(tsnr2_allpetals).sum() > 0:
                    exposures["TSNR2_{}".format(key)][iexp] = np.nanmean(
                        tsnr2_allpetals
                    )
            # AR getting header from the ext=1 of the raw data
            fn = os.path.join(
                rawdir,
                "{}".format(night),
                "{:08d}".format(expid),
                "desi-{:08d}.fits.fz".format(expid),
            )
            hdr = fits.getheader(fn, 1)
            # AR header informations
            for key in hdrkeys:
                if key == "MJDOBS":
                    exposures[key][iexp] = hdr["MJD-OBS"]
                else:
                    exposures[key][iexp] = hdr[key]
            # AR Arizona time of observation
            time_utc = pytz.utc.localize(
                Time(hdr["MJD-OBS"], format="mjd").to_datetime()
            )
            exposures["ARIZONA_TIMEOBS"][iexp] = time_utc.astimezone(
                pytz.timezone("US/Arizona")
            ).strftime("%Y-%m-%dT%H:%M:%S")
            # AR field, targets
            it = np.where(tiles["TILEID"] == hdr["TILEID"])[0][0]
            exposures["FIELD"][iexp] = tiles["FIELD"][it]
            exposures["TARGETS"][iexp] = tiles["TARGETS"][it]
            # AR ebv (excluding sky fibers where ebv=0)
            tmpd = fitsio.read(tiles["FN"][it], columns=["EBV", "OBJTYPE"])
            exposures["EBV"][iexp] = float(
                "{:.3f}".format(np.median(tmpd["EBV"][tmpd["OBJTYPE"] == "TGT"]))
            )
            # AR number of targets per tracer
            for targname, targkey in zip(targnames, targkeys):
                exposures[targkey][iexp] = tiles[targname][it]
            # AR SKY_{GRZ}MAG_AB from integrating the sky model over the decam gzr-bands
            # AR we convert to linear flux units (nMgy/arcsec**2)
            if exposures["SPECPROD_BITSKYFN"][iexp] != 0:
                tmpgmag, tmprmag, tmpzmag = get_sky_grzmag_ab(
                    night, expid, specprod, exposures["EXPTIME"][iexp], "model",
                )
                exposures["SPECMODEL_SKY_GFLUX"][iexp] = 10.0 ** (
                    (22.5 - tmpgmag) / 2.5
                )
                exposures["SPECMODEL_SKY_RFLUX"][iexp] = 10.0 ** (
                    (22.5 - tmprmag) / 2.5
                )
                exposures["SPECMODEL_SKY_ZFLUX"][iexp] = 10.0 ** (
                    (22.5 - tmpzmag) / 2.5
                )
            # AR SKY MONITOR measurement
            # AR we keep those in "native" units (DJS email from 02Mar2021)
            keep = skymon["MJD"] >= exposures["MJDOBS"][iexp]
            keep &= (
                skymon["MJD"]
                <= exposures["MJDOBS"][iexp]
                + exposures["EXPTIME"][iexp] / 3600.0 / 24.0
            )
            keep &= np.isfinite(skymon["SKYCAM0"])
            keep &= np.isfinite(skymon["SKYCAM1"])
            keep &= np.isfinite(skymon["AVERAGE"])
            exposures["SKYMON_NEXP"][iexp] = keep.sum()
            for quant in ["SKYCAM0", "SKYCAM1", "AVERAGE"]:
                if keep.sum() > 0:
                    exposures["SKYMON_{}_MEAN".format(quant)][iexp] = skymon[quant][
                        keep
                    ].mean()
                    exposures["SKYMON_{}_MEAN_ERR".format(quant)][iexp] = skymon[quant][
                        keep
                    ].std() / np.sqrt(keep.sum())
            # AR GFA information
            # AR    first trying matched_coadd
            # AR    if not present, trying acq
            # AR    else no gfa...
            gfa_origin = "matched_coadd"
            ii = np.where(gfa[gfa_origin]["EXPID"] == expid)[0]
            if len(ii) == 0:
                gfa_origin = "acq"
                ii = np.where(gfa[gfa_origin]["EXPID"] == expid)[0]
            if len(ii) > 1:
                sys.exit("More than 1 GFA match: exiting")
            if len(ii) == 1:
                exposures["GFA_ORIGIN"][iexp] = gfa_origin
                i = ii[0]
                for key in gfakeys:
                    # AR TRANSPARENCY x FIBER_FRACFLUX
                    if key == "TRANSPFRAC":
                        exposures["GFA_{}".format(key)][iexp] = (
                            gfa[gfa_origin]["TRANSPARENCY"][i]
                            * gfa[gfa_origin]["FIBER_FRACFLUX"][i]
                        )
                    else:
                        # AR FIBERFAC, FIBERFAC_ELG, FIBERFAC_BGS not in *acq*fits
                        #if (key not in ["FIBERFAC", "FIBERFAC_ELG", "FIBERFAC_BGS"]) or (
                        if (key not in ["FIBERFAC", "FIBERFAC_ELG", "FIBERFAC_BGS", "FRACFLUX_NOMINAL_BGS", "FIBER_FRACFLUX_BGS"]) or (
                            gfa_origin == "matched_coadd"
                        ):
                            exposures["GFA_{}".format(key)][iexp] = gfa[gfa_origin][
                                key
                            ][i]
                # AR/DK exposure depths (needs gfa information)
                # AR/DK adding also depth including ebv+airmass
                if exposures["SPECPROD_BITSKYFN"][iexp] > 0:
                    depths_i = determine_tile_depth2(
                        night,
                        expid,
                        specprod,
                        exposures["EXPTIME"][iexp],
                        exposures["GFA_TRANSPFRAC"][iexp],
                    )
                    for camera in ["B", "R", "Z"]:
                        exposures["{}_DEPTH".format(camera)][iexp] = depths_i[
                            camera.lower()
                        ]
                        ebv = exposures["EBV"][iexp]
                        fact_ebv = 10.0 ** (
                            -2 * 0.4 * depth_coeffs["EBV"][camera] * ebv
                        )
                        airmass = exposures["GFA_AIRMASS"][iexp]
                        fact_air = 10.0 ** (
                            -2 * 0.4 * depth_coeffs["AIRMASS"][camera] * (airmass - 1.0)
                        )
                        exposures["{}_DEPTH_EBVAIR".format(camera)][iexp] = (
                            depths_i[camera.lower()] * fact_ebv * fact_air
                        )
                    # AR/DJS effective exposure times
                    (
                        efftime_dark,
                        efftime_bright,
                        efftime_backup,
                        speed_dark,
                        speed_bright,
                        speed_backup,
                    ) = get_efftime_speed(
                        exposures["EXPTIME"][iexp],
                        # exposures["SKYMON_AVERAGE_MEAN"][iexp],
                        exposures["SPECMODEL_SKY_RFLUX"][iexp],
                        exposures["EBV"][iexp],
                        exposures["GFA_TRANSPARENCY"][iexp],
                        exposures["GFA_AIRMASS"][iexp],
                        exposures["GFA_FIBER_FRACFLUX"][iexp],
                        exposures["GFA_FIBER_FRACFLUX_ELG"][iexp],
                        exposures["GFA_FIBER_FRACFLUX_BGS"][iexp],
                        exposures["GFA_FRACFLUX_NOMINAL_POINTSOURCE"][iexp],
                        exposures["GFA_FRACFLUX_NOMINAL_ELG"][iexp],
                        exposures["GFA_FRACFLUX_NOMINAL_BGS"][iexp],
                        exposures["GFA_KTERM"][iexp],
                    )
                    exposures["EFFTIME_DARK"][iexp] = efftime_dark
                    exposures["EFFTIME_BRIGHT"][iexp] = efftime_bright
                    exposures["EFFTIME_BACKUP"][iexp] = efftime_backup
                    exposures["SPEED_DARK"][iexp] = speed_dark
                    exposures["SPEED_BRIGHT"][iexp] = speed_bright
                    exposures["SPEED_BACKUP"][iexp] = speed_backup
        # AR mjd at the middle of the exposures
        midexp_mjds = exposures["MJDOBS"] + exposures["EXPTIME"] / 2.0 / 3600.0 / 24.0
        # AR ephem 1d-quantities
        nightind = get_nightind(midexp_mjds)
        for key in ephkeys:
            exposures["EPHEM_{}".format(key)] = ephem.table[
                key.lower().replace("lst", "LST")
            ][nightind]
        # AR obsconditions (DARK=1, GRAY=2, BRIGHT=4, else==-1)
        exposures["OBSCONDITIONS"] = get_conditions(midexp_mjds)
        # AR building/writing fits
        cols = []
        for key, fmt, unit in zip(allkeys, allfmts, allunits):
            cols += [fits.Column(name=key, format=fmt, unit=unit, array=exposures[key])]
        h = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
        h.writeto(outfn, overwrite=True)
    return "{}-{}".format(night, len(expids))


# AR per exposure information
if args.exposures == "y":
    # AR listing existing nights
    # AR using daily for all nights
    # AR as it should contain any existing exposures
    nights = np.array(
        [
            int(fn.split("/")[-1])
            for fn in np.sort(
                glob(
                    os.path.join(
                        os.getenv("DESI_ROOT"),
                        "spectro",
                        "redux",
                        "daily",
                        "exposures",
                        "202?????",
                    )
                )
            )
            if int(fn.split("/")[-1]) >= firstnight
        ]
    )
    # nights = [20210101]
    # nights = nights[nights<20201231]
    print(nights)
    # AR update only the ~latest nights?
    # AR TODO: make some sanity check that all relevant nights
    # AR TODO:     not processed has a fits file
    if args.update == "y":
        if not os.path.isfile(outfns["exposures"]):
            print("{} is missing; exiting".format(outfns["exposures"]))
            sys.exit()
        else:
            d = fits.open(outfns["exposures"])[1].data
            firstnight_to_process = d["NIGHT"].max().astype(str)
        nights_to_process = nights[nights >= firstnight_to_process]
    else:
        nights_to_process = nights.copy()
    # AR processing nights
    # AR wrapper on process_night() given an night
    def _process_night(night):
        nightnexp = process_night(night, nightoutdir, skymon, gfa, ephem, rawdir, tiles)
        return nightnexp

    # AR parallel processing?
    if args.numproc > 1:
        pool = sharedmem.MapReduce(np=args.numproc)
        with pool:
            nightnexps = pool.map(_process_night, nights_to_process)
    else:
        nightnexps = []
        for night in nights_to_process:
            nightnexps += [_process_night(night)]
    print(nightnexps)
    # AR merging all nights
    # AR for now, just merging all existing night files
    # AR TODO: make some sanity check to verify that no night is missing
    fns = np.sort(glob("{}/sv1-exposures-202?????.fits".format(nightoutdir)))
    hs = [fits.open(fn) for fn in fns]
    ns = [h[1].data.shape[0] for h in hs]
    hmerge = fits.BinTableHDU.from_columns(hs[0][1].columns, nrows=np.sum(ns))
    start = 0
    for h, n in zip(hs, ns):
        for key in hmerge.columns.names:
            hmerge.data[key][start : start + n] = h[1].data[key]
        start += n
    # AR storing thrufn + ephem infos
    for camera in ["b", "r", "z"]:
        hmerge.header["THRU{}".format(camera)] = thrufn[camera]
    hmerge.header["EPHEMFN"] = ephem_fn
    hmerge.header["DMAXSA"] = ephem_config.programs.DARK.max_sun_altitude().value
    hmerge.header["GMAXMI"] = ephem_config.programs.GRAY.max_moon_illumination()
    hmerge.header[
        "GMAXMIAP"
    ] = ephem_config.programs.GRAY.max_moon_illumination_altitude_product().value
    hmerge.header["BMAXSA"] = ephem_config.programs.BRIGHT.max_sun_altitude().value
    # AR writing to fits
    hmerge.writeto(outfns["exposures"], overwrite=True)


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
    # AR observed exposures
    d = fits.open(outfns["exposures"])[1].data
    d = d[(d["EXPTIME"] > 100) & (d["OBSCONDITIONS"] != -1) ] # AR remove crap
    #
    ramw, decmw = get_radec_mw(tiles["TILERA"], tiles["TILEDEC"], org)
    for radec in np.unique(tiles["radec"]):
        ii = np.where(tiles["radec"] == radec)[0]
        for i in ii:
            # AR is observed?
            if tiles["TILEID"][i] in d["TILEID"]:
                marker, s, alpha = "o", 50, 0.5
            else:
                marker, s, alpha = "+", 20, 1
            ax.scatter(ramw[i], decmw[i], c=tiles["COLOR"][i], marker=marker, s=s, alpha=alpha)
    # AR
    for key in list(flavdict.keys()):
        ax.scatter(100, 100, marker="o", s=50, c=flavdict[key]["COLOR"], label=key)
    ax.legend(loc=3, ncol=2)
    plt.savefig(outfns["skymap"], bbox_inches="tight")
    plt.close()

    # AR per-night on-sky exptime
    # AR  and per-tile efftime_dark
    d = fits.open(outfns["exposures"])[1].data
    # AR remove CRAP...
    d = d[(d["EXPTIME"] > 100) & (d["OBSCONDITIONS"] != -1)]
    nights, ii = np.unique(d["NIGHT"], return_index=True)
    mjds = 0.5 + d["EPHEM_NOON"][ii]
    #
    names = ["DARK", "GRAY", "BRIGHT"]
    vals = [1, 2, 4]
    cols = ["k", "g", "y"]
    #
    # AR per night
    #
    # AR getting stats
    exptimes = np.zeros((len(nights), 4))
    depths = np.zeros((len(nights), 4))
    for i in range(len(nights)):
        for j in range(3):
            keep = (d["NIGHT"] == nights[i]) & (d["OBSCONDITIONS"] == vals[j])
            exptimes[i, j] = d["EXPTIME"][keep].sum() / 3600.0
            depths[i, j] = d["R_DEPTH_EBVAIR"][keep].sum() / 3600.0
    # AR
    fig, ax = plt.subplots(figsize=(20, 5))
    for i in range(len(nights)):
        start = 0
        if i == 0:
            labels = [
                "{}={:.0f} hrs".format(names[j], exptimes[:, j].sum()) for j in range(3)
            ]
        else:
            labels = [None, None, None]
        for j in range(3):
            tmpx = mjds[i] + np.array([-0.25, 0.25, 0.25, -0.25])
            tmpy = start + np.array([0, 0, exptimes[i, j], exptimes[i, j]])
            ax.fill(tmpx, tmpy, fill=True, color=cols[j], alpha=0.5, label=labels[j])
            start += exptimes[i, j]
            if j == 0:
                ax.text(
                    mjds[i],
                    0,
                    "{}".format(nights[i]),
                    rotation=45,
                    ha="right",
                    va="top",
                )
    # AR
    ax.grid(True)
    ax.set_axisbelow(True)
    ax.set_title(
        "SV1 on-sky EXPTIME from {} to {} ({} exposures with EXPTIME>100 and OBSCONDITIONS!=-1)".format(
            d["NIGHT"].min(), d["NIGHT"].max(), len(d),
        )
    )
    ax.set_xlabel("MJD-OBS")
    ax.set_ylabel("EXPTIME per NIGHT [hours]")
    ax.set_ylim(-2, 10)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.legend(loc=2)
    plt.savefig(outfns["onsky"], bbox_inches="tight")
    plt.close()
    #
    # AR per tile
    #
    for targets in ["QSO+LRG", "ELG", "QSO+ELG", "BGS+MWS"]:
        if targets == "BGS+MWS":
            effkey = "EFFTIME_BRIGHT"
            expnom = 150
            ylim = (-0.25 * 3600, 1 * 3600)
            ytxt = -100
        else:
            effkey = "EFFTIME_DARK"
            expnom = 1000
            ylim = (-1 * 3600, 5 * 3600)
            ytxt = -500
        tileids = np.unique(d["TILEID"][d["TARGETS"] == targets])
        # AR getting stats
        efftimes = np.zeros((len(tileids), 4))
        for i in range(len(tileids)):
            for j in range(3):
                keep = (d["TILEID"] == tileids[i]) & (d["OBSCONDITIONS"] == vals[j])
                efftimes[i, j] = d[effkey][keep].sum()
        # AR
        fig, ax = plt.subplots(figsize=(20, 5))
        for i in range(len(tileids)):
            start = 0
            if i == 0:
                labels = [
                    "{}={:.0f} hrs cumulated".format(
                        names[j], efftimes[:, j].sum() / 3600.0
                    )
                    for j in range(3)
                ]
            else:
                labels = [None, None, None]
            for j in range(3):
                tmpx = i + np.array([-0.25, 0.25, 0.25, -0.25])
                tmpy = start + np.array([0, 0, efftimes[i, j], efftimes[i, j]])
                ax.fill(
                    tmpx, tmpy, fill=True, color=cols[j], alpha=0.5, label=labels[j]
                )
                start += efftimes[i, j]
                if j == 0:
                    ax.text(
                        i,
                        ytxt,
                        "{}".format(tileids[i]),
                        rotation=45,
                        ha="right",
                        va="top",
                    )
        # AR
        ax.axhline(
            4.0 * expnom,
            color="r",
            ls="--",
            label="SV goal 4x{}s nominal".format(expnom),
        )
        ax.grid(True)
        ax.set_axisbelow(True)
        ax.set_title(
            "{} : {} from {} to {} ({} exposures with EXPTIME>100 and OBSCONDITIONS!=-1)".format(
                targets,
                effkey,
                d["NIGHT"].min(),
                d["NIGHT"].max(),
                (d["TARGETS"] == targets).sum(),
            )
        )
        ax.set_xlabel("")
        ax.set_ylabel("{} per TILEID [s]".format(effkey))
        ax.set_ylim(ylim)
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_locator(MultipleLocator(1000))
        ax.set_yticklabels(
            ["{:.0f}".format(x) if x > 0 else "" for x in ax.get_yticks()]
        )
        ax.legend(loc=1)
        plt.savefig(
            outfns["efftime-{}".format(targets.lower().replace("+", ""))],
            bbox_inches="tight",
        )
        plt.close()

    # AR observing conditions : cumulative distributions
    # targetss = ["BGS+MWS", "QSO+LRG", "ELG", "QSO+ELG"]
    # cols = ["g", "r", "b", "c"]
    d = fits.open(outfns["exposures"])[1].data
    d = d[d["GFA_ORIGIN"] != "-"]
    keys = [
        "EBV",
        "GFA_AIRMASS",
        "GFA_TRANSPARENCY",
        "GFA_FWHM_ASEC",
        "GFA_SKY_MAG_AB",
        "GFA_FIBER_FRACFLUX",
        "GFA_FIBER_FRACFLUX_ELG",
    ]
    xmins = [0, 1.0, 0, 0.5, 17, 0, 0]
    xmaxs = [0.2, 2.5, 1.05, 3.5, 22, 0.8, 0.8]
    nx, ny = 2, 4
    fig = plt.figure(figsize=(20, 10))
    gs = gridspec.GridSpec(nx, ny, wspace=0.25, hspace=0.20)
    ip = 0
    for key, xmin, xmax in zip(keys, xmins, xmaxs):
        ax = plt.subplot(gs[ip])
        bins = np.linspace(xmin, xmax, 101)
        # for targets, col in zip(["ALL"] + targetss, ["k"] + cols):
        for targ, col in zip(
            ["ALL", "BGS+MWS", "QSO+LRG", "ELG", "QSO+ELG"], ["k", "g", "r", "b", "c"]
        ):
            if targ == "ALL":
                keep = np.ones(len(d), dtype=bool)
            else:
                keep = d["TARGETS"] == targ
            _ = ax.hist(
                d[key][keep],
                bins=bins,
                cumulative=True,
                density=True,
                histtype="step",
                color=col,
                alpha=0.8,
                label="{} (median={:.2f})".format(targ, np.median(d[key][keep])),
            )
        ax.set_xlabel(key)
        if ip % ny == 0:
            ax.set_ylabel("Cumulative fraction")
        ax.grid(True)
        ax.set_axisbelow(True)
        if key in ["EBV", "GFA_AIRMASS", "GFA_FWHM_ASEC"]:
            ax.legend(loc=4)
        else:
            ax.legend(loc=2)
        ax.axhline(0.5, ls="--", c="k")
        ip += 1
    plt.savefig(outfns["obscond"]["cumul"], bbox_inches="tight")
    plt.close()

    # AR observing conditions = f(MJD) (one plot per month)
    keys = [
        "GFA_MOON_ILLUMINATION",
        "GFA_MOON_ZD_DEG",
        "GFA_MOON_SEP_DEG",
        "GFA_AIRMASS",
        "GFA_TRANSPARENCY",
        "GFA_FWHM_ASEC",
        "GFA_SKY_MAG_AB",
        "GFA_FIBER_FRACFLUX",
        "EFFTIME_DARK / EXPTIME",
    ]
    mlocs = [0.20, 25, 25, 0.20, 0.20, 0.50, 1.0, 0.20, 0.5]
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
                if keys[i] == "EFFTIME_DARK / EXPTIME":
                    y = d["EFFTIME_DARK"] / d["EXPTIME"]
                else:
                    y = d["{}".format(keys[i])]
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

    # AR exptime, depth, efftime
    # targetss = ["BGS+MWS", "QSO+LRG", "ELG", "QSO+ELG"]
    # cols = ["g", "r", "b", "c"]
    d = fits.open(outfns["exposures"])[1].data
    d = d[(d["EXPTIME"] > 100) & (d["OBSCONDITIONS"] != -1) & (d["EFFTIME_DARK"] > 0)]
    title = "SV1 observations from {} to {}\n({} exposures with EXPTIME>100 and OBSCONDITIONS!=-1 and EFFTIME_DARK>0)".format(
        d["NIGHT"].min(), d["NIGHT"].max(), len(d),
    )
    keys = ["EXPTIME", "EFFTIME_DARK"]
    xlabels = ["EXPTIME [s]", "EFFTIME_DARK [s]"]
    ylabels = [
        "Cumulative number of tiles",
        "Cumulative number of tiles",
    ]
    xlim = (0, 25000)
    fig = plt.figure(figsize=(10, 5))
    gs = gridspec.GridSpec(2, 1, hspace=0.3)
    ax = {
        "EXPTIME": plt.subplot(gs[0]),
        "EFFTIME_DARK": plt.subplot(gs[1]),
    }
    # for targets, col in zip(targetss, cols):
    for targ, col in zip(
        ["BGS+MWS", "QSO+LRG", "ELG", "QSO+ELG"], ["orange", "r", "b", "g"]
    ):
        mydict = {}
        mydict["TILEID"] = np.unique(d["TILEID"][d["TARGETS"] == targ])
        ntiles = len(mydict["TILEID"])
        for key in keys:
            mydict[key] = np.zeros(ntiles)
        for i in range(ntiles):
            keep = d["TILEID"] == mydict["TILEID"][i]
            for key in keys:
                mydict[key][i] = int(d[key][keep].sum())
        for key in keys:
            bins = np.linspace(
                0, 1 + mydict[key].max(), int(mydict[key].max() / 100.0)
            ).astype(
                int
            )  # 1 bin per 100 seconds...
            _ = ax[key].hist(
                mydict[key],
                bins=bins,
                cumulative=True,
                color=col,
                alpha=0.3,
                label="{} ({} tiles, median={:.0f}s)".format(
                    targ, ntiles, np.median(mydict[key])
                ),
            )
    for key, xlabel, ylabel in zip(keys, xlabels, ylabels):
        if key == "EXPTIME":
            ax[key].set_title(title)
        ax[key].set_xlabel(xlabel)
        ax[key].set_ylabel(ylabel)
        ax[key].set_xlim(0, 25000)
        ax[key].xaxis.set_major_locator(MultipleLocator(2000))
        ax[key].yaxis.set_major_locator(MultipleLocator(5))
        ax[key].grid(True)
        ax[key].set_axisbelow(True)
        ax[key].legend()
    plt.savefig(outfns["depth"]["cumul"], bbox_inches="tight")
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
    htmlmain.write("\t\t<li><a href='#onsky' >Per-night on-sky EXPTIME</a></li>\n")
    htmlmain.write("\t\t<li><a href='#depth' >Tiles: exptime and depth</a></li>\n")
    htmlmain.write(
        "\t\t<li><a href='#tile-nexp-design' >Tiles: NEXP and design</a></li>\n"
    )
    htmlmain.write(
        "\t\t<li><a href='#per-exposure-properties' > Per-exposure properties</a> ({} exposure(s) over {} night(s))</li>\n".format(
            len(d), len(np.unique(d["NIGHT"]))
        )
    )
    htmlmain.write(
        "\t\t<li><a href='#obsconds-cumul' >Observing conditions: cumulative distributions</a></li>\n"
    )
    for month in months:
        if os.path.isfile(outfns["obscond"][month]):
            htmlmain.write(
                "\t\t<li><a href='#obsconds-{}' > Observing conditions for {}</a></li>\n".format(
                    month, month
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
        "<p style='font-size:1vw'><a href='../sv1-exposures.fits' target='external'> sv1-exposures.fits </a> : main file with all exposures. Column content: {}.</p>\n".format(
            ", ".join(d.dtype.names),
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

    # AR Per-night on-sky EXPTIME
    htmlmain.write("<h2><a id='onsky' href='#onsky' > Per-night on-sky EXPTIME </a>\n")
    htmlmain.write(
        "<a href='#top' style='position: absolute; right: 0;'>Top of the page</a></h2>\n"
    )
    tmppng = outfns["onsky"].replace(args.outdir, "../")
    htmlmain.write("<tr>\n")
    htmlmain.write(
        "<td align=center><a href='{}'><img SRC='{}' width=80% height=auto></a></td>\n".format(
            tmppng, tmppng
        )
    )
    htmlmain.write("</tr>\n")
    htmlmain.write("\n")

    # AR Tiles: exptime and depth
    htmlmain.write("<h2><a id='depth' href='#depth' > Tiles: exptime and depth </a>\n")
    htmlmain.write(
        "<a href='#top' style='position: absolute; right: 0;'>Top of the page</a></h2>\n"
    )
    for png in [outfns["depth"]["cumul"]] + [
        outfns["efftime-{}".format(fsh)]
        for fsh in ["qsolrg", "qsoelg", "elg", "bgsmws"]
    ]:
        tmppng = png.replace(args.outdir, "../")
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

    # AR Observing conditions: cumulative distributions
    htmlmain.write(
        "<h2><a id='obsconds-cumul' href='#obsconds-cumul' > Observing conditions: cumulative distributions</a>\n"
    )
    htmlmain.write(
        "<a href='#top' style='position: absolute; right: 0;'>Top of the page</a></h2>\n"
    )
    tmppng = outfns["obscond"]["cumul"].replace(args.outdir, "../")
    htmlmain.write("<tr>\n")
    htmlmain.write(
        "<td align=center><a href='{}'><img SRC='{}' width=50% height=auto></a></td>\n".format(
            tmppng, tmppng
        )
    )
    htmlmain.write("</tr>\n")
    htmlmain.write("\n")

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
