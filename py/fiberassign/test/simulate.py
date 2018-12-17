"""
Simulation utilities for fiberassign tests.
"""
import os

import numpy as np

import fitsio

from desitarget.targetmask import desi_mask

from fiberassign import __version__

from fiberassign.hardware import load_hardware

from fiberassign.targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                                 TARGET_TYPE_STANDARD)


def sim_data_dir():
    dir = "test_fiberassign_output"
    if not os.path.isdir(dir):
        os.makedirs(dir)
    return dir


def sim_tiles(path):
    tile_dtype = np.dtype([
        ("TILEID", "i4"),
        ("RA", "f8"),
        ("DEC", "f8"),
        ("IN_DESI", "i4"),
        ("PROGRAM", "S6"),
        ("OBSCONDITIONS", "i4")
    ])
    ntile = 7
    fdata = np.zeros(ntile, dtype=tile_dtype)
    fdata[0] = (1165, 150.69, 33.86, 1, "DARK", 1)
    fdata[1] = (6927, 151.78, 33.84, 1, "DARK", 1)
    fdata[2] = (11108, 150.87, 31.23, 1, "DARK", 1)
    fdata[3] = (16870, 151.96, 31.21, 1, "DARK", 1)
    fdata[4] = (18465, 150.47, 33.20, 1, "DARK", 1)
    fdata[5] = (24227, 151.56, 33.18, 2, "DARK", 1)
    fdata[6] = (28408, 150.73, 30.52, 2, "DARK", 1)

    if os.path.isfile(path):
        os.remove(path)
    fd = fitsio.FITS(path, "rw")

    header = dict()
    header["FBAVER"] = __version__
    fd.write(fdata, header=header)
    return


def sim_targets(path, tgtype, tgoffset, density=5000.0):
    ramin = 148.0
    ramax = 154.0
    decmin = 28.0
    decmax = 37.0
    target_dtype = np.dtype([
        ("TARGETID", "i8"),
        ("RA", "f8"),
        ("DEC", "f8"),
        ("DESI_TARGET", "i8"),
        ("PRIORITY", "i4"),
        ("SUBPRIORITY", "f8"),
        ("OBSCONDITIONS", "i4"),
        ("NUMOBS_MORE", "i4")
    ])
    ndim = np.sqrt(density)
    nra = int(ndim * (ramax - ramin))
    ndec = int(ndim * (decmax - decmin))
    ntarget = nra * ndec

    fdata = np.zeros(ntarget, dtype=target_dtype)
    fdata["TARGETID"] = tgoffset + np.arange(ntarget)
    fdata["RA"] = np.random.uniform(low=ramin, high=ramax, size=ntarget)
    fdata["DEC"] = np.random.uniform(low=decmin, high=decmax, size=ntarget)
    fdata["OBSCONDITIONS"] = np.ones(ntarget, dtype=np.int32)
    fdata["NUMOBS_MORE"] = np.ones(ntarget, dtype=np.int32)
    fdata["SUBPRIORITY"] = np.random.uniform(low=0.0, high=1.0, size=ntarget)

    sky_mask = desi_mask["SKY"].mask
    std_mask = desi_mask["STD_BRIGHT"].mask
    # We could be fancier and set the DESI_TARGET bits for science targets
    # to exactly match the priority class.  For now we just set a single
    # bit that will all the fiberassign code to determine which targets are
    # science targets.
    sci_mask = desi_mask["ELG"].mask

    if tgtype == TARGET_TYPE_SKY:
        fdata["PRIORITY"] = np.zeros(ntarget, dtype=np.int32)
        fdata["DESI_TARGET"] |= sky_mask
    elif tgtype == TARGET_TYPE_STANDARD:
        fdata["PRIORITY"] = 1500 * np.ones(ntarget, dtype=np.int32)
        fdata["DESI_TARGET"] |= std_mask
    else:
        fdata["DESI_TARGET"] |= sci_mask
        # These are the fractions of each target type:
        dist = {
            1600: 0.001,
            2000: 0.160,
            2100: 0.190,
            2998: 0.001,
            3000: 0.440,
            3200: 0.090,
            3400: 0.118
        }
        ndist = {x: int(y*ntarget) for x, y in dist.items()}
        ntot = np.sum([y for x, y in ndist.items()])
        if ntot < ntarget:
            ndist[3400] += ntarget - ntot
        priority = np.concatenate([(x * np.ones(y, dtype=np.int32))
                                   for x, y in ndist.items()])
        np.random.shuffle(priority)
        fdata["PRIORITY"] = priority

    if os.path.isfile(path):
        os.remove(path)
    fd = fitsio.FITS(path, "rw")

    header = dict()
    header["FBAVER"] = __version__
    fd.write(fdata, header=header)
    return ntarget
