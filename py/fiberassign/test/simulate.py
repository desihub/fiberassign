"""
Simulation utilities for fiberassign tests.
"""
import os
import shutil
from datetime import datetime

from collections import OrderedDict

import numpy as np

import fitsio

from desitarget.targetmask import desi_mask

import desimodel.io as dmio

from fiberassign import __version__

from fiberassign.hardware import (load_hardware, FIBER_STATE_STUCK,
                                  FIBER_STATE_BROKEN)

from fiberassign.targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                                 TARGET_TYPE_STANDARD)


def sim_data_dir():
    dir = "test_fiberassign_output"
    if not os.path.isdir(dir):
        os.makedirs(dir)
    return dir


def test_subdir_create(name):
    test_dir = os.path.join(sim_data_dir(), name)
    if os.path.isdir(test_dir):
        shutil.rmtree(test_dir)
    os.makedirs(test_dir)
    return test_dir


def sim_focalplane(runtime=None, fakepos=False):
    if runtime is None:
        runtime = datetime.utcnow()

    # First get the starting focalplane from desimodel
    fp, exclude, state, tmstr = dmio.load_focalplane(runtime)

    npos = len(fp)

    # Are we replacing the arms and theta lengths with the nominal values?
    # This is useful in the field rotation test.
    if fakepos:
        phys_t = 380.0
        phys_p = 200.0
        t_min = -0.5 * phys_t
        t_max = 0.5 * phys_t
        p_min = 185.0 - phys_p
        p_max = 185.0
        for row in range(npos):
            petalrot_deg = (float(7 + fp["PETAL"][row]) * 36.0) % 360.0
            fp["OFFSET_T"][row] = -170.0 + petalrot_deg
            fp["OFFSET_P"][row] = -5.0
            fp["MIN_T"][row] = t_min
            fp["MAX_T"][row] = t_max
            fp["MIN_P"][row] = p_min
            fp["MAX_P"][row] = p_max
            fp["LENGTH_R1"][row] = 3.0
            fp["LENGTH_R2"][row] = 3.0
            state["STATE"][row] = 0
    else:
        # Now set some fibers to stuck / broken
        mods = {
            95: FIBER_STATE_BROKEN,
            62: FIBER_STATE_BROKEN,
            102: FIBER_STATE_STUCK,
            82: FIBER_STATE_STUCK,
            131: FIBER_STATE_STUCK
        }
        for row, loc in enumerate(state["LOCATION"]):
            if loc in mods:
                state["STATE"][row] = mods[loc]

    return fp, exclude, state


def sim_tiles(path, selectfile=None):
    tile_dtype = np.dtype([
        ("TILEID", "i4"),
        ("RA", "f8"),
        ("DEC", "f8"),
        ("IN_DESI", "i4"),
        ("PROGRAM", "S6"),
        ("OBSCONDITIONS", "i4")
    ])
    # ntile = 7
    ntile = 2
    fdata = np.zeros(ntile, dtype=tile_dtype)
    fdata[0] = (11108, 150.87, 31.23, 1, "DARK", 1)
    fdata[1] = (1165, 150.69, 33.86, 1, "DARK", 1)
    # fdata[2] = (18465, 150.47, 33.20, 1, "DARK", 1)
    # fdata[3] = (6927, 151.78, 33.84, 1, "DARK", 1)
    # fdata[4] = (24227, 151.56, 33.18, 2, "DARK", 1)
    # fdata[5] = (16870, 151.96, 31.21, 1, "DARK", 1)
    # fdata[6] = (28408, 150.73, 30.52, 2, "DARK", 1)

    if os.path.isfile(path):
        os.remove(path)
    fd = fitsio.FITS(path, "rw")

    header = dict()
    header["FBAVER"] = __version__
    fd.write(fdata, header=header)

    if selectfile is not None:
        if os.path.isfile(selectfile):
            os.remove(selectfile)
        tindx = 0
        select = list()
        for td in fdata:
            if tindx % 2 == 0:
                select.append(td[0])
        with open(selectfile, "w") as f:
            for t in select:
                f.write("{}\n".format(t))
    return


def sim_targets(path, tgtype, tgoffset, density=5000.0):
    ramin = 148.0
    ramax = 154.0
    decmin = 28.0
    decmax = 37.0
    target_cols = OrderedDict([
        ("TARGETID", "i8"),
        ("RA", "f8"),
        ("DEC", "f8"),
        ("DESI_TARGET", "i8"),
        ("BRICKID", "i4"),
        ("BRICK_OBJID", "i4"),
        ("BRICKNAME", "a8"),
        ("PRIORITY", "i4"),
        ("SUBPRIORITY", "f8"),
        ("OBSCONDITIONS", "i4"),
        ("NUMOBS_MORE", "i4"),
        ("FIBERFLUX_G", "f4"),
        ("FIBERFLUX_R", "f4"),
        ("FIBERFLUX_Z", "f4"),
        ("FIBERFLUX_IVAR_G", "f4"),
        ("FIBERFLUX_IVAR_R", "f4"),
        ("FIBERFLUX_IVAR_Z", "f4")
    ])

    target_dtype = np.dtype([(x, y) for x, y in target_cols.items()])

    ndim = np.sqrt(density)
    nra = int(ndim * (ramax - ramin))
    ndec = int(ndim * (decmax - decmin))
    ntarget = nra * ndec

    fdata = np.zeros(ntarget, dtype=target_dtype)
    fdata["TARGETID"] = tgoffset + np.arange(ntarget)
    fdata["RA"] = np.random.uniform(low=ramin, high=ramax, size=ntarget)
    fdata["DEC"] = np.random.uniform(low=decmin, high=decmax, size=ntarget)
    fdata["OBSCONDITIONS"] = np.ones(ntarget, dtype=np.int32)
    fdata["NUMOBS_MORE"] = 5 * np.ones(ntarget, dtype=np.int32)
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


def petal_rotation(npos, reverse=False):
    """Return a dictionary that implements the petal rotation.
    """
    rot_petal = dict()
    for p in range(10):
        if reverse:
            newp = p - npos
            if newp < 0:
                newp += 10
        else:
            newp = p + npos
            if newp > 9:
                newp -= 10
        rot_petal[p] = newp
    return rot_petal
