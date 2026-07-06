"""
Simulation utilities for fiberassign tests.
"""
import os
import shutil
from datetime import datetime, timezone

from collections import OrderedDict

import numpy as np

import fitsio

from desitarget.targetmask import desi_mask

import desimodel.io as dmio

from fiberassign import __version__

from fiberassign.hardware import (load_hardware, FIBER_STATE_STUCK,
                                  FIBER_STATE_BROKEN)

from fiberassign.targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                                 TARGET_TYPE_SUPPSKY,
                                 TARGET_TYPE_SUPPSKY, TARGET_TYPE_STANDARD)

sim_assign_date = "2020-01-01T00:00:00+00:00"

def sim_data_dir():
    dir = "test_fiberassign_output"
    if not os.path.isdir(dir):
        os.makedirs(dir)
    return dir


def sim_data_subdir_create(name):
    test_dir = os.path.join(sim_data_dir(), name)
    if os.path.isdir(test_dir):
        shutil.rmtree(test_dir)
    os.makedirs(test_dir)
    return test_dir


def sim_science_fractions():
    """
    Elements of the tuple are (priority, numobs, fraction, desitarget_type).
    The fractions should add up to 1.0...
    """
    return [
        (3000, 1, 0.200, "ELG"),
        (3200, 1, 0.066, "LRG"),
        (3210, 1, 0.034, "LGE"),
        (3400, 1, 0.500, "QSO"),
        (3400, 3, 0.200, "QSO")
    ]


def sim_focalplane(rundate=None, fakepos=False):
    runtime = None
    if rundate is None:
        runtime = datetime.now(tz=timezone.utc)
    else:
        runtime = datetime.strptime(rundate, "%Y-%m-%dT%H:%M:%S%z")

    # First get the starting focalplane from desimodel
    fp, exclude, state, tmstr = dmio.load_focalplane(runtime)

    # Make a copy since desimodel is caching the original
    fp = fp.copy()
    state = state.copy()

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

def sim_stuck_sky(dirnm, hw, tiles):
    from desiutil.brick import Bricks
    bricks = Bricks(bricksize=1)
    tab = bricks.to_table()
    skybricks_dir = os.path.join(dirnm, 'skybricks')
    os.makedirs(skybricks_dir)
    skyfn = os.path.join(skybricks_dir, 'skybricks-exist.fits')
    tab.write(skyfn)

    # create bricks touched by 'tiles'.
    rad = hw.focalplane_radius_deg

    def radec2xyz(r, d):
        r = np.deg2rad(r)
        d = np.deg2rad(d)
        x = np.cos(r) * np.cos(d)
        y = np.sin(r) * np.cos(d)
        z = np.sin(d)
        return np.vstack((x,y,z)).T

    skyxyz = radec2xyz(tab['RA'], tab['DEC'])
    ns = skyxyz.shape[0]
    tilexyz = radec2xyz(tiles.ra, tiles.dec)
    nt = tilexyz.shape[0]
    maxr2 = (np.deg2rad(rad + np.sqrt(2)/2.))**2

    skybricks = set()
    for i in range(nt):
        r2 = np.sum((tilexyz[i,:] - skyxyz)**2, axis=1)
        J = np.flatnonzero(r2 < maxr2)
        skybricks.update(J)
    for i in skybricks:
        print('Writing fake SKYBRICK', tab['BRICKNAME'][i])
        # skybrick size in pixels
        sz = 3672
        # fraction of pixels to say are good sky
        goodsky_frac = 0.75
        randsky = (16*(np.random.uniform(size=(sz,sz)) > goodsky_frac)).astype(np.uint8)
        # Add WCS header
        hdr = fitsio.FITSHDR()
        hdr['CTYPE1'] = 'RA---TAN'
        hdr['CTYPE2'] = 'DEC--TAN'
        hdr['CRVAL1'] = tab['RA'][i]
        hdr['CRVAL2'] = tab['DEC'][i]
        hdr['CRPIX1'] = (sz+1)/2.
        hdr['CRPIX2'] = (sz+1)/2.
        hdr['CD1_1'] = -1./3600.
        hdr['CD1_2'] =  0.
        hdr['CD2_1'] =  0.
        hdr['CD2_2'] =  1./3600.
        outfn = os.path.join(skybricks_dir, 'sky-%s.fits.gz' % tab['BRICKNAME'][i])
        fitsio.write(outfn, randsky, header=hdr)

    os.environ['SKYBRICKS_DIR'] = skybricks_dir

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
    ntile = 3
    fdata = np.zeros(ntile, dtype=tile_dtype)
    fdata[0] = (1234, 150.0, 31.0, 1, "DARK", 1)
    fdata[1] = (1235, 150.2, 31.0, 1, "DARK", 1)
    fdata[2] = (1236, 150.1, 31.2, 1, "DARK", 1)

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


def sim_targets(path, tgtype, tgoffset, density=5000.0, science_frac=None):
    ramin = 147.0
    ramax = 153.0
    decmin = 28.0
    decmax = 34.0
    target_cols = OrderedDict([
        ("TARGETID", "i8"),
        ("RA", "f8"),
        ("DEC", "f8"),
        ("DESI_TARGET", "i8"),
        ("BRICKID", "i4"),
        ("BRICK_OBJID", "i4"),
        ("BRICKNAME", "S8"),
        ("PRIORITY", "i4"),
        ("SUBPRIORITY", "f8"),
        ("OBSCONDITIONS", "i4"),
        ("NUMOBS_INIT", "i4"),
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
    fdata["SUBPRIORITY"] = np.random.uniform(low=0.0, high=1.0, size=ntarget)

    sky_mask = desi_mask["SKY"].mask
    suppsky_mask = desi_mask["SUPP_SKY"].mask
    std_mask = desi_mask["STD_BRIGHT"].mask

    if tgtype == TARGET_TYPE_SKY:
        fdata["PRIORITY"] = np.zeros(ntarget, dtype=np.int32)
        fdata["DESI_TARGET"] |= sky_mask
    elif tgtype == TARGET_TYPE_STANDARD:
        fdata["PRIORITY"] = 1500 * np.ones(ntarget, dtype=np.int32)
        fdata["DESI_TARGET"] |= std_mask
    elif tgtype == TARGET_TYPE_SUPPSKY:
        fdata["PRIORITY"] = np.zeros(ntarget, dtype=np.int32)
        fdata["DESI_TARGET"] |= suppsky_mask
    else:
        # These are the fractions of each target type:
        dist = science_frac
        if dist is None:
            dist = sim_science_fractions()
        ntot = 0
        for priority, numobs, frac, tgbits in dist:
            ntg = int(frac * ntarget)
            fdata["PRIORITY"][ntot:ntot+ntg] = priority
            fdata["NUMOBS_INIT"][ntot:ntot+ntg] = numobs
            fdata["DESI_TARGET"][ntot:ntot+ntg] = desi_mask[tgbits].mask
            ntot += ntg
        if ntot < ntarget:
            # Add extra targets of the final type to make up the difference
            fdata["PRIORITY"][ntot:ntarget] = dist[-1][0]
            fdata["NUMOBS_INIT"][ntot:ntarget] = dist[-1][1]
            fdata["DESI_TARGET"][ntot:ntarget] = desi_mask[dist[-1][3]].mask

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
