import os
import fitsio
import numpy as np
from glob import glob
from collections import Counter
from astropy.table import Table
import sys

sys.path.append('desicode')
from desimodel.footprint import radec2pix

_sec_cache = {}
def read_secondary(fn):
    global _sec_cache
    try:
        return _sec_cache[fn]
    except KeyError:
        pass
    f = fitsio.read(fn)
    _sec_cache[fn] = f
    return f

def patch(infn = 'desi/target/fiberassign/tiles/trunk/001/fiberassign-001000.fits.gz',
          desiroot = '/global/cfs/cdirs/desi',
          outfn = None,
          ):
    #'desi/target/fiberassign/tiles/trunk/000/fiberassign-000001.fits.gz'
    #'desi/target/fiberassign/tiles/trunk/020/fiberassign-020004.fits.gz'
    patched_rows = set()

    print('Reading', infn)
    F = fitsio.FITS(infn)
    primhdr = F[0].read_header()
    tilera  = primhdr['TILERA']
    tiledec = primhdr['TILEDEC']

    mtldir = primhdr['MTL']
    targdir = primhdr['TARG']
    scndfn = primhdr['SCND']
    toofn = primhdr.get('TOO', '')

    mtldir = mtldir.replace('DESIROOT', desiroot)
    targdir = targdir.replace('DESIROOT', desiroot)
    scndfn = scndfn.replace('DESIROOT', desiroot)
    toofn = toofn.replace('DESIROOT', desiroot)

    hduname = 'FIBERASSIGN'
    ff = F[hduname]
    tab = ff.read()

    targetid = tab['TARGETID']
    ra = tab['TARGET_RA']
    dec = tab['TARGET_DEC']

    # Some objects have NaN TARGET_RA; skip these when looking up the healpix
    Iok = np.flatnonzero(np.isfinite(ra))
    nside = 8
    hps = radec2pix(nside, ra[Iok], dec[Iok])

    alltargets = []
    tidset = set(targetid)
    for hp in np.unique(hps):
        pat = os.path.join(targdir, '*-hp-%03i.fits' % hp)
        #print(pat)
        fns = glob(pat)
        #print(fns)
        assert(len(fns) == 1)
        fn = fns[0]
        T = fitsio.read(fn)
        tid = T['TARGETID']
        I = np.flatnonzero([t in tidset for t in tid])
        print('Read targets', fn, '->', len(T), 'targets,', len(I), 'matching TARGETIDs')
        if len(I) == 0:
            continue
        alltargets.append(T[I])    

    targets = np.hstack(alltargets)
    tidmap = dict([(tid,i) for i,tid in enumerate(targets['TARGETID'])])
    I = np.array([tidmap.get(t, -1) for t in targetid])
    J = np.flatnonzero(I >= 0)
    I = I[J]

    def equal(a, b):
        bothnan = False
        try:
            bothnan = np.logical_not(np.isfinite(a)) * np.logical_not(np.isfinite(b))
        except:
            pass
        return np.logical_or(a == b, bothnan)

    for col in ['BRICK_OBJID', 'BRICKID', 'RELEASE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z',
                'REF_CAT', 'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG', 'MASKBITS', 'REF_ID',
                'MORPHTYPE']:
        old = tab[col][J]
        new = targets[col][I]
        eq = equal(old, new)
        #print('Target catalogs:', col, 'equal:', Counter(eq))
        if not np.all(eq):
            diff = np.flatnonzero(old != new)
            #print(old[diff], 'vs', new[diff])
            print('Target catalogs: col', col, 'patching', len(diff), 'rows')
            tab[col][J[diff]] = new[diff]
            patched_rows.update(J[diff])

    scnd = read_secondary(scndfn)
    print('Read', len(scnd), 'secondaries from', scndfn)
    sidmap = dict([(tid,i) for i,tid in enumerate(scnd['TARGETID'])])

    I = np.array([sidmap.get(t, -1) for t in targetid])
    J = np.flatnonzero(I >= 0)
    I = I[J]

    for col in ['FLUX_G', 'FLUX_R', 'FLUX_Z',
                'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG',]:
        old = tab[col][J]
        new = scnd[col][I]
        eq = equal(old, new)
        #print(col, Counter(eq))
        if not np.all(eq):
            diff = np.flatnonzero(eq)
            print('Secondary targets: col', col, 'patching', len(diff), 'rows')
            tab[col][J[diff]] = new[diff]
            #print(old[diff], 'vs', new[diff])
            #print('Difference:', old[diff] - new[diff])
            patched_rows.update(J[diff])

    toomap = {
        '/global/cfs/cdirs/desi/survey/ops/staging/mtl/main/ToO/ToO.ecsv' :
        '/global/cfs/cdirs/desi/target/ToO/ToO.ecsv',

        '/data/afternoon_planning/surveyops/trunk/mtl/sv3/ToO/ToO.ecsv' :
        '/global/cfs/cdirs/desi/target/ToO/sv3/ToO.ecsv',
        }
    toofn = toomap.get(toofn, toofn)

    if len(toofn):
        too = Table.read(toofn)
        print('Read', len(too), 'ToO from', toofn)
        toomap = dict([(tid,i) for i,tid in enumerate(too['TARGETID'])])

        I = np.array([toomap.get(t, -1) for t in targetid])
        J = np.flatnonzero(I >= 0)
        I = I[J]

        for col in ['FLUX_G', 'FLUX_R', 'FLUX_Z',
                    'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG',]:
            old = tab[col][J]
            new = too[col][I]
            eq = equal(old, new)
            if not np.all(eq):
                diff = np.flatnonzero(eq)
                print('ToO targets: col', col, 'patching', len(diff), 'rows')
                tab[col][J[diff]] = new[diff]
                patched_rows.update(J[diff])

    if len(patched_rows) == 0:
        print('No need to patch', infn)
        return 0
    print('Patched', len(patched_rows), 'rows')

    outdir = os.path.dirname(outfn)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    assert(outfn != infn)

    if os.path.exists(outfn):
        os.remove(outfn)
    Fout = fitsio.FITS(outfn, 'rw')
    for ext in F:
        #print(ext.get_extname())
        extname = ext.get_extname()
        hdr = ext.read_header()
        data = ext.read()
        if extname == 'PRIMARY':
            newhdr = fitsio.FITSHDR()
            for r in hdr.records():
                if r['name'] == 'COMMENT':
                    #print('Skipping comment:', r['name'], r['value'])
                    continue
                newhdr.add_record(r)
            hdr = newhdr
        if extname == hduname:
            data = tab
        Fout.write(data, header=hdr, extname=extname)
    Fout.close()
    print('Wrote', outfn)
    return len(patched_rows)

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('files', type=str, nargs='*', help='Fiberassign files to patch')
    parser.add_argument('--in-base', default='/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk',
                        help='Input base path')
    parser.add_argument('--outdir', default='patched-fa')
    opt = parser.parse_args()

    fa_base = opt.in_base
    for infn in opt.files:
        #infn = os.path.join(fa_base, '001/fiberassign-001000.fits.gz')
        if not infn.startswith(fa_base):
            print('All input filenames must start with --in-base = ', fa_base)
            return -1
        outfn = infn.replace(fa_base, 'patched-fa/')
        patch(infn=infn, outfn=outfn)
        print()

if __name__ == '__main__':
    main()
