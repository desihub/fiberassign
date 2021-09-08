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
    '''Read the given secondary targets file (FITS), with caching.'''
    global _sec_cache
    try:
        return _sec_cache[fn]
    except KeyError:
        pass
    scnd = fitsio.read(fn)
    sidmap = dict([(tid,i) for i,tid in enumerate(scnd['TARGETID']) if tid > 0])
    val = scnd,sidmap
    # Only cache one!!
    _sec_cache = {fn: val}
    #_sec_cache[fn] = val
    return val

def patch(infn = 'desi/target/fiberassign/tiles/trunk/001/fiberassign-001000.fits.gz',
          desiroot = '/global/cfs/cdirs/desi',
          outfn = None,
          ):
    '''Repair the data corruption reported in https://github.com/desihub/fiberassign/pull/375
    for a given single fiberassign file.

    The patching proceeds by reading the FIBERASSIGN hdu, while checking the headers
    for the locations of targeting files used.

    Based on TARGETID, we match objects in the FIBERASSIGN file and
    the targeting files.  For a subset of columns, we check for cases
    where the targeting files contain different values, and patch them
    in to the FIBERASSIGN table.

    First we read the primary target files themselves, then secondary
    targets, then ToO targets.

    If any rows are changed, then an updated file is written out, preserving all other HDUs.

    '''
    log_to_string = True
    logstr = ''
    def log(*a):
        nonlocal logstr
        if log_to_string:
            logstr += ' '.join([str(s) for s in a]) + '\n'
        else:
            print(*a)

    stats = {}
    patched_rows = set()
    patched_neg_tids = False

    log('Reading', infn)
    F = fitsio.FITS(infn)
    primhdr = F[0].read_header()
    #tilera  = primhdr['TILERA']
    #tiledec = primhdr['TILEDEC']
    tileid = primhdr['TILEID']
    fa_surv = primhdr['FA_SURV']
    stats.update(infn=infn, tileid=tileid, fa_surv=fa_surv)

    # Find targeting files / directories in the headers.
    #mtldir = primhdr['MTL']
    targdirs = [primhdr['TARG']]
    for i in range(2, 100):
        key = 'TARG%i' % i
        if not key in primhdr:
            break
        targdirs.append(primhdr[key])
    scndfn = primhdr.get('SCND', '').strip()
    toofn = primhdr.get('TOO', '')

    # Swap in the local file locations. (Header values may start with the string "DESIROOT")
    #mtldir = mtldir.replace('DESIROOT', desiroot)
    targdirs = [targdir.replace('DESIROOT', desiroot)
                for targdir in targdirs]
    scndfn = scndfn.replace('DESIROOT', desiroot)
    toofn = toofn.replace('DESIROOT', desiroot)

    # Tile 80611 has /data as the prefix.
    pre = '/data/'
    targdirs = [targdir.replace(pre, desiroot+'/')
                if targdir.startswith(pre)
                else targdir
                for targdir in targdirs]
    if scndfn.startswith(pre):
        scndfn = scndfn.replace(pre, desiroot+'/')
    if toofn.startswith(pre):
        toofn = toofn.replace(pre, desiroot+'/')

    # Tile 80687 has SCND a directory, not a file.
    if os.path.isdir(scndfn):
        # Tile 80687 has SCND a directory that contains two files, with the same number of rows
        # but one includes DR9 photometry columns.
        fn = os.path.join(scndfn, 'sv1targets-dark-secondary-dr9photometry.fits')
        if os.path.exists(fn):
            log('Special-case updated SCND from', scndfn, 'to', fn)
            scndfn = fn
        # Similar with 80717.
        fn = os.path.join(scndfn, 'sv1targets-bright-secondary-dr9photometry.fits')
        if os.path.exists(fn):
            log('Special-case updated SCND from', scndfn, 'to', fn)
            scndfn = fn

    stats.update(scndfn=scndfn, toofn=toofn)
    for i,t in enumerate(targdirs):
        stats['targdir%i' % i] = t

    # Read the original FIBERASSIGN table.
    hduname = 'FIBERASSIGN'
    ff = F[hduname]
    tab = ff.read()

    # We're going to match based on (positive) TARGETID.
    targetid = tab['TARGETID']
    ra = tab['TARGET_RA']
    dec = tab['TARGET_DEC']
    # Per https://github.com/desihub/fiberassign/issues/385,
    # Swap in unique values for negative TARGETIDs.
    I = np.flatnonzero(targetid < 0)
    stats.update(neg_targetids=len(I), pos_targetids=len(np.flatnonzero(targetid>0)))
    if len(I):
        tab['TARGETID'][I] = -(tileid * 10000 + tab['LOCATION'][I])
        log('Updated', len(I), 'negative TARGETID values')
        patched_neg_tids = True

    # Start with primary target files.  These are split into healpixes, so look up which healpixes
    # contain targets used in this tile.

    # Some objects have NaN TARGET_RA; skip these when looking up the healpix
    Iok = np.flatnonzero(np.isfinite(ra))
    nside = 8
    hps = radec2pix(nside, ra[Iok], dec[Iok])

    # Read relevant healpixes.
    alltargets = []
    tidset = set(targetid[targetid > 0])
    for hp in np.unique(hps):
        foundhp = False
        for targdir in targdirs:
            pat = os.path.join(targdir, '*-hp-%i.fits' % hp)
            #log(pat)
            fns = glob(pat)
            #log(fns)
            if len(fns) != 1:
                log('Searching for pattern', pat, 'found files', fns, '; expected 1 file.')
            assert(len(fns) <= 1)
            if len(fns) == 0:
                continue
            foundhp = True
            fn = fns[0]
            T = fitsio.read(fn)
            tid = T['TARGETID']
            I = np.flatnonzero([t in tidset for t in tid])
            log('Read targets', fn, '->', len(T), 'targets,', len(I), 'matching (positive) TARGETIDs')
            if len(I) == 0:
                continue
            alltargets.append(T[I])
        assert(foundhp)

    # Merge targets and match on TARGETID.
    targets = np.hstack(alltargets)
    tidmap = dict([(tid,i) for i,tid in enumerate(targets['TARGETID']) if tid>0])
    I = np.array([tidmap.get(t, -1) for t in targetid])
    J = np.flatnonzero(I >= 0)
    I = I[J]
    log('Checking', len(I), 'matched TARGETIDs')
    stats.update(matched_targets=len(I))

    # We use this function to test for equality between arrays -- both
    # being non-finite (eg NaN) counts as equal.
    def equal(a, b):
        bothnan = False
        try:
            bothnan = np.logical_not(np.isfinite(a)) * np.logical_not(np.isfinite(b))
        except:
            pass
        return np.logical_or(a == b, bothnan)

    # This is the subset of columns we patch, based on the bug report.
    for col in ['BRICK_OBJID', 'BRICKID', 'RELEASE', 'FLUX_G', 'FLUX_R', 'FLUX_Z',
                'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z',
                'REF_CAT', 'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG',
                'MASKBITS', 'REF_ID', 'MORPHTYPE']:
        old = tab[col][J]
        new = targets[col][I]
        eq = equal(old, new)
        #log('Target catalogs:', col, 'equal:', Counter(eq))
        if not np.all(eq):
            diff = np.flatnonzero(np.logical_not(eq))
            log('Target catalogs: col', col, 'patching', len(diff), 'rows')
            log('  rows', J[diff][:5])
            log('  old vals', tab[col][J[diff]][:5])
            log('  new vals', new[diff][:5])
            tab[col][J[diff]] = new[diff]
            patched_rows.update(J[diff])

    stats.update(patched_rows_targets=len(patched_rows))

    # Were secondary targets used for this tile?
    if scndfn == '-':
        # eg fiberassign-081000.fits.gz
        scndfn = ''
    if len(scndfn):
        log('Reading secondary targets file:', scndfn)
        scnd,sidmap = read_secondary(scndfn)
        log('Read', len(scnd), 'secondaries from', scndfn)

        # Match on targetid
        I = np.array([sidmap.get(t, -1) for t in targetid])
        J = np.flatnonzero(I >= 0)
        I = I[J]
        log('Matched', len(I), 'secondary targets on TARGETID')
        stats.update(matched_scnd=len(I))

        for col in ['FLUX_G', 'FLUX_R', 'FLUX_Z',
                    'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG',]:
            # if not col in tab.dtype.fields:
            #     log('No column', col, 'in original FIBERASSIGN table')
            #     continue
            if not col in scnd.dtype.fields:
                log('No column', col, 'in secondary table')
                continue
            old = tab[col][J]
            new = scnd[col][I]
            eq = equal(old, new)
            #log(col, Counter(eq))
            if not np.all(eq):
                diff = np.flatnonzero(np.logical_not(eq))
                log('Secondary targets: col', col, 'patching', len(diff), 'rows')
                log('  rows', J[diff][:5])
                log('  old vals', tab[col][J[diff]][:5])
                log('  new vals', new[diff][:5])
                tab[col][J[diff]] = new[diff]
                patched_rows.update(J[diff])

        stats.update(patched_rows_scnd=len(patched_rows))

    # ToO files -- we're going to ignore the header values and look
    # up the file location based on FA_SURV.
    toomap = {
        'sv3': '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/sv3/ToO/ToO.ecsv',
        'main': '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/main/ToO/ToO.ecsv',
        'cmx': '',
        'sv2': '',
        'sv1': '',
    }
    if not fa_surv in toomap:
        log('Failed to find ToO file for FA_SURV', fa_surv)
    toofn = toomap[fa_surv]

    # Were ToO targets used in this tile?
    if len(toofn):
        log('Reading ToO file', toofn)
        too = Table.read(toofn)
        log('Read', len(too), 'ToO from', toofn)
        toomap = dict([(tid,i) for i,tid in enumerate(too['TARGETID']) if tid > 0])

        I = np.array([toomap.get(t, -1) for t in targetid])
        J = np.flatnonzero(I >= 0)
        I = I[J]
        log('Matched', len(I), 'ToO targets on TARGETID')
        stats.update(matched_too=len(I))

        for col in ['FLUX_G', 'FLUX_R', 'FLUX_Z',
                    'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG',]:
            old = tab[col][J]
            new = too[col][I]
            eq = equal(old, new)
            if not np.all(eq):
                diff = np.flatnonzero(np.logical_not(eq))
                log('ToO targets: col', col, 'patching', len(diff), 'rows')
                tab[col][J[diff]] = new[diff]
                patched_rows.update(J[diff])

        stats.update(patched_rows_too=len(patched_rows))

    if len(patched_rows) == 0 and not patched_neg_tids:
        log('No need to patch', infn)
        stats.update(patched=0)
    else:
        log('Patched', len(patched_rows), 'data rows')
        stats.update(patched=1)

        # Make sure output directory exists
        outdir = os.path.dirname(outfn)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        assert(outfn != infn)

        # Write out output file, leaving other HDUs unchanged.
        Fout = fitsio.FITS(outfn, 'rw', clobber=True)
        for ext in F:
            #log(ext.get_extname())
            extname = ext.get_extname()
            hdr = ext.read_header()
            data = ext.read()
            if extname == 'PRIMARY':
                # fitsio will add its own headers about the FITS format, so trim out all COMMENT cards.
                newhdr = fitsio.FITSHDR()
                for r in hdr.records():
                    if r['name'] == 'COMMENT':
                        #log('Skipping comment:', r['name'], r['value'])
                        continue
                    newhdr.add_record(r)
                hdr = newhdr
            if extname == hduname:
                # Swap in our updated FIBERASSIGN table!
                data = tab
            Fout.write(data, header=hdr, extname=extname)
        Fout.close()
        log('Wrote', outfn)
    if log_to_string:
        print(logstr)
    return stats

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('files', type=str, nargs='*', help='Fiberassign files to patch')
    parser.add_argument('--in-base', default='/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk',
                        help='Input base path')
    parser.add_argument('--outdir', default='patched-fa')
    parser.add_argument('--stats', default='patch-fa-stats.fits',
                        help='Output filename for patching statistics.')
    parser.add_argument('--threads', type=int, help='Multiprocessing threads')
    opt = parser.parse_args()

    args = []
    fa_base = opt.in_base
    for infn in opt.files:
        #infn = os.path.join(fa_base, '001/fiberassign-001000.fits.gz')
        if not infn.startswith(fa_base):
            print('All input filenames must start with --in-base = ', fa_base)
            return -1
        outfn = infn.replace(fa_base, 'patched-fa/')
        args.append(dict(infn=infn, outfn=outfn))

    if opt.threads:
        from astrometry.util.multiproc import multiproc
        mp = multiproc(opt.threads)
        allstats = mp.map(_bounce_patch, args)
    else:
        allstats = []
        for kw in args:
            stats = patch(**kw)
            allstats.append(stats)
            print()

    # Collate and write out stats.
    allcols = set()
    for s in allstats:
        allcols.update(s.keys())
    allcols = list(allcols)
    allcols.sort()
    print('allcols:', allcols)
    data = dict([(c,[]) for c in allcols])
    # if not ''
    blankvals = dict(
        tileid=-1,
        neg_targetids=-1,
        pos_targetids=-1,
        matched_targets=-1,
        patched_rows_targets=-1,
        matched_scnd=-1,
        patched_rows_scnd=-1,
        matched_too=-1,
        patched_rows_too=-1,
        patched=-1,
        )
    for s in allstats:
        for c in allcols:
            try:
                data[c].append(s[c])
            except KeyError:
                data[c].append(blankvals.get(c, ''))
    for c in allcols:
        #print('Column', c, 'length', len(data[c]))
        data[c] = np.array(data[c])
        #print('Np array:', len(data[c]))
        #print(data[c])
    fitsio.write(opt.stats, data, names=allcols, clobber=True)

def _bounce_patch(x):
    return patch(**x)

if __name__ == '__main__':
    main()
