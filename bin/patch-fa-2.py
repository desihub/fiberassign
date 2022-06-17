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

# We use this function to test for equality between arrays -- both
# being non-finite (eg NaN) counts as equal.
def arrays_equal(a, b):
    bothnan = False
    try:
        bothnan = np.logical_not(np.isfinite(a)) * np.logical_not(np.isfinite(b))
    except:
        pass
    return np.logical_or(a == b, bothnan)

def patch(infn = 'desi/target/fiberassign/tiles/trunk/001/fiberassign-001000.fits.gz',
          desiroot = '/global/cfs/cdirs/desi',
          outfn = None,
          save_logs = False,
          ):
    '''Repair the data corruption reported in https://github.com/desihub/fiberassign/pull/432
    for a given single fiberassign file.

    The patching proceeds by reading the FIBERASSIGN and TARGETS hdus,
    and uses the intermediate fiberassign files to pull in correct
    values.
    
    FIBERASSIGN and TARGETS extensions: check against the TILEID-{targ,sky,scnd,too}.fits files;
    SKY_MONITOR extension: check against the TILEID-sky.fits files;
    GFA_TARGETS extension: check against the TILEID-gfa.fits files;

    where those files are in, eg, $DESI_ROOT/survey/fiberassign/main/001/001099-targ.fits

    Targets are matched on TARGETID, and values updated as necessary.
    
    If any rows are changed, then an updated file is written out,
    preserving all other HDUs.
    '''
    log_to_string = save_logs
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
    tileid = primhdr['TILEID']
    fa_surv = primhdr['FA_SURV']
    stats.update(infn=infn, tileid=tileid, fa_surv=fa_surv)

    # Find targeting files
    tilestr = '%06i'%tileid
    targdir = os.path.join(desiroot, 'survey', 'fiberassign',
                           fa_surv, tilestr[:3])
    targfns = [os.path.join(targdir, '%s-%s.fits' % (tilestr, tag))
               for tag in ['targ', 'sky', 'scnd', 'too']]

    keepfns = []
    for fn in targfns:
        if not os.path.exists(fn):
            log('Does not exist:', fn)
        else:
            keepfns.append(fn)
    targfns = keepfns

    #stats.update(scndfn=scndfn, toofn=toofn)
    #for i,t in enumerate(targdirs):
    #    stats['targdir%i' % i] = t

    # Read the original FIBERASSIGN table.
    for hduname in ['FIBERASSIGN', 'TARGETS']:
        log('Reading FA hdu', hduname)
        ff = F[hduname]
        tab = ff.read()

        # We're going to match based on (positive) TARGETID.
        targetid = tab['TARGETID']
        #ra = tab['TARGET_RA']
        #dec = tab['TARGET_DEC']
        # Per https://github.com/desihub/fiberassign/issues/385,
        # Swap in unique values for negative TARGETIDs.
        # I = np.flatnonzero(targetid < 0)
        # stats.update(neg_targetids=len(I), pos_targetids=len(np.flatnonzero(targetid>0)))
        # if len(I):
        #     tab['TARGETID'][I] = -(tileid * 10000 + tab['LOCATION'][I])
        #     log('Updated', len(I), 'negative TARGETID values')
        #     patched_neg_tids = True

        # Targetids we're looking for
        tidset = set(targetid[targetid > 0])

        # Read in the intermediate fiberassign files for targets
        for fn in targfns:
            T = fitsio.read(fn)
            log('Read', len(T), 'from', fn)
            if len(T) == 0:
                continue
            tids = T['TARGETID']
            I = np.flatnonzero([t in tidset for t in tids])
            log('Found', len(I), 'matching TARGETIDs')
            if len(I) == 0:
                continue
            T = T[I]

            # "targetid" are the ones from the FA file
            # "tids" are the ones from the fiberassign targets files
            tidmap = dict([(tid,i) for i,tid in enumerate(T['TARGETID']) if tid>0])
            I = np.array([tidmap.get(t, -1) for t in targetid])
            J = np.flatnonzero(I >= 0)
            I = I[J]
            log('Checking', len(I), 'matched TARGETIDs')

            tabcols = set(tab.dtype.names)
            
            # This is the subset of columns we patch, based on the bug report.
            # for col in ['BRICK_OBJID', 'BRICKID', 'RELEASE', 'FLUX_G', 'FLUX_R', 'FLUX_Z',
            #             'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z',
            #             'REF_CAT', 'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG',
            #             'MASKBITS', 'REF_ID', 'MORPHTYPE']:
            #for col in T.columns:
            for col in T.dtype.names:
                log('Column', col)
                if not col in tabcols:
                    log('Not in FA table; skipping')
                    continue
                old = tab[col][J]
                new = T[col][I]
                eq = arrays_equal(old, new)
                if not np.all(eq):
                    diff = np.flatnonzero(np.logical_not(eq))
                    log('Target catalogs: col', col, 'patching', len(diff), 'rows')
                    log('  rows', J[diff][:5])
                    log('  old vals', tab[col][J[diff]][:5])
                    log('  new vals', new[diff][:5])
                    tab[col][J[diff]] = new[diff]
                    patched_rows.update(J[diff])

    #stats.update(patched_rows_targets=len(patched_rows))
    # if len(patched_rows) == 0 and not patched_neg_tids:
    #     log('No need to patch', infn)
    #     stats.update(patched=0)
    # else:
    #     log('Patched', len(patched_rows), 'data rows')
    #     stats.update(patched=1)
    # 
    #     # Make sure output directory exists
    #     outdir = os.path.dirname(outfn)
    #     if not os.path.exists(outdir):
    #         os.makedirs(outdir)
    #     assert(outfn != infn)
    # 
    #     # Write out output file, leaving other HDUs unchanged.
    #     Fout = fitsio.FITS(outfn, 'rw', clobber=True)
    #     for ext in F:
    #         #log(ext.get_extname())
    #         extname = ext.get_extname()
    #         hdr = ext.read_header()
    #         data = ext.read()
    #         if extname == 'PRIMARY':
    #             # fitsio will add its own headers about the FITS format, so trim out all COMMENT cards.
    #             newhdr = fitsio.FITSHDR()
    #             for r in hdr.records():
    #                 if r['name'] == 'COMMENT':
    #                     #log('Skipping comment:', r['name'], r['value'])
    #                     continue
    #                 newhdr.add_record(r)
    #             hdr = newhdr
    #         if extname == hduname:
    #             # Swap in our updated FIBERASSIGN table!
    #             data = tab
    #         Fout.write(data, header=hdr, extname=extname)
    #     Fout.close()
    #     log('Wrote', outfn)
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
    save_logs = False
    if opt.threads:
        save_logs = True
    for infn in opt.files:
        #infn = os.path.join(fa_base, '001/fiberassign-001000.fits.gz')
        if not infn.startswith(fa_base):
            print('All input filenames must start with --in-base = ', fa_base)
            return -1
        outfn = infn.replace(fa_base, opt.outdir+'/')
        args.append(dict(infn=infn, outfn=outfn, save_logs=save_logs))

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
