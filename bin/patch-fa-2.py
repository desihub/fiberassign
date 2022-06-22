import os
import fitsio
import numpy as np
from glob import glob
from collections import Counter
from astropy.table import Table
import sys
import tempfile

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

    log('Reading', infn)
    F = fitsio.FITS(infn)
    primhdr = F[0].read_header()
    tileid = primhdr['TILEID']
    fa_surv = primhdr['FA_SURV']
    survey = primhdr.get('SURVEY')
    fa_prog = primhdr.get('FAPRGRM')
    stats.update(infn=infn, tileid=tileid, fa_surv=fa_surv, survey=survey, fa_prog=fa_prog)

    # Find targeting files
    tilestr = '%06i'%tileid

    if fa_surv == 'cmx':
        log('Skipping CMX tile')
        if log_to_string:
            print(logstr)
        return stats,None

    if survey == 'special':
        if fa_prog.startswith('tertiary'):
            targdir_pat = os.path.join(desiroot, 'survey', 'fiberassign',
                                       survey, 'tertiary', '*', '*')
        else:
            targdir_pat = os.path.join(desiroot, 'survey', 'fiberassign',
                                       survey, '*')
    elif fa_surv in ['sv1', 'sv2', 'sv3']:
        targdir_pat = os.path.join(desiroot, 'survey', 'fiberassign',
                                   fa_surv.upper(), '*')
    else:
        targdir_pat = os.path.join(desiroot, 'survey', 'fiberassign',
                                   fa_surv, tilestr[:3])

    targfns = []
    for tag in ['targ', 'sky', 'scnd', 'too']:
        pat = os.path.join(targdir_pat, '%s-%s.fits' % (tilestr, tag))
        fns = glob(pat)
        assert((len(fns) == 0) or (len(fns) == 1))
        if len(fns) == 0:
            stats.update({'exists_'+tag: False})
            log('Does not exist:', pat)
        else:
            stats.update({'exists_'+tag: True})
            fn = fns[0]
            targfns.append(fn)

    patched_tables = {}
    all_patched_data = []

    # Read the original FIBERASSIGN table.
    for hduname in ['FIBERASSIGN', 'TARGETS']:
        log('Reading FA hdu', hduname)
        ff = F[hduname]
        tab = ff.read()

        # We're going to match based on (positive) TARGETID.
        targetid = tab['TARGETID']
        # Targetids we're looking for
        tidset = set(targetid[targetid > 0])

        patched_col_rows = {}

        # This is the subset of columns we patch, per Anand's suggestion
        fixcols = {'FIBERASSIGN':
                   ['SV3_DESI_TARGET', 'SV3_SCND_TARGET', 'MASKBITS',
                    'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG',
                    'FLUX_IVAR_G', 'FLUX_IVAR_Z', 'FLUX_IVAR_R',
                    'REF_CAT', 'REF_ID',
                    'RELEASE', 'BRICKID', 'BRICK_OBJID',
                    'MORPHTYPE', 'BRICKNAME', 'HPXPIXEL',
                    'FLUX_Z', 'FLUX_R', 'FLUX_G',
                    'PRIORITY_INIT',
                   ],
                   'TARGETS':
                   ['SV3_DESI_TARGET', 'SV3_SCND_TARGET',
                    'BRICKID', 'BRICKNAME', 'BRICK_OBJID', 'RELEASE', 'HPXPIXEL']}

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

            destcols = set(tab.dtype.names)
            srccols = set(T.dtype.names)

            for col in fixcols[hduname]:
                #log('Column', col)
                if not col in destcols:
                    log('Column', col, '- Not in FA table; skipping')
                    continue
                old = tab[col][J]


                if not col in srccols:
                    log('Column', col, '- Not in target table; treating as zeros?')
                    new = np.zeros_like(old)
                    zeroed = True
                    #continue
                else:
                    new = T[col][I]
                    zeroed = False
                #log('Old:', len(old), old[:5])
                #log('New:', len(new), new[:5])
                eq = arrays_equal(old, new)
                if not np.all(eq):
                    diff = np.flatnonzero(np.logical_not(eq))
                    log('Column', col, '- patching', len(diff), 'rows')
                    log('  rows', J[diff][:5])
                    log('  old vals', tab[col][J[diff]][:5])
                    log('  new vals', new[diff][:5])
                    tab[col][J[diff]] = new[diff]
                    #patched_rows.update(J[diff])

                    if col in patched_col_rows:
                        already_patched = patched_col_rows[col]
                        for i in J[diff]:
                            if i in already_patched:
                                log('Already patched col', col, 'row', i)
                                sys.exit(-1)
                    patched_col_rows[col] = patched_col_rows.get(col, set()).union(J[diff])

                    all_patched_data.append((infn, hduname, tileid, fa_surv, fn, col,
                                             targetid[J[diff]], old[diff], new[diff], zeroed))
                else:
                    log('Column', col, '- All equal!')

        patched_tables[hduname] = tab
        for col in fixcols[hduname]:
            hdunamemap = {'FIBERASSIGN': 'fa',
                          'TARGETS': 'targ'}
            stats.update({'patched_%s_%s' % (hdunamemap[hduname], col.lower()):
                          len(patched_col_rows.get(col, []))})

    if len(all_patched_data) > 0:
        # Make sure output directory exists
        outdir = os.path.dirname(outfn)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        assert(outfn != infn)
    
        # Write out output file, leaving other HDUs unchanged.
        f,tempout = tempfile.mkstemp(dir=outdir, suffix='.fits')
        os.close(f)
        Fout = fitsio.FITS(tempout, 'rw', clobber=True)
        for ext in F:
            extname = ext.get_extname()
            hdr = ext.read_header()
            data = ext.read()
            if extname == 'PRIMARY':
                # fitsio will add its own headers about the FITS format, so trim out all COMMENT cards.
                newhdr = fitsio.FITSHDR()
                for r in hdr.records():
                    if r['name'] == 'COMMENT':
                        continue
                    newhdr.add_record(r)
                hdr = newhdr
            if extname in patched_tables:
                # Swap in our updated FIBERASSIGN table!
                data = patched_tables[extname]
            Fout.write(data, header=hdr, extname=extname)
        Fout.close()
        os.rename(tempout, outfn)
        log('Wrote', outfn)
    else:
        log('No changes for tile', tileid, 'file', infn)

    if log_to_string:
        print(logstr)
    return stats, all_patched_data

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
        R = mp.map(_bounce_patch, args)
        allstats    = [r1 for r1,r2 in R]
        all_patched = [r2 for r1,r2 in R]
    else:
        allstats = []
        all_patched = []
        for kw in args:
            s,p = patch(**kw)
            allstats.append(s)
            all_patched.append(p)
            print()

    allstats = [s for s in allstats if s is not None]
    all_patched = [p for p in all_patched if p is not None]

    # Collect and write out all patched data.
    infns = []
    hdunames = []
    tileids = []
    fasurvs = []
    fns = []
    cols = []
    targetids = []
    oldvals = []
    newvals = []
    zeroed = []
    for P in all_patched:
        for infn,hdu,tileid,fasurv,fn,col,tid,oldv,newv,zero in P:
            N = len(oldv)
            infns.append([infn]*N)
            hdunames.append([hdu]*N)
            tileids.append([tileid]*N)
            fasurvs.append([fasurv]*N)
            fns.append([fn]*N)
            cols.append([col]*N)
            targetids.append(tid)
            oldvals.append(oldv.astype('U10'))
            newvals.append(newv.astype('U10'))
            zeroed.append([zero]*N)
    array_list=[np.hstack(x) for x in [tileids, targetids, hdunames, fns, cols,
                                       newvals, oldvals, zeroed]]
    names=['TILEID', 'TARGETID', 'EXTENSION', 'ORIGFN', 'KEY', 'ORIGVAL', 'FAVAL', 'ZEROED']
    fitsio.write(os.path.join(opt.outdir, 'all-patched-rows.fits'), array_list, names=names,
                 clobber=True)
    
    # Collate and write out stats.
    allcols = set()
    for s in allstats:
        allcols.update(s.keys())
    allcols = list(allcols)
    allcols.sort()
    print('allcols:', allcols)
    data = dict([(c,[]) for c in allcols])
    for s in allstats:
        for c in allcols:
            data[c].append(s.get(c, 0))
    for c in allcols:
        data[c] = np.array(data[c])
    fitsio.write(opt.stats, data, names=allcols, clobber=True)

def _bounce_patch(x):
    return patch(**x)

if __name__ == '__main__':
    main()
