'''
Functions for working with DESI mocks and fiberassignment

TODO (maybe):
This contains hardcoded hacks, especially wrt priorities and
interpretation of object types
'''

from __future__ import print_function, division

import sys, os
import numpy as np
from astropy.table import Table, Column
from fiberassign import io
from desitarget import desi_mask
import desitarget
import desispec.brick

def load_rdzipn(infile):
    """Read rdzipn infile and return target and truth tables
    """
    ra, dec, z, itype, priority, numobs = io.read_rdzipn(infile)
    n = len(ra)

    #- Martin's itype is 1 to n, while Bob's is 0 to n-1
    itype -= 1
    assert np.min(itype >= 0)

    #- rdzipn has float32 ra, dec, but it should be float64
    ra = ra.astype('float64') % 360     #- enforce 0 <= ra < 360
    dec = dec.astype('float64')

    #- MTL z is float64, which is probably overkill but ok
    z = ra.astype('float64')

    #- Hardcoded in rdzipn format
    # 0 : 'QSO',      #- QSO-LyA
    # 1 : 'QSO',      #- QSO-Tracer
    # 2 : 'LRG',      #- LRG
    # 3 : 'ELG',      #- ELG
    # 4 : 'STAR',     #- QSO-Fake
    # 5 : 'UNKNOWN',  #- LRG-Fake
    # 6 : 'STAR',     #- StdStar
    # 7 : 'SKY',      #- Sky
    
    qso_lya    = (itype==0)
    qso_tracer = (itype==1)
    qso_fake   = (itype==4)
    qso = qso_lya | qso_tracer | qso_fake
    lrg_real = (itype==2)
    lrg_fake = (itype==5)
    lrg      = lrg_real | lrg_fake
    elg      = (itype==3)
    std      = (itype==6)
    sky      = (itype==7)

    if not np.any(std):
        print("WARNING: no standard stars found")
    if not np.any(sky):
        print("WARNING: no sky locations found")
    if not np.any(~(std | sky)):
        print("WARNING: no science targets found")

    #- Priorities before we have observed the target
    prio_pre = np.zeros(n, dtype='i4')
    prio_pre[qso]  = desi_mask.QSO.priorities['UNOBS']
    prio_pre[lrg]  = desi_mask.LRG.priorities['UNOBS']
    prio_pre[elg]  = desi_mask.ELG.priorities['UNOBS']
    prio_pre[std] = -1
    prio_pre[sky]  = -1
    assert np.all(prio_pre != 0)

    #- Priorities after we know it is
    prio_post = -99 * np.ones(n, dtype='i4')
    prio_post[qso_lya] = desi_mask.QSO.priorities['MORE_ZGOOD']
    prio_post[qso_tracer] = desi_mask.QSO.priorities['DONE']
    prio_post[qso_fake] = desi_mask.QSO.priorities['DONE']
    prio_post[lrg_real] = desi_mask.LRG.priorities['MORE_ZGOOD']
    prio_post[lrg_fake] = desi_mask.LRG.priorities['DONE']
    prio_post[elg]  = desi_mask.ELG.priorities['DONE']
    prio_post[std] = -1
    prio_post[sky]  = -1
    assert np.all(prio_post != -99)

    #- Hard coded NUMOBS
    numobs_pre = np.zeros(n, dtype='i4')
    numobs_pre[qso] = 4
    numobs_pre[lrg] = 2
    numobs_pre[elg] = 1
    numobs_pre[std] = -1
    numobs_pre[sky]  = -1
    assert np.all(numobs_pre != 0)

    numobs_post = 999 * np.ones(n, dtype='i4')
    numobs_post[qso_lya] = 4
    numobs_post[qso_tracer] = 1
    numobs_post[qso_fake] = 0
    numobs_post[lrg_real] = 2
    numobs_post[lrg_fake] = 0
    numobs_post[elg]  = 1
    numobs_post[std] = -1
    numobs_post[sky]  = -1
    assert np.all(numobs_post != 999)
    
    #- Create a DESI_TARGET mask
    desi_target = np.zeros(n, dtype='i8')
    desi_target[qso] |= desi_mask.QSO
    desi_target[elg] |= desi_mask.ELG
    desi_target[lrg] |= desi_mask.LRG
    desi_target[sky] |= desi_mask.SKY
    desi_target[std] |= desi_mask.STD_FSTAR
    bgs_target = np.zeros(n, dtype='i8')    #- TODO
    mws_target = np.zeros(n, dtype='i8')    #- TODO

    #- True type
    truetype = np.zeros(n, dtype='S8')
    assert np.all(truetype == '')
    truetype[qso_lya | qso_tracer] = 'QSO'
    truetype[qso_fake] = 'STAR'
    truetype[elg] = 'ELG'
    truetype[lrg_real] = 'LRG'
    truetype[lrg_fake] = 'UNKNOWN'
    truetype[std] = 'STDSTAR'
    truetype[sky] = 'SKY'
    assert np.all(truetype != '')

    #- Misc other
    targetid = np.random.randint(2**62, size=n)
    lastpass = elg.astype('i4')
    ### brickname = np.zeros(n, dtype='S8')
    brickname = desispec.brick.brickname(ra, dec)

    mtl = Table()
    mtl['TARGETID'] = targetid
    mtl['BRICKNAME'] = brickname
    mtl['RA'] = ra
    mtl['DEC'] = dec
    mtl['NUMOBS'] = numobs_pre
    mtl['PRIORITY'] = prio_pre
    mtl['LASTPASS'] = lastpass
    mtl['DESI_TARGET'] = desi_target
    mtl['BGS_TARGET'] = bgs_target
    mtl['MWS_TARGET'] = mws_target

    truth = Table()
    truth['TARGETID'] = targetid
    truth['BRICKNAME'] = brickname
    truth['TRUEZ'] = z
    truth['TRUETYPE'] = truetype
    truth['PRIO_PRE'] = prio_pre
    truth['PRIO_POST'] = prio_post
    truth['NUMOBS_PRE'] = numobs_pre
    truth['NUMOBS_POST'] = numobs_post

    return mtl, truth

def rdzipn2mtl(infile='objects_ss_sf0.rdzipn', basename='mtl', clobber=False):
    """
    Converts input rdzipn file into MTL files:
    targets_mtl.fits, truth_mtl.fits, stdstars_mtl.fits, sky_mtl.fits
    and *_mtl_lite.fits versions for testing
    """
    truthfile = 'truth_{}.fits'.format(basename)
    targetfile = 'targets_{}.fits'.format(basename)
    stdstarfile = 'stdstars_{}.fits'.format(basename)
    skyfile = 'sky_{}.fits'.format(basename)
    truthfile_lite = 'truth_{}_lite.fits'.format(basename)
    targetfile_lite = 'targets_{}_lite.fits'.format(basename)
    stdstarfile_lite = 'stdstars_{}_lite.fits'.format(basename)
    skyfile_lite = 'sky_{}_lite.fits'.format(basename)

    #- Check if output files already exist
    ioerror = False
    for filename in (
        truthfile, targetfile, stdstarfile, skyfile,
        truthfile_lite, targetfile_lite, stdstarfile_lite, skyfile_lite ):
        if os.path.exists(filename):
            if clobber:
                os.remove(filename)
            else:
                print('{} already exists; use clobber=True to overwrite'.format(filename))
                ioerror = True
    if ioerror:
        raise IOError('Output files already exist; use clobber=True to overwrite')
                
    #- Read input rdzipn file
    mtl, truth = load_rdzipn(infile)
    iistd = truth['TRUETYPE'] == 'STDSTAR'
    iisky = truth['TRUETYPE'] == 'SKY'
    iitgt = ~(iistd | iisky)
    if not np.any(iistd):
        print("WARNING: no standard stars found")
    if not np.any(iisky):
        print("WARNING: no sky locations found")
    if not np.any(iitgt):
        print("WARNING: no science targets found")

    #- Write output files
    truth[iitgt].write(truthfile)
    mtl[iitgt].write(targetfile)
    mtl[iistd].write(stdstarfile)
    mtl[iisky].write(skyfile)

    #- Write lite version for testing
    lite = (mtl['RA'] <= 10) & (mtl['DEC'] <= 10) & (mtl['DEC'] >= -10)
    truth[iitgt & lite].write(truthfile_lite)
    mtl[iitgt & lite].write(targetfile_lite)
    mtl[iistd & lite].write(stdstarfile_lite)
    mtl[iisky & lite].write(skyfile_lite)

        
