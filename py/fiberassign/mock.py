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
from desitarget import desi_mask as M
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

    #- Also promote itype i4 -> i8
    itype = itype.astype('i8')

    #- Hack mapping of priorities
    newpriority = np.zeros(n, dtype='i4')
    qso = (itype==0) | (itype==1) | (itype==4)
    lrg = (itype==2) | (itype==5)
    elg = (itype==3)
    star = (itype==6)
    sky = (itype==7)

    #- TODO: replace these priorities with desitarget priorities
    newpriority[qso] = 2000
    newpriority[lrg] = 3000
    newpriority[elg] = 4000
    newpriority[star] = 9900
    newpriority[sky] = 9800
    assert np.all(newpriority > 0)

    #- Create a DESI_TARGET mask
    desi_target = np.zeros(n, dtype='i8')
    desi_target[qso] |= M.QSO
    desi_target[elg] |= M.ELG
    desi_target[lrg] |= M.LRG
    desi_target[sky] |= M.SKY
    desi_target[star] |= M.STD_FSTAR
    bgs_target = np.zeros(n, dtype='i8')    #- TODO
    mws_target = np.zeros(n, dtype='i8')    #- TODO

    targetid = np.random.randint(2**62, size=n)
    lastpass = elg.astype('i4')
    brickname = desispec.brick.brickname(ra, dec)

    mtl = Table()
    mtl.add_column(Column(np.arange(n), name='TARGETID'))
    mtl.add_column(Column(brickname,    name='BRICKNAME'))
    mtl.add_column(Column(ra,           name='RA'))
    mtl.add_column(Column(dec,          name='DEC'))
    mtl.add_column(Column(numobs,       name='NUMOBS'))
    mtl.add_column(Column(newpriority,  name='PRIORITY'))
    mtl.add_column(Column(lastpass,     name='LASTPASS'))
    mtl.add_column(Column(desi_target,  name='DESI_TARGET'))
    mtl.add_column(Column(bgs_target,   name='BGS_TARGET'))
    mtl.add_column(Column(mws_target,   name='MWS_TARGET'))

    truth = Table()
    truth.add_column(Column(np.arange(n), name='TARGETID'))
    truth.add_column(Column(brickname, name='BRICKNAME'))
    truth.add_column(Column(z, name='Z'))
    truth.add_column(Column(itype, name='TYPE'))

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
    iitgt = truth['TYPE'] < 6
    iistd = truth['TYPE'] == 6
    iisky = truth['TYPE'] == 7

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

        
