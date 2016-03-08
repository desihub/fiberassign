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

def rdzipn2targets(infile):
    """Read rdzipn infile and return target and truth tables
    """
    ra, dec, z, itype, priority, numobs = io.read_rdzipn(infile)
    n = len(ra)

    #- Martin's itype is 1 to n, while Bob's fiberassign is 0 to n-1
    itype -= 1
    assert np.min(itype >= 0)

    #- rdzipn has float32 ra, dec, but it should be float64
    ra = ra.astype('float64') % 360     #- enforce 0 <= ra < 360
    dec = dec.astype('float64')

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
    truetype = np.zeros(n, dtype='S10')
    assert np.all(truetype == '')
    truetype[qso_lya | qso_tracer] = 'QSO'
    truetype[qso_fake] = 'STAR'
    truetype[elg] = 'GALAXY'
    truetype[lrg_real] = 'GALAXY'
    truetype[lrg_fake] = 'UNKNOWN'
    truetype[std] = 'STAR'
    truetype[sky] = 'SKY'
    assert np.all(truetype != '')

    #- Misc other
    targetid = np.random.randint(2**62, size=n)
    ### brickname = np.zeros(n, dtype='S8')
    brickname = desispec.brick.brickname(ra, dec)
    subpriority = np.random.uniform(0, 1, size=n)

    targets = Table()
    targets['TARGETID'] = targetid
    targets['BRICKNAME'] = brickname
    targets['RA'] = ra
    targets['DEC'] = dec
    targets['DESI_TARGET'] = desi_target
    targets['BGS_TARGET'] = bgs_target
    targets['MWS_TARGET'] = mws_target
    targets['SUBPRIORITY'] = subpriority

    truth = Table()
    truth['TARGETID'] = targetid
    truth['BRICKNAME'] = brickname
    truth['RA'] = ra
    truth['DEC'] = dec
    truth['TRUEZ'] = z
    truth['TRUETYPE'] = truetype
    truth['CATEGORY'] = itype

    return targets, truth
