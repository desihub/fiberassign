#!/usr/bin/env python

"""
I/O routines for fiber assignment
"""

import sys, os
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column

def write_rdzipn(filename, ra, dec, z, itype, priority, numobs):
    """
    Writes the custom "rdzipn" file format
    
    Args (all of which are numpy arrays):
        ra: Right Ascension [degrees]
        dec : declination [degrees]
        z : redshift (float)
        itype : integer true object type [0-7]
        priority : integer observation priority
        numobs : requested number of observations (integer)
        
    TODO: the current rdzipn file format uses float32 for RA, dec, which isn't
    high enough precision.  Implementing this for backwards compatiblity right
    now, but this needs to be fixed.
    """
    fout = open(filename, "w")
    Nt      = np.array([ra.size],dtype='i4')
    Nt.tofile(fout)
    ra.astype('f4').tofile(fout)
    dec.astype('f4').tofile(fout)
    z.astype('f4').tofile(fout)
    itype.astype('i4').tofile(fout)
    priority.astype('f4').tofile(fout)
    numobs.astype('i4').tofile(fout)
    fout.close()
    
def read_rdzipn(filename):
    """Read filename in custom rdzipn format
    
    Returns arrays of ra, dec, z, itype, priority, numobs
    """
    fx = open(filename, "r")
    Nt = np.fromfile(fx, dtype='i4', count=1)
    ra = np.fromfile(fx, dtype='f4', count=Nt)
    dec = np.fromfile(fx, dtype='f4', count=Nt)
    z = np.fromfile(fx, dtype='f4', count=Nt)
    itype = np.fromfile(fx, dtype='i4', count=Nt)
    priority = np.fromfile(fx, dtype='f4', count=Nt)
    numobs = np.fromfile(fx, dtype='i4', count=Nt)
    fx.close()
    return ra, dec, z, itype, priority, numobs

def rdzipn2mtl(rdzipn, mtl):
    """
    Convert a ...
    """
    raise NotImplementedError
# mtl = Table()
# n = len(ra)
# mtl.add_column(Column(np.arange(n), name='TARGETID'))
# mtl.add_column(Column(np.zeros(n, dtype='S8'), name='BRICKNAME'))  #- TODO
# mtl.add_column(Column(ra, name='RA'))
# mtl.add_column(Column(ra, name='DEC'))
# mtl.add_column(Column(numobs, name='NUMOBS'))
# mtl.add_column(Column(priority, name='PRIORITY'))
# mtl.add_column(Column(np.zeros(n, dtype='i4'), name='LASTPASS'))     #- TODO
# mtl.add_column(Column(np.zeros(n, dtype='i8'), name='DESI_TARGET'))  #- TODO
# mtl.add_column(Column(np.zeros(n, dtype='i8'), name='MWS_TARGET'))   #- TODO
# mtl.add_column(Column(np.zeros(n, dtype='i8'), name='BGS_TARGET'))   #- TODO

def write_mtl(filename, table):
    raise NotImplementedError

def read_mtl(filename):
    return table.read(filename, 1)

def write_truth(filename):
    raise NotImplementedError

def read_truth(filename):
    raise NotImplementedError
    