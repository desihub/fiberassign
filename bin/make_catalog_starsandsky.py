#!/usr/bin/env python
# revised 4/10/15
# writes output in FITS format 
#
# revised 7/29/14
# Generates a catalog by combining existing mock catalogs
# for each type into a binary format accepted by "assign".
#    based on code by Martin White
# modified to achieve densities used in CDR and to include
# standard stars and sky fibers
# goal densities:
# QSOI (z<2.15) = 120/sq deg
# QSOII (z>2,15) = 50/sq deg
# bad QSO   =90 /sq deg
# LRG  300/sq deg good, 50/sq deg bad
# ELG  2400/dq deg
# 03/05/15  need to add standard stars and sky fibers
# can take as random, with correct density
# Schlegel suggest 20 SS per petal so 20 x 5.3(multiple coverage)/0.75 sq deg=140 /sq deg
# for sky fibers Scheleg suggests 200, i.e. 1400/sq deg
#
# densities in Martin's mocks v2 on 06/25/14
# ELGs 2607/sq deg  > 2611
# LRGs 301/sq deg   >  323
# QSOI  110/sq deg  > 118
# QSO II 63.8/sq deg  >
# random 632 /sq deg

import numpy  as N
import fitsio as F
import sys
from astropy.io import fits
import os.path
from desitarget import desi_mask, bgs_mask, mws_mask
import numpy.lib.recfunctions as rfn

__author__ = "Martin White modified by Robert Cahn  6/30/14"
__version__ = "1.0"
__email__  = "mjwhite@lbl.gov/rncahn@lbl.gov"

dbase = "/project/projectdirs/desi/mocks/preliminary/"

footprint_area=20.*45.*N.sin(45.*N.pi/180.)/(45.*N.pi/180.)
print("footprint area %f"%(footprint_area))

total_area=0.
def read_objects(fn):
    """
    read_objects(fn):
    Reads the RA, DEC and Z of the objects from a FITS file.
    """
    data = F.read(fn,columns=['RA','DEC','Z'])
    ra   = data[ 'RA'].astype('f4')
    dec  = data['DEC'].astype('f4')
    zz   = data[  'Z'].astype('f4')
    """
        use small region to determine densities of mocks
    """
    smalldata=data[(data['RA']>170.) & (data['RA']<190.) & (data['DEC']>0.) & (data['DEC']<45.)]
    galaxies=len(smalldata)

    density=galaxies/footprint_area
    return( (ra,dec,zz,density) )
    #

def read_objects_noZ(fn):
    """
        read_objects(fn):
        Reads the RA, DEC of the objects from a FITS file.
        """
    data = F.read(fn,columns=['RA','DEC'])
    ra   = data[ 'RA'].astype('f4')
    dec  = data['DEC'].astype('f4')
    zz    = N.zeros(len(data.ra),dtype='f4')
    smalldata=data[(data['RA']>170.) & (data['RA']<190.) & (data['DEC']>0.) & (data['DEC']<45.)]
    galaxies=len(smalldata)

    density=galaxies/footprint_area
    return( (ra,dec,zz),density )
#

def reduce(ra,dec,zz,frac):
    xra=N.array(ra)
    xdec=N.array(dec)
    xzz=N.array(zz)
   
    keepornot=N.random.uniform(0.,1.,len(ra))
    limit=N.zeros(len(xra)) +frac
    #create boolean array of which to keep
    #find new length
    kept=keepornot<limit
    yra=xra[kept]
    ydec=xdec[kept]
    yzz=xzz[kept]
    
    return( (yra,ydec,yzz))




def write_mtl(ra, dc, no, pp, lastpass, zz, types, icat=0, basename="targets"):
    fitsname="{}_{}_mtl_0.fits".format(basename,icat)
    print("writing to file {}".format(fitsname))
    if(os.path.exists(fitsname)):
        os.remove(fitsname)
        # structure for output data
    type_table = [
        ('TARGETID', '>i4'), 
        ('BRICKNAME', '|S8'),
        ('RA', '>f4'), 
        ('DEC', '>f4'),
        ('NUMOBS', '>i4'), 
        ('PRIORITY', '>i4'),
        ('LASTPASS', '>i4')
        ]
    

    data = N.ndarray(shape=ra.size, dtype=type_table) 
    data['TARGETID'] = N.arange(ra.size)
    data['RA'] = ra
    data['DEC'] = dc
    data['BRICKNAME'][:] = "00000000"
    data['NUMOBS'] = no
    data['PRIORITY'] = pp
    data['LASTPASS'] = lastpass
    
    
    desi_target = N.int_(types)
    bgs_target = N.zeros(ra.size, dtype=N.int64)
    mws_target = N.zeros(ra.size, dtype=N.int64)    


    data = rfn.append_fields(data,
                             ['DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET'],
                             [desi_target, bgs_target, mws_target], usemask=False)
    
        #- Create header to include versions, etc.
    hdr = F.FITSHDR()
    hdr['DEPNAM00'] = 'makecatalog'
    F.write(fitsname, data, extname='MTL', header=hdr, clobber=True)
    print('wrote {} items to target file'.format(len(ra)))
    
        #structure for truth file
    type_table = [
        ('TARGETID', '>i8'),
        ('BRICKNAME', '|S20'),
        ('Z', '>f8'),
        ('TYPE', '>i8')
        ]
    data = N.ndarray(shape=ra.size, dtype=type_table) 
    data['TARGETID'] = N.arange(ra.size)
    data['BRICKNAME'][:] = "00000000"
    data['Z'] = zz
    data['TYPE'] = types

    hdr = F.FITSHDR()
    hdr['DEPNAM00'] = 'makecatalog-truth'
    F.write('truth_'+fitsname, data, extname='MTL', header=hdr, clobber=True)
    print('wrote {} items to truth target file'.format(len(ra)))
    return 


def write_catalog(icat=0):
    """
    write_catalog(icat=0):
    Writes the catalog information to a file.
    This fills in the number of observations required and the
    priorities as well as the position information from the mock
    catalog files.
    """
    


    goal_qsoI=120.
    goal_qsoII=50.
    goal_badqso=90.
    goal_lrg=300.
    goal_badqso=90.
    goal_badlrg=50.
    goal_elg=2400.
    goal_standardstar=140.
    goal_skyfiber=1400.
    
    possible_types = ['ELG', 'LRG', 'QSO', 'STDSTAR', 'GAL', 'OTHER']
    types = N.empty((0))

    # LyA QSOs are priority 1, other QSOs are priority 2.
    xra,xdc,xzz,density = read_objects(dbase+"v2_qso_%d.fits"%icat)

    low=xzz<2.1
    high=xzz>2.1
    ra_low=xra[low]
    dc_low=xdc[low]
    zz_low=xzz[low]
    total_area=len(xra)/density
    print 'total area in sq deg',total_area
    ra_high=xra[high]
    dc_high=xdc[high]
    zz_high=xzz[high]
    mockdensity_qsoI=density*len(ra_low)/(len(ra_high)+len(ra_low))
    mockdensity_qsoII=density*len(ra_high)/(len(ra_high)+len(ra_low))

    print ('mock density qsoI  %8.1f  mock density qsoII %8.1f\n'%(mockdensity_qsoI,mockdensity_qsoII))
    # need to increase number of qsois probably'
    frac_qsoI=goal_qsoI/mockdensity_qsoI
    qsoI_needed=frac_qsoI-1.
    print 'need %8.3f more qsoI \n' %qsoI_needed
    print 'first qsoIs',len(ra_low)
    
    id       = N.zeros(ra_low.size,dtype='i4') + 2
    pp       = N.zeros(ra_low.size,dtype='f4') + desi_mask['QSO'].priorities['UNOBS'] 
    no       = N.zeros(ra_low.size,dtype='i4') + 1

    frac_qsoII=goal_qsoII/mockdensity_qsoII
    print ('goals  qsoI   %8.1f   qsoII  %8.1f\n'%(goal_qsoI,goal_qsoII))
    nra,ndc,nzz = reduce(ra_high,dc_high,zz_high,frac_qsoII)
    #combine all qsoIs with qsoIIs
    nid       = N.zeros(nra.size,dtype='i4') + 1
    npp       = N.zeros(nra.size,dtype='f4') + desi_mask['QSO'].priorities['UNOBS'] 
    nno       = N.zeros(nra.size,dtype='i4') + 5
    ra       = N.append(ra_low,nra)
    dc       = N.append(dc_low,ndc)
    zz       = N.append(zz_low,nzz)
    id       = N.append(id,nid)
    pp       = N.append(pp,npp)
    no       = N.append(no,nno)
    print 'qsoIIs',len(nra)

    tmp_type = N.ones(ra.size, dtype='i8') * desi_mask.QSO
    types = N.append(types, tmp_type)
    print  'total ra, types', N.size(ra), N.size(types) 
    
    # get extra qsos for qsoI
    icatp=icat+1
    xra,xdc,xzz,density = read_objects(dbase+"v2_qso_%d.fits"%icatp)
    
    low=xzz<2.1
    nra_low=xra[low]
    ndc_low=xdc[low]
    nzz_low=xzz[low]
    nra,ndc,nzz=reduce(nra_low,ndc_low,nzz_low,qsoI_needed)
    nid       = N.zeros(nra.size,dtype='i4')+2
    npp       = N.zeros(nra.size,dtype='f4') + desi_mask['QSO'].priorities['UNOBS'] 
    nno       = N.zeros(nra.size,dtype='i4')+1
    ra       = N.append(ra,nra)
    dc       = N.append(dc,ndc)
    zz       = N.append(zz,nzz)
    id       = N.append(id,nid)
    pp       = N.append(pp,npp)
    no       = N.append(no,nno)
    

    tmp_type = N.ones(nra.size, dtype='i8') * desi_mask.QSO
    types = N.append(types, tmp_type)
    print  'total ra, types', N.size(ra), N.size(types) 
    print  'total ra, types', len(ra), N.size(types) 
    print' added qsoIs', len(nra)

    # LRGs are priority 3.
    # density is just right
    nra,ndc,nzz,density=read_objects(dbase+"v2_lrg_%d.fits"%icat)
    print 'lrg mock density', density

    nid      = N.zeros(nra.size,dtype='i4') + 3
    npp      = N.zeros(nra.size,dtype='f4') + desi_mask['LRG'].priorities['UNOBS'] 
    nno      = N.zeros(nra.size,dtype='i4') + 2   #only use 2 exposures for LRG
    ra       = N.append(ra,nra)
    dc       = N.append(dc,ndc)
    zz       = N.append(zz,nzz)
    id       = N.append(id,nid)
    pp       = N.append(pp,npp)
    no       = N.append(no,nno)


    tmp_type = N.ones(nra.size, dtype='i8') * desi_mask.LRG
    types = N.append(types, tmp_type)

    print 'lrgs added',len(nra)
    print 'lrg density',len(nra)/total_area
    print  'total ra, types', len(ra), N.size(types) 

    # ELGs are priority 4.
    
    mra,mdc,mzz,density=read_objects(dbase+"v2_elg_%d.fits"%icat)
    print 'mock density elgs', density
    
    nra,ndc,nzz=reduce(mra,mdc,mzz,(goal_elg/density))
    
    
    nid      = N.zeros(nra.size,dtype='i4') + 4
    npp      = N.zeros(nra.size,dtype='f4') + desi_mask['ELG'].priorities['UNOBS'] 
    nno      = N.zeros(nra.size,dtype='i4') + 1
    ra       = N.append(ra,nra)
    dc       = N.append(dc,ndc)
    zz       = N.append(zz,nzz)
    id       = N.append(id,nid)
    pp       = N.append(pp,npp)
    no       = N.append(no,nno)

    tmp_type = N.ones(nra.size, dtype='i8') * desi_mask.ELG
    types = N.append(types, tmp_type)    

    print 'elgs added',len(nra)
    print 'new total nid'
    print 'elg density',len(nra)/total_area
    print  'total ra, types', len(ra), len(types)     

    
    # and now we have "fake" qsos, placed randomly.
    data     = F.read(dbase+"v2_randoms_big.fits",columns=['RA','DEC'])
    smalldata=data[(data['RA']>170.) & (data['RA']<190.) & (data['DEC']>0.) & (data['DEC']<45.)]
    galaxies=len(smalldata)

    mockdensity_randoms=galaxies/footprint_area
    print 'randoms density', mockdensity_randoms
    end_badqso=int(goal_badqso*total_area)
    end_badlrg=int(end_badqso+goal_badlrg*total_area)
    end_standardstar=int(end_badlrg+goal_standardstar*total_area)
    end_skyfiber=int(end_standardstar+goal_skyfiber*total_area)
    print 'end_badqso %d end_badlrg %d\n'%(end_badqso,end_badlrg)
    nra      = data[ 'RA'][:end_badqso].astype('f4')
    ndc      = data['DEC'][:end_badqso].astype('f4')
    nzz      = N.zeros(nra.size,dtype='f4')
    nid      = N.zeros(nra.size,dtype='i4')+5
    npp      = N.zeros(nra.size,dtype='f4') + desi_mask['QSO'].priorities['UNOBS'] 
    nno      = N.zeros(nra.size,dtype='i4')+1
    ra       = N.append(ra,nra)
    dc       = N.append(dc,ndc)
    zz       = N.append(zz,nzz)
    id       = N.append(id,nid)
    pp       = N.append(pp,npp)
    no       = N.append(no,nno)
    density=len(nra)/total_area

    tmp_type = N.ones(nra.size, dtype='i8') * desi_mask.QSO
    types = N.append(types, tmp_type)    
    print 'fake qso density',density

    print 'fake qsos',len(nra), len(types)
    print  'total ra, types', len(ra), len(types) 

    #now need bad lrgs at 50/sq deg
    nra      = data[ 'RA'][end_badqso+1:end_badlrg].astype('f4')
    ndc      = data['DEC'][end_badqso+1:end_badlrg].astype('f4')
    nzz      = N.zeros(nra.size,dtype='f4')
    nid      = N.zeros(nra.size,dtype='i4') + 6
    npp      = N.zeros(nra.size,dtype='f4') + desi_mask['LRG'].priorities['UNOBS'] 
    nno      = N.zeros(nra.size,dtype='i4') + 1
    ra       = N.append(ra,nra)
    dc       = N.append(dc,ndc)
    zz       = N.append(zz,nzz)
    id       = N.append(id,nid)
    pp       = N.append(pp,npp)
    no       = N.append(no,nno)
    density=len(nra)/total_area

    tmp_type = N.ones(nra.size, dtype='i8') * desi_mask.LRG
    types = N.append(types, tmp_type)
    
    print 'fake qso density',density

    print 'fake lrg density',density
    
    print 'fake lrgs added', len(nra), len(types)
    print  'total ra, types', len(ra), len(types) 
    lastpass = N.zeros(ra.size, dtype='int')    

    #now need standard stars at 140/sq deg
    star_ra      = data[ 'RA'][end_badlrg+1:end_standardstar].astype('f4')
    star_dc      = data['DEC'][end_badlrg+1:end_standardstar].astype('f4')
    star_zz      = N.zeros(star_ra.size,dtype='f4')
    star_id      = N.zeros(star_ra.size,dtype='i4')+7
    star_pp      = N.zeros(star_ra.size,dtype='f4') + 100
    star_no      = N.zeros(star_ra.size,dtype='i4') + 10
    star_type = N.ones(star_ra.size, dtype='i8') * desi_mask.STD_FSTAR
    star_lastpass = N.zeros(star_ra.size, dtype='int')    

    density=len(star_ra)/total_area
    print 'standardstar density',density
    
    #now need sky fibers at 1400/sq deg
    sky_ra      = data[ 'RA'][end_standardstar+1:end_skyfiber].astype('f4')
    sky_dc      = data['DEC'][end_standardstar+1:end_skyfiber].astype('f4')
    sky_zz      = N.zeros(sky_ra.size, dtype='f4')
    sky_id     = N.zeros(sky_ra.size, dtype='i4')+8
    sky_pp      = N.zeros(sky_ra.size, dtype='f4') + 200
    sky_no      = N.zeros(sky_ra.size, dtype='i4')+10
    sky_type   = N.ones(sky_ra.size, dtype='i8') * desi_mask.SKY
    sky_lastpass = N.zeros(sky_ra.size, dtype='int')    

    density=len(sky_ra)/total_area
    print 'sky fiber density',density    

   
    #
    print "Writing information for ",ra.size," objects."
    fout = open("objects_ss_sf%d.rdzipn"%icat,"w")
    Nt      = N.array([ra.size],dtype='i4')
    Nt.tofile(fout)
    ra.tofile(fout)
    dc.tofile(fout)
    zz.tofile(fout)
    id.tofile(fout)
    pp.tofile(fout)
    no.tofile(fout)
    fout.close()
    #
    print  'total ra, types', len(ra), len(types) 
    write_mtl(ra, dc, no, pp, lastpass, zz, types, basename="targets", icat=icat)
    write_mtl(star_ra, star_dc, star_no, star_pp, star_lastpass, star_zz, star_type, basename="stdstar", icat=icat)
    write_mtl(sky_ra, sky_dc, sky_no, sky_pp, sky_lastpass, sky_zz, sky_type, basename="sky", icat=icat)

if __name__=="__main__":
    write_catalog()
