



#!/usr/bin/env python
import sys
import os
import numpy as N
from pylab import *
import matplotlib.pyplot as P
import matplotlib.patches as mpatches


def readmulticolumn(f,names,types):
    for line in f:
        line=line.strip()
        if line.startswith('#') or len(line)==0:
            continue
        items=line.split()
        for j in range(len(names)):
            names[j].append(types[j](items[j]))
    return

ra=list()
dec=list()
av_tfs1=list()
av_tfs2=list()
newfig=list()
names=[ra,dec,av_tfs1,av_tfs2]
types=[float,float,int,int]


tiling=input("Standard tiling=0   Randomized tiling=1")
if(tiling==0):
    fstandard=open('/project/projectdirs/desi/mocks/preliminary/new_random_map.txt','r')
else:
    fstandard=open('/project/projectdirs/desi/mocks/preliminary/random_tiles_new_random_map.txt','r')

readmulticolumn(fstandard,names,types)


mra=N.array(ra)
mdec=N.array(dec)
mav_tfs1=N.array(av_tfs1)
mav_tfs2=N.array(av_tfs2)

while (1==1):

    ra_min=input("ra_min = ")
    ra_max=input("ra_max = ")
    dec_min=input("dec_min = ")
    dec_max=input("dec_max = ")
    spotsize=input("spotsize = ")
    ra_center=(ra_min+ra_max)/2.
    dec_center=(dec_min+dec_max)/2.

    overlap=input("Don't include overlap of positioners=0  Do include overlap of positioners")
    if (ra_min>0 & ra_max>0):
        ii=(mra<ra_max)&(mra>ra_min)&(mdec<dec_max)&(mdec>dec_min)
    if (ra_min<0 & ra_max>0):
        ii=(mdec<dec_max)&(mdec>dec_min)&(mra<ra_max|mra>ra_mix+360)
    if (ra_min<0 & ra_max<0):
        ii=(360+mra<ra_max)&(360+mra>ra_min)&(mdec<dec_max)&(mdec>dec_min)
    nra=mra[ii]
    ndec=mdec[ii]
    if(overlap==0):
        nav_tfs=mav_tfs2[ii]
    else:
        nav_tfs=mav_tfs2[ii]

    max_av_tfs=N.amax(nav_tfs)
    min_av_tfs=N.amin(nav_tfs)
    values=N.arange(min_av_tfs,max_av_tfs+1)
    for i in range(100): print nra[i],ndec[i],nav_tfs[i]
    print len(nra)
    length=len(nra)
    #convert to radians
    convert=3.1415296/180
    nra=nra*convert
    ndec=ndec*convert
    ra_center=ra_center*convert
    dec_center=dec_center*convert
    ndec_center=N.empty(length); ndec_center.fill(dec_center)
    nra_center=N.empty(length); nra_center.fill(ra_center)
    delta_ra=nra-nra_center
    delta_dec=ndec-ndec_center
    sindec=N.sin(ndec)
    cosdec=N.cos(ndec)
    sindeltar=N.sin(delta_ra)
    cosdeltar=N.cos(delta_ra)
    cosdzero=N.cos(dec_center)
    sindzero=N.sin(dec_center)
    
    a=cosdec*sindeltar
    b=sindec*cosdzero-sindzero*cosdec*cosdeltar
    c=sindzero*sindec+cosdzero*cosdec*cosdeltar
    
    x=a/c
    y=b/c
    x=x/convert+ra_center/convert
    y=y/convert+dec_center/convert
    f,(ax1,ax2)=plt.subplots(1,2,figsize=(14,7))

    P.subplot(1,2,1)
    
    def cmap_discretize(cmap, N):
        if type(cmap) == str:
            cmap = get_cmap(cmap)
        colors_i = concatenate((linspace(0, 1., N), (0.,0.,0.,0.)))
        colors_rgba = cmap(colors_i)
        indices = linspace(0, 1., N+1)
        cdict = {}
        for ki,key in enumerate(('red','green','blue')):
            cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]
        # Return colormap object.
        return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

    cmap=cmap_discretize('jet',max_av_tfs-min_av_tfs+1)

    P.scatter(x,y,c=nav_tfs-min_av_tfs,s=spotsize,edgecolor='none',cmap=cmap)
    P.axis("scaled")
    P.xlim(ra_min,ra_max)
    P.ylim(dec_min,dec_max)
    P.gca().invert_xaxis()
    P.xlabel("RA",fontsize=16)
    P.ylabel("dec",fontsize=16)
    mappable = cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(min_av_tfs-0.5,max_av_tfs+0.5)
    bar=P.colorbar(mappable)
    bar.set_ticks(N.linspace(min_av_tfs,max_av_tfs,len(values)))
    bar.set_ticklabels(values)
    

    P.subplot(1,2,2)
    counts=N.zeros(15)
    for i in range(15):
        counts[i]=sum(nav_tfs==i)
    norm_counts=counts/length
    print norm_counts
    nsq=nav_tfs*nav_tfs

    summ=float(sum(nav_tfs))
    mean=summ/float(length)
    rms=(float(sum(nsq)-summ*summ/length))/length
    print length, mean, rms

    index=N.arange(15)
    P.bar(index,norm_counts,align='center')
    P.xlim([0,14])
    P.text(8,0.2,"mean = {:10.3f}".format(mean))
    P.text(8,0.17,"rms   = {:10.3f}".format(rms))
    P.xlabel("Fraction of Area Covered n Times",fontsize=16)
    show()

