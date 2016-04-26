# this produces plots of the tiles periodically in an equal area Sanson Flamsteed projection
# steps gives the number of plots produced, with equal intervals
# the order is initially just the ordering of the tiles (plates) in desi-tiles-full.par
# but can be changed with the function reorder
# the output files are numbered tiling00.png, tiling01.png, ...tiling10.png,...tiling11.png...
# makes is easy to make a movie from them using ImageMagick, for example with

# convert -delay 100 -loop 1 tiling*.png tiling.gif
# the code expects to find the tiling file desi-tiles.par in the same directory
# to run, need to give file listing order of tiles, e.g.
# python universal_tiling_sequence_images.py default_survey_list.txt
# the last is  a sequence of tiles between 1 and 28810 from the full list desi-tiles.par
# with typically 10666 lines





#!/usr/bin/env python
import sys
import os
import numpy as N
from pylab import *
import matplotlib.pyplot as P



def readmulticolumn(f,names,types):
    for line in f:
        line=line.strip()
        if line.startswith('#') or len(line)==0:
            continue
        items=line.split()
        for j in range(len(names)):
            names[j].append(types[j](items[j]))
    return

def reorder(j):
    #stub for now.  allows us to use an order of the tiles other than the one in
    #desi-tiles.par
    return tile_order[j]

tileid=list()
ra=list()
dec=list()
passs=list()
in_desi=list()
ebv_med=list()
airmass=list()
exposefac=list()
no=list()

newfig=list()
tilenames=[tileid,no,ra,dec,passs,in_desi]
tiletypes=[str,int,float,float,int,int]
fstandard=open('/project/projectdirs/desi/software/cori/desimodel/0.4desi-tiles.par','r')
readmulticolumn(fstandard,tilenames,tiletypes)

# get ordering of tile from external file specified in execution command
forder=open(str(sys.argv[1]),'r')
tile_order=list()
for line in forder:
    tile_order.append([int(x) for x in line.split()][0])

for i in range(len(tile_order)):
    if(i<1000):
        print i,tile_order[i]

const=N.pi/180.
print len(tile_order)
len_tile_order=len(tile_order)
steps=30
y=list()
x=list()
color=list()
color_choices=['black','blue','red','green','yellow','magenta']
for k in range(steps):
 
    for j in range(k*len_tile_order/steps,(k+1)*len_tile_order/steps):
        i=reorder(j)-1
        if 1==1:
            this_ra=ra[i]
            if this_ra>300:
                this_ra-=360.
            y.append(dec[i])
            if(dec[i]<-20):
                print " error", k,j,no[i],ra[i],dec[i],in_desi[i]
                exit()
            x.append((-this_ra+120.)*N.cos(const*dec[i]))

            color.append(color_choices[passs[i]])
    fig=P.scatter(x,y,c=color,s=5,edgecolor='none')
    angle=[-60,30,120,210,300]
    for i in range(5):
        this_angle=angle[i]
        yy=linspace(-20,90,100)
        xx=(-this_angle+120.)*N.cos(const*yy)
 
        newfig.append(plot(xx,yy,color='0.03'))

    text(-20,-30,'RA=120')
    text(-180,-30,'RA=300')
    text(-100,-30,'RA=210')
    text(60,-30,'RA=30')
    text(140,-30,'RA=-60')
    fig.axes.get_xaxis().set_ticks([])
    ylabel('DEC (deg)')
    if(k<10):
        savefig('tiling00'+str(k)+'.png')
    elif(k<100):
        savefig('tiling0'+str(k)+'.png')
    else:
        savefig('tiling'+str(k)+'.png')





