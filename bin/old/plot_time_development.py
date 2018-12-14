#!/usr/bin/env python

# this takes the file time2.dat from fiberassign/doc/figs and plots the number of various
# galaxy types observed as a function of the number of tiles observed
# it is therefore sensitive to the survey strategy defined by surveyFile in
# the surveysim code

# to run:  python plot_time_development.py  file.ext
# the last argument is the name of the time2.dat file to be used

import sys
import os
import numpy as N
import pylab as P
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

def readmulticolumn(f,names,types):
    for line in f:
        line=line.strip()
        if line.startswith('#') or len(line)==0:
            continue
        items=line.split()
        for j in range(len(names)):
            names[j].append(types[j](items[j]))
    return

ntile=list()
qso1=list()
qso2=list()
qso3=list()
qso4=list()
qso5=list()
lrg1=list()
lrg2=list()
qsotr=list()
elg=list()
other=list()
tilenames=[ntile,qso1,qso2,qso3,qso4,qso5,lrg1,lrg2,qsotr,elg]
tiletypes=[int,int,int,int,int,int,int,int,int,int]
#f=open('time2.dat','r')
f=open(str(sys.argv[1]),'r')
readmulticolumn(f,tilenames,tiletypes)
plt.figure(1)
figa, ax = plt.subplots()
ax.set_yscale("log")
ax.set_ylim(1e4,1e8)
plt.xlabel("tiles observed")
plt.ylabel("galaxies observed")

plt.plot(ntile,qso1,'r-',linewidth=2.0)
plt.plot(ntile,qso2,'r-',linewidth=2.0)
plt.plot(ntile,qso3,'r-',linewidth=2.0)
plt.plot(ntile,qso4,'r-',linewidth=2.0)
plt.plot(ntile,qso5,'r-',linewidth=2.0)
plt.plot(ntile,lrg1,'b--',linewidth=2.0)
plt.plot(ntile,lrg2,'g--',linewidth=2.0)
plt.plot(ntile,qsotr,'c-',linewidth=2.0)
plt.plot(ntile,elg,'k-',linewidth=2.0)
plt.draw()
plt.savefig(sys.argv[1]+".linear.pdf")
plt.figure(2)
figb, bx = plt.subplots()
bx.set_ylim(0,5e6)
bx.set_xlim(0,4000)
bx.set_yscale("linear")
plt.xlabel("tiles observed")
plt.ylabel("galaxies observed")
plt.plot(ntile,qso1,'r-',linewidth=2.0)
plt.plot(ntile,qso2,'r-',linewidth=2.0)
plt.plot(ntile,qso3,'r-',linewidth=2.0)
plt.plot(ntile,qso4,'r-',linewidth=2.0)
plt.plot(ntile,qso5,'r-',linewidth=2.0)
plt.plot(ntile,lrg1,'b--',linewidth=2.0)
plt.plot(ntile,lrg2,'g--',linewidth=2.0)
plt.plot(ntile,qsotr,'c-',linewidth=2.0)
plt.plot(ntile,elg,'k-',linewidth=2.0)
plt.draw()
plt.savefig(sys.argv[1]+".log.pdf")


    
