#!/usr/bin/env python
import sys
import os
import numpy as N
import pylab as P
import matplotlib as M

f=open('desi-tiles.par','r')
fout=open('survey_tiling_order.txt','w')
fsurvey=open('default_survey_list.txt','w')
def readmulticolumn(f,names,types):
    for line in f:
        line=line.strip()
        if line.startswith('#') or len(line)==0:
            continue
        items=line.split()
        for j in range(len(names)):
            names[j].append(types[j](items[j]))
    return
tileid=list()
ra=list()
dec=list()
passs=list()
in_desi=list()
ebv_med=list()
airmass=list()
exposefac=list()
no=list()
tilenames=[tileid,no,ra,dec,passs,in_desi,ebv_med,airmass,exposefac]
tiletypes=[str,int,float,float,int,int,float,float,float]
#read canonical list of tiles
readmulticolumn(f,tilenames,tiletypes)
ten_by_ten=list()
no_tiles=len(ra)
ten=list()
print "  number of tiles = %d "%no_tiles
#cover full sky i for dec j for ra
ten_by_ten=[[[] for i in range(18)] for j in range(36)]
    
count_in_desi_x=0    
count_in_desi=0                   
for i in range(no_tiles):
    #assign tiles to 10 deg by 10 deg patches
    Nra=int(ra[i])
    if Nra<0:
        Nra+=360
    Ndec=int(dec[i])+90
    Nra_ten=Nra/10
    Ndec_ten=Ndec/10
    if Ndec_ten==18:
        Ndec_ten-=1
    if Ndec_ten<18:
        if Nra_ten<36:
            #only consider tiles in desi footprint
            if in_desi[i]==1:
                count_in_desi_x+=1
                ten_by_ten[Nra_ten][Ndec_ten].append(no[i]-1)

print " tiles used %d "%(count_in_desi_x)
count=0
total=0
for i in range(36):
    for j in range(18):
        #if there are tiles in this patch
        if len(ten_by_ten[i][j])>0:
            print "%d %d  %d " %(i,j,len(ten_by_ten[i][j]))
            count+=1
            total+=len(ten_by_ten[i][j])
print " count = %d total = %d "%(count,total)

## distribution of tiles by dec

for j in range(18):
    tiles_by_dec=0
    for i in range(36):
        tiles_by_dec+=len(ten_by_ten[i][j])
    print " j= %d %d  tiles in this dec interval dec= %d to dec %d  "%(j,tiles_by_dec,10*(j-9),10*(j-8))  




def add_list(passs,ipass,dec,ra):
    add=[]

    current=ten_by_ten[ra][dec]
    for m in current:
        #this selects gray-time targets by their tileid mod 4
        if (passs[m]==ipass or (passs[m]==5 and m  % 4 ==ipass-1)):
            add.append(m)
    return  add          

full_tile_list=[]
for ipass in range(1,5):
    dec_list=[7,8,9,10,11,12,13,14,15,16,17]
    for xdec in dec_list:
        for xra in range(36):
            to_add=add_list(passs,ipass,xdec,xra)
            full_tile_list.extend(to_add) 
        print "dec %d tile so far %d"%(xdec,len(full_tile_list))

print "  total tiles = %d"%len(full_tile_list)

## write new list of tiles

for i in range(len(full_tile_list)):
    k=no[full_tile_list[i]]-1
    fout.write( "STRUCT1 %5d %9.5f  %8.6f %2d %2d   %9.7f  %6.5f %6.5f\n" %( no[k],ra[k],dec[k],passs[k],in_desi[k],ebv_med[k],airmass[k],exposefac[k]))
    fsurvey.write("%5d\n"%no[k])
fout.close()
fsurvey.close()
  
    
    
