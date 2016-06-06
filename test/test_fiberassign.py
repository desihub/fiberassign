#!/usr/bin/env python

"""
Test fiberassign code at NERSC
"""

from __future__ import print_function, division
import sys, os
from glob import glob

import optparse

parser = optparse.OptionParser(usage = "%prog [options]")
parser.add_option("--skipcompile", action="store_true", help="Don't test compilation")
opts, args = parser.parse_args()

#-------------------------------------------------------------------------
#- Basic setup of directories and environment
if 'FIBERASSIGN' in os.environ:
    fiberassign_dir = os.environ['FIBERASSIGN']
else:
    script = os.path.abspath(sys.argv[0])
    fiberassign_dir = os.path.dirname(os.path.dirname(script))

testdir = os.path.join(fiberassign_dir, 'test')

#- output directory
if 'FIBERASSIGN_OUTDIR' in os.environ:
    outdir = os.environ['FIBERASSIGN_OUTDIR']
elif 'SCRATCH' in os.environ:
    outdir = os.path.join(os.environ['SCRATCH'], 'desi', 'test_fiberassign')
else:
    outdir = os.path.join(fiberassign_dir, 'test', 'output')

print('Test output dir '+outdir)    

if not os.path.exists(outdir):
    os.makedirs(outdir)

#- dangerous
def remove_output(outdir):
    for fitsfile in glob(outdir+'/tile*.fits'):
        os.remove(fitsfile)
    
#-------------------------------------------------------------------------
#- test compilation
if not opts.skipcompile:
    print('---------------------------------------------------------')
    print('-- Testing compiling')
    os.chdir(fiberassign_dir)
    for command in ('make clean', 'make', 'make install'):
        err = os.system(command)
        assert err==0, "FAILED: " + command

#-------------------------------------------------------------------------
#- write parameter file ("features file") with scratch output directory
params = ''.join(open(testdir+'/template_fiberassign.txt').readlines())
paramfile = os.path.join(testdir, 'params_fiberassign.txt')
fx = open(paramfile, 'w')
fx.write(params.format(outdir=outdir))
fx.close()

#-------------------------------------------------------------------------
#- test fiberassign    
print('---------------------------------------------------------')
print('-- Testing fiberassign')
remove_output(outdir)
### command = 'export OMP_NUM_THREADS=24; srun -n 1 fiberassign '+paramfile
command = './bin/fiberassign '+paramfile
print(command)
err = os.system(command)
assert err==0, 'FAILED: '+command

#-------------------------------------------------------------------------
#- Test fiberassign_surveysim
print('---------------------------------------------------------')
print('-- Testing fiberassign_surveysim')
remove_output(outdir)
### command = 'export OMP_NUM_THREADS=24; srun -n 1 fiberassign_surveysim '+paramfile
command = './bin/fiberassign_surveysim '+paramfile
print(command)
err_surveysim = os.system(command)
assert err==0, 'FAILED: '+command

print('---------------------------------------------------------')
print('-- SUCCESS')
print('Output files in '+outdir)
