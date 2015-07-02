Fiber assignment code for DESI
------------------------------------------------------

If the master version does not work, you may want to use the Vtalk0 branch version, which works

------------------------------------------------------

On Nersc, load modules PrgEnv-gnu and cfitsio like that (to compile/run) :

load module PrgEnv-gnu

load module cfitsio

------------------------------------------------------

Compilation/execution : (after having loaded right modules)

cd src

make all

cd ..

./assign features.txt (or on NERSC : "qsub run", script that does it on your NERSC account)

------------------------------------------------------

Modify parameters in features.txt (input files, etc...)

------------------------------------------------------

"talk0" corresponds to version Vtalk0 (branched on github on June 25th 2015) which works

To create Tikz graphs, process .dat files like in makepdf.tex
