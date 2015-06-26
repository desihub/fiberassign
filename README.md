Fiber assignment code for DESI
------------------------------------------------------
In a nutshell (to run it) :

cd src

make all

cd ..

./assign features.txt (or on NERSC : "qsub run", script that does it on your NERSC account)

------------------------------------------------------
 On Nersc, load modules PrgEnv-gnu and cfitsio

------------------------------------------------------

Modify parameters in features.txt (input files, etc...)


------------------------------------------------------

"talk0" corresponds to version Vtalk0 (branched on github on June 25th 2015) which works

To create Tikz graphs, process .dat files like in makepdf.tex

