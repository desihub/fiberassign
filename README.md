Fiber assignment code for DESI
------------------------------------------------------
In a nutshell (to run it) :

cd src

make all

cd ..

./assign features.txt (or on NERSC : "qsub run", script that does it on your NERSC account)

------------------------------------------------------

Modify parameters in features.txt (input files, etc...)


------------------------------------------------------

"talk0" corresponds to version Vtalk0 (branched on github) which works

Compile doc.tex with luatex

To create Tikz graphs, process .dat files like in makepdf.tex

