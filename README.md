# fiberassign
Fiber assignment code for DESI

------------------------------------------------------
In a nutshell :

make all
./assign features.txt

------------------------------------------------------

"make all" creates an executable "assign", run it with "./assign features.txt". An example of running on NERSC is with the command "qsub run". "features.txt" contains all the parameters, input files included


Compile tex file with luatex :

lualatex -shell-escape doc.tex
bibtex.doc


