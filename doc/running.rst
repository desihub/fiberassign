.. _running:


Running
===============


The main executable ``fiberassign`` is a python wrapper around the
C++ ``fiberassign.exec`` code.  Run ``fiberassign --help`` to see the
full set of command line options.

An example would be:: 

  fiberassign  \
    --mtl targets/mtl.fits \
    --stdstar targets/standards-dark.fits \
    --sky targets/sky.fits \
    --surveytiles fiberassign/dark-tiles.txt \
    --fibstatusfile fiberstatus.ecsv \
    --outdir $SCRATCH/temp/


``doc/Guide_to_FiberAssignment.tex`` contains more details about underlying
algorithms but is out of date for the details of running fiberassign.
A pdf snapshot is available to DESI collaborators at
https://desi.lbl.gov/DocDB/cgi-bin/private/ShowDocument?docid=2742 .
