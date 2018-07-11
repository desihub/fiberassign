===========
fiberassign
===========

Introduction
------------

This repository provides code for DESI fiber assignment, *i.e.* assigning
which fibers of which telescope pointings (tiles) should be assigned to
observe which objects in a target catalog.

Compiling
---------

To compile the code, set ``$PLATFORM`` to one of the recipes in the
``platforms/`` and then run ``make install``;  *e.g.* at NERSC::

    make PLATFORM=nersc_desiconda install

The version of fiberassign can be updated with::

    make version

or::

    make TAG=1.2.3 version

to set the version when tagging.

Running
-------

The main executable ``fiberassign`` is a python wrapper around the
C++ ``fiberassign_exec`` code.  Run ``fiberassign --help`` to see the
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
