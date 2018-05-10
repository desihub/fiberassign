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
``platforms/`` and then run ``make install``;  *e.g.* on Cori Haswell::

    make PLATFORM=nersc_cori_haswell install

This will create the ``bin/fiberassign`` executable.

The version of fiberassign can be updated with::

    make version

or::

    make TAG=1.2.3 version

to set the version when tagging.

Running
-------

The main executable ``fiberassign`` takes its arguments (eight in total) through the command line.

An example would be

```
fiberassign  --mtl mtl.fits \
    --stdstar standards-dark.fits \
    --sky sky.fits \
    --surveytiles dark-tiles.txt \
    --footprint $DESIMODEL/data/footprint/desi-tiles.fits \
    --positioners $DESIMODEL/data/focalplane/fiberpos.txt \
    --fibstatusfile fiberstatus.ecsv \
    --outdir ./
```

``doc/Guide_to_FiberAssignment.tex`` contains more details.  A pdf snapshot
is available to DESI collaborators at
https://desi.lbl.gov/DocDB/cgi-bin/private/ShowDocument?docid=2742 .
