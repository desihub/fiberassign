.. _running:


Running
===============


The main executable ``fiberassign`` is a python wrapper around the C++
``fiberassign.exec`` code.  

Running ``fiberassign --help`` gives the full set of command line
options::

  usage: fiberassign [-h] --mtl MTL --sky SKY --stdstar STDSTAR --fibstatusfile FIBSTATUSFILE 
                   [--footprint FOOTPRINT]
                   [--positioners POSITIONERS] [--surveytiles SURVEYTILES]
                   [--telra TELRA] [--teldec TELDEC] [--tileid TILEID]
                   [--tileobsconditions TILEOBSCONDITIONS] [--outdir OUTDIR]
                   [--starmask STARMASK] [--rundate RUNDATE]
                   [--gfafile GFAFILE] [--nstarpetal NSTARPETAL]
                   [--nskypetal NSKYPETAL] [--nocleanup]

  optional arguments:
  -h, --help            show this help message and exit
  --mtl MTL             input targets (FITS file)
  --sky SKY             input sky positions (FITS file)
  --stdstar STDSTAR     input std stars (FITS file)
  --fibstatusfile FIBSTATUSFILE
                        list of positioners and its status (ECSV file)
  --footprint FOOTPRINT
                        list of tiles defining the footprint (FITS file)
  --positioners POSITIONERS
                        list of positioners on the focal plane (FITS file)
  --surveytiles SURVEYTILES
                        set of tiles to run fiberassign on (text file)
  --telra TELRA         Right Ascension of arbitrary pointing - overrides
                        --surveytiles
  --teldec TELDEC       Declination of arbitrary pointing - overrides
                        --surveytiles
  --tileid TILEID       Integer ID of arbitrary pointing - overrides
                        --surveytiles
  --tileobsconditions TILEOBSCONDITIONS
                        Mask describing observing program (DARK:1, GRAY:2,
                        BRIGHT:4) - overrides --surveytiles
  --outdir OUTDIR       output directory (default = ./)
  --starmask STARMASK   integer mask defining standard stars
  --rundate RUNDATE     run date [YYYY-MM-DD]
  --gfafile GFAFILE     GFA file (FITS tile)
  --nstarpetal NSTARPETAL
                        number of standard stars per petal (default=10)
  --nskypetal NSKYPETAL
                        number of sky fibers per petal (default=40)
  --nocleanup

   


An example that provides the minimal set of required arguments would be::

  fiberassign  --mtl mtl.fits --stdstar std.fits --sky sky.fits
  --fibstatusfile fiberstatus.ecsv --outdir $SCRATCH/temp/

In this example there are four files that must be **explicitly**
provided:

- ``mtl.fits``:  DESI Merged Target List files contain a single binary
  table covering the entire footprint. They contain the variables in
  the Targets files plus other variables that define the priority and
  number of observations as required by fiber assignment. These
  variables are computed using the available information both in the
  target and the DESI redshift catalogs. The full datamodel can be
  found `here
  <https://desidatamodel.readthedocs.io/en/latest/DESI_TARGET/mtl.html>`_.  

- ``std.fits``: DESI standard star locations contain a single binary
  table covering the entire footprint. They contain the variables in
  the Targets files with the objects that were flagged as
  standards. This file follows the `mtl` datamodel.

- ``sky.fits``:  DESI sky locations contain a single binary table
  covering the entire Legacy Surveys footprint. The imaging “blob 
  maps” are bisected to achieve a requisite number of sky locations
  per sq. deg. Sky locations are placed within the bisected grid as
  far from blobs that contain sources as is possible. Flux is measured
  in an aperture at each sky location. The full datamodel can be found
  `here
  <https://desidatamodel.readthedocs.io/en/latest/DESI_TARGET/skies.html>`_.  


- ``fiberstatus.ecsv``: 


