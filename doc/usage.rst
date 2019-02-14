.. _usage:

Usage
===============

The fiberassign package has some main scripts as well as lower-level python tools that can be used in custom scripts.

Basic Tools
---------------------

There are several command-line programs which can be used to drive the fiber assignment process.  The first is the `fba_run` script which takes one or more target files and other options::

    %> fba_run --targets mtl.fits --targets standards.fits \
       --targets sky.fits --targets other_targets.fits \
       --outdir out_test

The names of the target files above are arbitrary and they can contain any targets.  The full list of targets is built up in memory from the contents of these files.  Running this command will produce an output directory of files (one per tile) containing the basic target information and assignment results.  For convenience, it is sometimes desired to copy extra columns of target information from the original input files into the output assignment files.  This can be done using a separate script::

    %> fba_merge_results --targets mtl.fits \
       --targets standards.fits \
       --targets sky.fits --targets other_targets.fits \
       --dir out_test

It is also frequently useful to plot the results of the assignment.  There are many customized plotting options possible using the low-level tools, but there
is also a command-line script to create a vector graphics (SVG) format plot of each tile.  Running this will require several minutes per tile, but multiple processes will be used to plot tiles in parallel::

    %> fba_plot_results --dir out_test

Some simple QA on the assignments can be run with::

    %> fba_run_qa --dir out_test

Which by default produces a JSON format named "qa.json" in the output directory.  To plot a simple sky representation of these results do::

    %> fba_plot_qa --qafile "out_test/qa.json"

If you are only plotting a few tiles and want to see the tile IDs on the plot, use the "--labels" option.


Interactive Debugging and Testing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fiberassign codebase respects the DESI_LOGLEVEL environment variable.  The default level if this is not set is "INFO".  Setting this environment variable to "DEBUG" will produce more detailed output.  There are other environment variables that can be set to dump even more details about the internal assignment process.  Examples::

    %> export DESI_DEBUG_TARGET=123456789
    %> export DESI_DEBUG_FIBER=4321
    %> export DESI_DEBUG_TILE=1111

These options are combined with a logical OR and any combination of tile, fiber, or target specified will have all possible info logged.

.. warning::
    Use of these "extra" debug variables will have a large impact on code
    performance.  Do not use in large production runs.


Legacy Compatibility Wrappers
---------------------------------------

In order to maintain backwards compatibility with the previous version of fiberassign, there is a replacement wrapper script called ``fiberassign`` which translates the command line options into arguments to pass to the newer scripts.

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
  found `here <https://desidatamodel.readthedocs.io/en/latest/DESI_TARGET/mtl.html>`_.

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
  `here <https://desidatamodel.readthedocs.io/en/latest/DESI_TARGET/skies.html>`_.

- ``fiberstatus.ecsv``:
