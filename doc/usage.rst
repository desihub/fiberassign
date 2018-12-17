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


Debugging and Testing
-----------------------------

The fiberassign codebase respects the DESI_LOGLEVEL environment variable.  The default level if this is not set is "INFO".  Setting this environment variable to "DEBUG" will produce more detailed output.  There are other environment variables that can be set to dump even more details about the internal assignment process.  Examples::

    %> export DESI_DEBUG_TARGET=123456789
    %> export DESI_DEBUG_FIBER=4321
    %> export DESI_DEBUG_TILE=1111

These options are combined with a logical OR and any combination of tile, fiber, or target specified will have all possible info logged.

.. warning::
    Use of these "extra" debug variables will have a large impact on code
    performance.  Do not use in large production runs.
