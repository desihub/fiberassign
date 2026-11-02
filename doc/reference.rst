.. _reference:

API Reference
========================

There are several key classes used by fiber assignment to load inputs, run the
calculation, and write outputs.  There are also some utility functions which
are useful on their own.


High Level Interface
-----------------------

The operations performed by the command-line scripts described in
:ref:`basictools` simply call functions from the high-level interface to parse
arguments and run the different steps.  Here we will duplicate the examples in
:ref:`basictools` but do these interactively from a python session rather than
calling scripts from the command line.  Each of the command line operations is
basically implemented by two functions:  one to parse the arguments and one to
run the code.  When running interactively, you can use a small helper function
to turn a dictionary of options into a list that can be passed to the parse
functions:

.. autofunction:: fiberassign.utils.option_list

For running the fiber assignment for a survey, we parse options with:

.. autofunction:: fiberassign.scripts.assign.parse_assign

And then run the survey with one of two functions depending on how we are doing
the assignment:

.. autofunction:: fiberassign.scripts.assign.run_assign_full

.. autofunction:: fiberassign.scripts.assign.run_assign_bytile

**EXAMPLE:** Run the full assignment interactively::

    from fiberassign.utils import option_list
    from fiberassign.scripts.assign import (parse_assign, run_assign_full)

    opts = {
        "targets": ["mtl.fits", "standards.fits", other.fits"],
        "write_all_targets": True,
        "dir": "out_raw"
    }
    optlist = option_list(opts)
    args = parse_assign(optlist)
    run_assign_full(args)

When merging results interactively, we can use the two functions:

.. autofunction:: fiberassign.scripts.merge.parse_merge

.. autofunction:: fiberassign.scripts.merge.run_merge

**EXAMPLE:** Run the output merging interactively::

    from fiberassign.utils import option_list
    from fiberassign.scripts.merge import (parse_merge,
                                           run_merge)
    opts = {
        "targets": ["mtl.fits", "standards.fits", other.fits"],
        "dir": "out_raw",
        "out": "out_merged"
    }
    optlist = option_list(opts)
    args = parse_merge(optlist)
    run_merge(args)

When running the QA interactively, we can use these functions:

.. autofunction:: fiberassign.scripts.qa.parse_qa

.. autofunction:: fiberassign.scripts.qa.run_qa

**EXAMPLE:** Run the output merging interactively::

    from fiberassign.utils import option_list
    from fiberassign.scripts.qa import (parse_qa, run_qa)

    opts = {
        "dir": "out_raw"
    }
    optlist = option_list(opts)
    args = parse_qa(optlist)
    run_qa(args)

When running the plotting interactively, these are the relevant functions:

.. autofunction:: fiberassign.scripts.plot.parse_plot

.. autofunction:: fiberassign.scripts.plot.run_plot

**EXAMPLE:** Run the output merging interactively::

    from fiberassign.utils import option_list
    from fiberassign.scripts.plot import (parse_plot, run_plot)

    opts = {
        "dir": "out_raw",
        "tiles": "plot_tiles.txt",
        "real_shapes": True
    }
    optlist = option_list(opts)
    args = parse_plot(optlist)
    run_plot(args)

And similarly when plotting the QA output:

.. autofunction:: fiberassign.scripts.qa_plot.parse_plot_qa

.. autofunction:: fiberassign.scripts.qa_plot.run_plot_qa

**EXAMPLE:** Run the output merging interactively::

    from fiberassign.utils import option_list
    from fiberassign.scripts.qa_plot import (parse_plot_qa,
                                             run_plot_qa)
    opts = {
        "qafile": "out_raw/qa.json"
    }
    optlist = option_list(opts)
    args = parse_plot_qa(optlist)
    run_plot_qa(args)


Hardware Properties
-----------------------

All of the physical properties of the DESI instrument are contained in the `Hardware` class.  Data from various files / sources must be read when constructing an instance of the Hardware class, so a helper function (`load_hardware`) is provided to make this easier.  In the future, there will likely be even more instrument details contained in the Hardware class and we may need to move to using a config file or some other mechanism to pass all this information.  Creating a hardware instance from Python should be done with the
``load_hardware`` function.  The full module API (exclusion zones, and the
CS5 / RA-DEC coordinate transforms) is documented below:

.. automodule:: fiberassign.hardware
    :members:

Once loaded, there are a variety of public methods available in the Hardware class:

.. autoclass:: fiberassign.hardware.Hardware
    :members:


Tile List (Survey Footprint)
---------------------------------

The properties of the current pointings / tiles in use during assignment is stored in a `Tiles` instance.  This just stores vectors of tile properties and the ordering of the tiles.  For convenience, the `load_tiles` function should be used from Python to construct a Tiles object.  This supports passing both a footprint (FITS format) file and also a list of tile IDs to select out of this larger footprint.

.. autofunction:: fiberassign.tiles.load_tiles

Once created, the tile properties can be accessed directly from the Tiles object:

.. autoclass:: fiberassign.tiles.Tiles
    :members:


Target Catalogs
------------------------------

Target catalogs are loaded one at a time and a subset of their data is appended
to a global dictionary of target properties.  Targets in a catalog are assigned
to one or more of the internal types (science, standard, sky, safe) based on
the DESI_TARGET (or other column) value and either default or custom bitmasks.

Every science target is assigned an integer priority value prior to fiber
assignment.  This integer value is typically the same for objects of the same
classification ("QSO", "LRG", "ELG", etc).  Sky targets effectively have a
priority of zero.  Each object is also assigned a random "subpriority"
value, which is a double precision number.  When two objects have the same
integer priority value, the subpriority is used to indicate which object has a
higher overall priority.  Sky targets are also given a random subpriority
value.

The target properties stored internally in fiberassign are found in the `Target` class for a single object:

.. autoclass:: fiberassign.targets.Target
    :members:

Before reading data, and empty `Targets` instance is created:

.. autoclass:: fiberassign.targets.Targets
    :members:

And then the `load_target_file` function is called one or more times to append target catalogs to the class instance.

After loading all the data, targets are indexed spatially with the C++ KD-tree (`src/kdTree.h`) for fast nearest-neighbor queries; a `scipy.spatial.KDTree` is only used for sky-source lookups (see `targets.py`).

Given a hardware model, our set of tiles, and this tree structure, we can now compute the target IDs available to every device location on every tile.  This is done with the `TargetsAvailable` class:

.. autoclass:: fiberassign.targets.TargetsAvailable
    :members:

We can also compute the inverse quantity:  the tile-location combinations that can reach each target ID.  This is done with the `LocationsAvailable` class:

.. autoclass:: fiberassign.targets.LocationsAvailable
    :members:

The full target-catalog module API (classification from DESI targeting
columns, the per-survey default bitmasks for main / SV1 / SV2 / SV3 / CMX,
and catalog loading) is:

.. automodule:: fiberassign.targets
    :members:


Assignment
---------------

The Assignment class is what drives the overall assignment process.  The algorithms of these functions are described below.  See the documentation for the individual methods for the mechanics of using them.

Assigning Targets to Unused Fibers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When assigning a particular class of target ("science", "standard", or "sky")
to unused locations, the same technique (and code) is used.

.. todo:: Algorithm discussion here, once relevant github issues and their associated changes are done.


Redistributing Science Targets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo:: Algorithm discussion here once relevant github issues and their associated changes are done.


Forced Assignment of Standards and Sky
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When calling this method, our goal is to displace science targets with either standards or sky targets in order to meet some required number per petal.  For each petal, we rank the available objects of interest (standards or sky) by total priority (priority + subpriority).  We also identify all science target priority values across the whole fiber assignment run and sort these from lowest to highest.

For each petal, we do the following pseudocode::

    for each science target priority "P" (lowest to highest)
        for each object (standard or sky) in total priority from high to low
            if object is reachable by devices on targets with priority "P"
                remove target and place positioner on object
                re-assign target to unused location on future tile if possible
                if enough objects on this petal
                    break

.. note:: This algorithm may be changed slightly in the near future.  See open github tickets for latest developments.


.. autoclass:: fiberassign.assign.Assignment
    :members:


Visualization and I/O
----------------------------

Output files are written, read and merged by the functions in
:mod:`fiberassign.assign` (including the top-level
:func:`fiberassign.assign.run` entry point), and plotting is done by
:mod:`fiberassign.vis` (used by the ``fba_plot`` and ``fba_plot_qa``
scripts):

.. automodule:: fiberassign.assign
    :members:

.. automodule:: fiberassign.vis
    :members:


Utilities
---------------

The fiberassign package includes a number of useful math and helper functions.

Geometric Shapes
~~~~~~~~~~~~~~~~~~~~~~~~

The geometry of the positioners is represented internally as a shape consisting of line segments and circles.

.. autoclass:: fiberassign.utils.Circle
    :members:

.. autoclass:: fiberassign.utils.Segments
    :members:

.. autoclass:: fiberassign.utils.Shape
    :members:


Logging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: fiberassign.utils.Logger
    :members:


Timing
~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: fiberassign.utils.Timer
    :members:

.. autoclass:: fiberassign.utils.GlobalTimers
    :members:


Rundate Helpers
~~~~~~~~~~~~~~~~~~~~~~~

These compare the rundate against the cutoff dates in
``py/fiberassign/data/cutoff-dates.yaml`` (see :ref:`usage`):

.. autofunction:: fiberassign.utils.get_date_cutoff

.. autofunction:: fiberassign.utils.get_obsdate

.. autofunction:: fiberassign.utils.get_fba_use_fabs

.. autofunction:: fiberassign.utils.get_whether_to_reorder_mtl

.. autofunction:: fiberassign.utils.get_whether_to_use_np_concatenate

.. autofunction:: fiberassign.utils.get_mjd

.. autofunction:: fiberassign.utils.get_default_static_obsdate

.. autofunction:: fiberassign.utils.assert_isoformat_utc

.. autofunction:: fiberassign.utils.get_obstheta_corr

.. autofunction:: fiberassign.utils.get_svn_version

.. autofunction:: fiberassign.utils.get_fa_gitversion

.. autofunction:: fiberassign.utils.get_last_line

.. autofunction:: fiberassign.utils.read_ecsv_keys

.. autofunction:: fiberassign.utils.get_revs_dates

.. autofunction:: fiberassign.utils.get_rev_fiberassign_changes


QA Metrics
----------

The per-tile and per-target metrics behind ``fba_run_qa``:

.. automodule:: fiberassign.qa
    :members:


GFA Targets
-----------

GFA / ETC target selection and flagging (the ``GUIDE_FLAG``, ``FOCUS_FLAG``
and ``ETC_FLAG`` bitmasks):

.. automodule:: fiberassign.gfa
    :members:


Stuck-Sky Handling
------------------

On-the-fly sky targets for stuck positioners (these count toward the
per-petal SKY requirement):

.. automodule:: fiberassign.stucksky
    :members:


Survey Launch, Patching, Reruns and Tertiary Targets
----------------------------------------------------

These modules back the ``fba_launch``, ``fba_patch``, ``fba_rerun`` and
``fba_tertiary_too`` entry points.  The module docstrings give an overview;
the key entry points are:

.. automodule:: fiberassign.fba_launch_io

.. autofunction:: fiberassign.fba_launch_io.launch_onetile_fa

.. autofunction:: fiberassign.fba_launch_io.assert_env_vars

.. autofunction:: fiberassign.fba_launch_io.create_tile

.. autofunction:: fiberassign.fba_launch_io.create_sky

.. autofunction:: fiberassign.fba_launch_io.create_mtl

.. autofunction:: fiberassign.fba_launch_io.create_too

.. autofunction:: fiberassign.fba_launch_io.update_fiberassign_header

.. automodule:: fiberassign.fba_patch_io

.. autofunction:: fiberassign.fba_patch_io.patch

.. autofunction:: fiberassign.fba_patch_io.diagnose_values_fafn

.. autofunction:: fiberassign.fba_patch_io.diagnose_columns_fafn

.. autofunction:: fiberassign.fba_patch_io.get_patching_params

.. automodule:: fiberassign.fba_rerun_io

.. autofunction:: fiberassign.fba_rerun_io.fba_rerun_fbascript

.. autofunction:: fiberassign.fba_rerun_io.fba_rerun_check

.. autofunction:: fiberassign.fba_rerun_io.get_fba_rerun_scriptname

.. automodule:: fiberassign.fba_tertiary_io

.. autofunction:: fiberassign.fba_tertiary_io.create_tertiary_too

.. autofunction:: fiberassign.fba_tertiary_io.get_targdir

.. autofunction:: fiberassign.fba_tertiary_io.read_targfn

.. autofunction:: fiberassign.fba_tertiary_io.read_priofn
