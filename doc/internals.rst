.. _internals:

Internals
========================

There are several key classes used by fiber assignment to load inputs, run the calculation, and write outputs.

Hardware Properties
-----------------------

.. autofunction:: fiberassign.hardware.load_hardware

.. autoclass:: fiberassign.hardware.Hardware
    :members:


Tile List (Survey Footprint)
---------------------------------

.. autofunction:: fiberassign.tiles.load_tiles

.. autoclass:: fiberassign.tiles.Tiles
    :members:


Target Catalogs
------------------------------

Target catalogs are loaded one at a time and appended to a global dictionary of target properties.  The bitmask used to ...

Every science target is assigned an integer priority value prior to fiber assignment.  This integer value is the same for objects of the same classification ("QSO", "LRG", "ELG", etc).  Sky targets effectively have a priority of zero.  Calibration standards may also be assigned a priority value.  If a calibration target is **also** a science target (which often means it is listed in multiple input target files), then this object is given the larger of its priority values.  Each object is also assigned a "subpriority" value, which is a double precision number.  When two objects have the same integer priority value, the subpriority is used to indicate which object has a higher overall priority.  Sky targets are also given a random subpriority value.

.. autofunction:: fiberassign.targets.load_target_file

.. autoclass:: fiberassign.targets.Target
    :members:

.. autoclass:: fiberassign.targets.Targets
    :members:

.. autoclass:: fiberassign.targets.TargetTree
    :members:

.. autoclass:: fiberassign.targets.TargetsAvailable
    :members:

.. autoclass:: fiberassign.targets.FibersAvailable
    :members:



Assignment
---------------

The fiberassign.Assignment class is what drives the overall assignment process.  The algorithms of these functions are described below.  See the documentation for the individual methods for the mechanics of using them.

.. autoclass:: fiberassign.assign.Assignment
    :members:


Assigning Targets to Unused Fibers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When assigning a particular class of target ("science", "standard", or "sky")
to unused fibers, the same technique (and code) is used.



Redistributing Science Targets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Forced Assignment of Standards and Sky
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When calling this method, our goal is to displace science targets with either standards or sky targets in order to meet some required number per petal.  For each petal, we rank the available objects of interest (standards or sky) by total priority (priority + subpriority).  We also identify all science target priority values across the whole fiber assignment run and sort these from lowest to highest.

For each petal, we do the following pseudocode::

    for each science target priority "P" (lowest to highest)
        for each object (standard or sky) in total priority from high to low
            if object is reachable by fibers on targets with priority "P"
                remove target and place fiber on object
                re-assign target to unused fiber on future tile if possible
                if enough objects on this petal
                    break
