.. _changes:

fiberassign change log
======================

1.2.1 (unreleased)
------------------

* No changes yet.

1.2.0 (2019-10-17)
------------------

* QA updates (PR `#216`_, `#230`_).
* Implement field rotation (PR `#219`_).
* Enforce sorting by fiber on output (PR `#223`_).
* fiberassign support for CMX targets + MAIN skies (PR `#224`_).
* Added cmx_science bits for first light targets (PR `#225`_).
* Use per-tile field rotations from desimodel.focalplane.fieldrot (PR `#226`_).
* Add GFA target quality cuts (PR `#227`_).
* Format updates to match ICS and some cleanup (PR `#228`_).

.. _`#216`: https://github.com/desihub/fiberassign/pull/216
.. _`#219`: https://github.com/desihub/fiberassign/pull/219
.. _`#223`: https://github.com/desihub/fiberassign/pull/223
.. _`#224`: https://github.com/desihub/fiberassign/pull/224
.. _`#225`: https://github.com/desihub/fiberassign/pull/225
.. _`#226`: https://github.com/desihub/fiberassign/pull/226
.. _`#227`: https://github.com/desihub/fiberassign/pull/227
.. _`#228`: https://github.com/desihub/fiberassign/pull/228
.. _`#230`: https://github.com/desihub/fiberassign/pull/230

1.1.0 (2019-09-25)
------------------

* Dynamic focalplane model (PR `#207`_).
* Add new bits to the cmx sciencemask and std mask (PR `#213`_).

.. _`#213`: https://github.com/desihub/fiberassign/pull/213
.. _`#207`: https://github.com/desihub/fiberassign/pull/207


1.0.4 (2019-06-24)
------------------

* Fix an issue with reproducibility of the ordering of available tile-fibers
  for each target (PR `#203`_).
* Switch to using device location (rather than fiber ID) as an indexing key
  throughout the code (PR `#204`_).
* Remove "short cut" when computing fiber collisions.  Always do the collision
  check (PR `#206`_).
* Restore sorting of output assignment in fiber ID order rather than device
  location (PR `#208`_).

.. _`#203`: https://github.com/desihub/fiberassign/pull/203
.. _`#204`: https://github.com/desihub/fiberassign/pull/204
.. _`#206`: https://github.com/desihub/fiberassign/pull/206
.. _`#208`: https://github.com/desihub/fiberassign/pull/208

1.0.3 (2019-05-30)
------------------

* PR `#202`_:

  * Gracefully allow fiberassign --stdstar to have duplicates with --mtl
  * Expose fba_run --sciencemask, --stdmask, etc. to fiberassign too
  * support fitsio 1.0.x
  * fix uninitialized variables bug

.. _`#202`: https://github.com/desihub/fiberassign/pull/202

1.0.1 (2019-05-13)
------------------

* Support different default masks for each program (PR `#193`_).
* Assign SAFE targets as backup if no SKY are available for sky monitor
  (PR `#191`_).
* Restored "safe" target type instead of just low priority science (PR `#189`_).
* Reorganized high-level code into package instead of script (PR `#188`_).

.. _`#188`: https://github.com/desihub/fiberassign/pull/188
.. _`#189`: https://github.com/desihub/fiberassign/pull/189
.. _`#191`: https://github.com/desihub/fiberassign/pull/191
.. _`#193`: https://github.com/desihub/fiberassign/pull/193

1.0.0 (2019-02-22)
------------------

* First tag of refactor/rewrite after merge (PR `#153`_).
* New C++ extension wrapped with pybind11.
* Python functions for I/O, visualization, QA.
* New commandline scripts for running assignment, merging input catalogs
  with output, making plots of outputs, etc.
* Overhaul of documentation.

.. _`#153`: https://github.com/desihub/fiberassign/pull/153

0.11.1 (2019-01-25)
-------------------

* Bug fix when using non-standard tiling (PR `#158`_).

.. _`#158`: https://github.com/desihub/fiberassign/pull/158

0.11.0 (2018-12-16)
-------------------

* Format updates to be closer to ICS fiberassign data model (PR `#157`_).
* Set `OBJTYPE='BAD'` and `DESI_TARGET=desi_mask.NO_TARGET` for broken, stuck,
  and unassigned fibers (PR `#154`_).
* Fix POTENTIAL target assignments HDU (broken in 0.10.2) (PR `#156`_).

.. _`#154`: https://github.com/desihub/fiberassign/pull/154
.. _`#156`: https://github.com/desihub/fiberassign/pull/156
.. _`#157`: https://github.com/desihub/fiberassign/pull/157

0.10.2 (2018-11-07)
-------------------

* Sort output by fiberid (PR `#147`_).
* Simplify required options (PR `#149`_).
* Add `--version` option (PR `#150`_).

.. _`#147`: https://github.com/desihub/fiberassign/pull/147
.. _`#149`: https://github.com/desihub/fiberassign/pull/149
.. _`#150`: https://github.com/desihub/fiberassign/pull/150

0.10.0 (2018-09-26)
-------------------

* Support both STD_FSTAR and STD bit names (PR `#139`_).
* Add more columns to output (PR `#141`_).
* Additional changes to try to match the data model (PR `#144`_).
* Fix collision calculation (PR `#146`_).

.. _`#139`: https://github.com/desihub/fiberassign/pull/139
.. _`#141`: https://github.com/desihub/fiberassign/pull/141
.. _`#144`: https://github.com/desihub/fiberassign/pull/144
.. _`#146`: https://github.com/desihub/fiberassign/pull/146


0.9.0 (2018-07-18)
------------------

* Standard star DESI_TARGET mask as input parameter (PR `#114`_).
* :command:`fiberassign` is now a python wrapper around the C++ executable (PR `#116`_).
* Adds sky monitor fiber assignments (PR `#119`_).
* Adds GFA targets HDU (PR `#122`_).
* Code format cleanup (PR `#123`_).
* Update build files; fix valgrind / compiler warnings (PR `#124`_).
* Bug fix: do not assume tileid is 5 digits long (PR `#126`_).
* Fixes sign flip in x,y <-> RA,dec conversions  (PR `#127`_).
* Checks for missing files (PR `#128`_).
* Fix unclosed file error (PR `#129`_).
* Bug fix: overflowing integer for SS flag (PR `#131`_).
* Show stuck/broken/unassigned fibers in :command:`qa-fiberassign` (PR `#132`_).

.. _`#114`: https://github.com/desihub/fiberassign/pull/114
.. _`#116`: https://github.com/desihub/fiberassign/pull/116
.. _`#119`: https://github.com/desihub/fiberassign/pull/119
.. _`#122`: https://github.com/desihub/fiberassign/pull/122
.. _`#123`: https://github.com/desihub/fiberassign/pull/123
.. _`#124`: https://github.com/desihub/fiberassign/pull/124
.. _`#126`: https://github.com/desihub/fiberassign/pull/126
.. _`#127`: https://github.com/desihub/fiberassign/pull/127
.. _`#128`: https://github.com/desihub/fiberassign/pull/128
.. _`#129`: https://github.com/desihub/fiberassign/pull/129
.. _`#131`: https://github.com/desihub/fiberassign/pull/131
.. _`#132`: https://github.com/desihub/fiberassign/pull/132

0.8.1 (2018-05-10)
------------------

* New FIBERMASK columns in fibermap files. (PR `#112`_).
* Computes RA+dec for unassigned, stuck, and broken fibers. (PR `#112`_).

.. _`#112`: https://github.com/desihub/fiberassign/pull/112


0.8.0 (2019-03-29)
------------------

* Clean up the command-line interface (PR `#105`_).
* Make fiberassign take more responsibility for installing itself (PR `#104`_).
* Allow fiberassign to report its version (PR `#104`_).

.. _`#105`: https://github.com/desihub/fiberassign/pull/105
.. _`#104`: https://github.com/desihub/fiberassign/pull/104

0.7.1 (2018-03-01)
------------------

* Fixed ``qa-fiberassign`` imports for desitarget 0.19.0 (PR `#102`_).

.. _`#102`: https://github.com/desihub/fiberassign/pull/102

0.7.0 (2018-02-23)
------------------

* Fill unassigned fibers with sky and stdstars if possible (PR `#100`_).
* Account for broken fibers and stuck positioners (PR `#101`_).

.. _`#101`: https://github.com/desihub/fiberassign/pull/101
.. _`#100`: https://github.com/desihub/fiberassign/pull/100

0.6.0 (2017-11-09)
------------------

* Guarantee that higher priority targets are placed first (PR `#84`_).
* Keep RA, Dec as double precision, not single precision (PR `#88`_).

.. _`#84`: https://github.com/desihub/fiberassign/pull/84
.. _`#88`: https://github.com/desihub/fiberassign/pull/88

0.5.3 (2017-09-30)
------------------

* ``bin/qa-fiberassign`` bug fixes.

0.5.2 (2017-09-30)
------------------

* Fixed indexing bug for ``LOCATION`` output.
* added WIP ``bin/qa-fiberassign``.
* Fixed missing collision checks (PR `#81`_).

.. _`#81`: https://github.com/desihub/fiberassign/pull/81

0.5.1 (2017-06-30)
------------------

* Reference tag.
