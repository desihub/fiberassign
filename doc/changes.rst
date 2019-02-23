.. _changes:

fiberassign change log
======================

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
