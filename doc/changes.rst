fiberassign change log
======================

0.8.2 (unreleased)
------------------

* Standard star DESI_TARGET mask as input parameter (PR `#114`_)
* `fiberassign` is now a python wrapper around the C++ executable (PR `#116`_)
* Adds sky monitor fiber assignments (PR `#119`_)

.. _`#114`: https://github.com/desihub/fiberassign/pull/114
.. _`#116`: https://github.com/desihub/fiberassign/pull/116
.. _`#119`: https://github.com/desihub/fiberassign/pull/119


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
