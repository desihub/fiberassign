.. _changes:

fiberassign change log
======================

5.4.1 (unreleased)
------------------

5.4.0 (2021-11-23)
------------------

This version now requires desitarget >= 2.2.1
Various updates, in particular to handle the Main/BACKUP program (but not only).

* Remove STD_WD from the per-petal STD counting (PR `# 415`_)
* Skyhealpixs gaia (PR `# 416`_)
* Gaia std (PR `# 417`_)
* Fba rerun5 (PR `# 420`_)
* Move some functions from fba_launch_io.py to utils.py (PR `# 419`_)
* change main/backup goaltime from 30 to 60 (PR `# 418`_)
* Populate EBV=0 rows in FIBERASSIGN with correct values (PR `# 414`_)

.. _`#416`: https://github.com/desihub/fiberassign/pull/415
.. _`#416`: https://github.com/desihub/fiberassign/pull/416
.. _`#417`: https://github.com/desihub/fiberassign/pull/417
.. _`#420`: https://github.com/desihub/fiberassign/pull/420
.. _`#419`: https://github.com/desihub/fiberassign/pull/419
.. _`#418`: https://github.com/desihub/fiberassign/pull/418
.. _`#414`: https://github.com/desihub/fiberassign/pull/414


5.3.0 (2021-09-20)
------------------

Allow several folders for inflating ledgers.

* allow several folders for inflating ledgers (PR `# 408`_)

.. _`#408`: https://github.com/desihub/fiberassign/pull/408

5.2.2 (2021-09-17)
------------------

Small fix to record the last change revision number in SVNDM and SVNMTL keywords.

* get_svn_version(): now grab 'Last Changed Rev' (PR `#407`_)

.. _`#407`: https://github.com/desihub/fiberassign/pull/407

5.2.1 (2021-09-17)
------------------

Small fixes related to fiberassign on-the-fly at KPNO and rerun at NERSC.

* Get SVN revision number in more compatible way (PR `#405`_)
* Fbarerun bugfix (PR `#406`_)

.. _`#405`: https://github.com/desihub/fiberassign/pull/405
.. _`#406`: https://github.com/desihub/fiberassign/pull/406

5.2.0 (2021-09-16)
------------------

Several changes during the 2021 summer shutdown:
- RUNDATE/MTLTIME/NOWTIME/ToO_MJD definitions to improve reproducibility
- functions for rerun
- propagate flux columns for secondaries
- bug fixes

* Try to remove files before overwriting them. (PR `#386`_)
* Also use touching-tile ledgers to set default MTLTIME (PR `#388`_)
* Set default RUNDATE to the latest TIME in the desi-state*.ecsv files (PR `#390`_)
* Add NOWTIME keyword in fiberassign-TILEID.fits.gz header (PR `#392`_)
* MTLTIME argument in create_too() (PR `#394`_)
* ToO: MJD definition and cut (PR `#396`_)
* Do not patch target tables from the gfa targets file (PR `#375`_)
* Avoid recomputing radec2xy when writing out results. (PR `#378`_)
* Fix duplicate negative TARGETID for stuck positioners (PR `#397`_)
* fba_cmx: remove GFAs when calling run_merge() (PR `#400`_)
* fba_rerun4 (PR `#398`_)
* update run_assign_bytile() with xycs5 (PR `#402`_)
* create_mtl(): add flux columns for secondaries (PR `#404`_)

.. _`#386`: https://github.com/desihub/fiberassign/pull/386
.. _`#388`: https://github.com/desihub/fiberassign/pull/388
.. _`#390`: https://github.com/desihub/fiberassign/pull/390
.. _`#392`: https://github.com/desihub/fiberassign/pull/392
.. _`#394`: https://github.com/desihub/fiberassign/pull/394
.. _`#396`: https://github.com/desihub/fiberassign/pull/396
.. _`#375`: https://github.com/desihub/fiberassign/pull/375
.. _`#378`: https://github.com/desihub/fiberassign/pull/378
.. _`#397`: https://github.com/desihub/fiberassign/pull/397
.. _`#400`: https://github.com/desihub/fiberassign/pull/400
.. _`#398`: https://github.com/desihub/fiberassign/pull/398
.. _`#402`: https://github.com/desihub/fiberassign/pull/402
.. _`#404`: https://github.com/desihub/fiberassign/pull/404


5.1.1 (2021-07-09)
------------------

No algorithmic changes, bug fixes.

* fix case if no ToO targets selected, for mv_temp2final() (PR `#382`_)
* Proper handling of --worldreadable in fba_rerun (PR `#383`_)

.. _`#382`: https://github.com/desihub/fiberassign/pull/382
.. _`#383`: https://github.com/desihub/fiberassign/pull/383

5.1.0 (2021-07-08)
------------------

Changes to ease fiberassign on-the-fly, and fba_rerun script.

* A couple of speed-ups (stuck-sky, hardware-loading) (PR `#373`_)
* fba_rerun script (PR `#376`_)
* Add features to fba_launch to support fiberassign on the fly (PR `#380`_)
* Quickread2: taking advantage of desitarget/1.2.2 speed-ups (PR `#381`_)

.. _`#373`: https://github.com/desihub/fiberassign/pull/373
.. _`#376`: https://github.com/desihub/fiberassign/pull/376
.. _`#380`: https://github.com/desihub/fiberassign/pull/380
.. _`#381`: https://github.com/desihub/fiberassign/pull/381

5.0.0 (2021-05-29)
------------------

Algorithmic changes to not change SUBPRIORITY when running fiberassign.

* Refactor internal dataflow for PLATE_RA/PLATE_DEC, without external
  changes to outputs except different POTENTIAL_TARGETS row order (PR `#353`_).
* Include desimeter in DEPNAMnn/DEPVERnn keywords (PR `#364`_).
* Don't override SUBPRIORITY while preparing files (PR `#366`_).
* Add Gaia-based variability bit 5 to ETC_FLAG (PR `#367`_).
* fba_launch options to run/exclude specific steps (PR `#368`_).
* use desitarget.gaiamatch.gaia_psflike() for PSF-like criterion (PR `#369`_).
* use np.nan_to_num() to avoid warnings: RuntimeWarning: invalid value encountered in greater (PR `#370`_).
* Only interpret exclusion regions on demand (PR `#371`_).
* adding desimeter path, version in log (PR `#372`_).

.. _`#353`: https://github.com/desihub/fiberassign/pull/353
.. _`#364`: https://github.com/desihub/fiberassign/pull/364
.. _`#366`: https://github.com/desihub/fiberassign/pull/366
.. _`#367`: https://github.com/desihub/fiberassign/pull/367
.. _`#368`: https://github.com/desihub/fiberassign/pull/368
.. _`#369`: https://github.com/desihub/fiberassign/pull/369
.. _`#370`: https://github.com/desihub/fiberassign/pull/370
.. _`#371`: https://github.com/desihub/fiberassign/pull/371
.. _`#372`: https://github.com/desihub/fiberassign/pull/372

4.0.1 (2021-05-18)
------------------

No algorithmic changes.

* Add timeout to wget fetch of imaging cutout for QA (PR `#361`_).

.. _`#361`: https://github.com/desihub/fiberassign/pull/361

4.0.0 (2021-05-14)
------------------

First release used for main survey observations.

Note: the format changed to add PLATE_RA, PLATE_DEC output columns, thus
bumping the major version number even though the results are algorithmically
identical to 3.0.0.

* Robust if target-of-opportunity (ToO) inputs don't exist (PR `#352`_).
* Don't set $SKYBRICKS_DIR in module file (desitarget does that now) (direct commit).
* Add PLATE_RA, PLATE_DEC columns while merging as placeholders for future
  chromatic offsets use (PR `#355`_).
* Add ``fba_launch --hdr_survey X --hdr_faprgrm Y`` options, defaulting to
  ``--survey`` and ``--program`` (PR `#356`_).

.. _`#352`: https://github.com/desihub/fiberassign/pull/352
.. _`#355`: https://github.com/desihub/fiberassign/pull/355
.. _`#356`: https://github.com/desihub/fiberassign/pull/356

3.0.0 (2021-05-13)
------------------

Major update to use desimeter for x,y <-> ra,dec transforms to include
airmass and ADC distortions.
Requires desimeter >= 3.6.5 and desitarget >= 1.0.0 .

* Simplify and improve ``bin/fba_plot`` (PR `#336`_).
* Use Gaia-based FLUX_R for GFA_TARGETS extension to avoid DR9 saturation
  (PR `#344`_).
* Record $DESI_SUREYOPS/mtl and $DESIMODEL/data svn revision numbers in
  output header keywords SVNMTL and SVNDM (PR `#346`_).
* Add inner exclusion ``|R1+R2|+100um`` (commits `01206c1`_ and `6e78851`_)
* ``fba_launch --mtltile`` default to latest timestamp in mtl file (PR `#347`_).
* ``fba_launch`` add support for main survey inputs (PR `#349`_).
* Use desimeter for x,y <-> ra,dec transforms (PR `#348`_).
* Expand default positioner polygons by 50 microns and edges by 400 microns,
  adjustable with options (PR `#350`_).

.. _`#336`: https://github.com/desihub/fiberassign/pull/336
.. _`#346`: https://github.com/desihub/fiberassign/pull/346
.. _`01206c1`: https://github.com/desihub/fiberassign/commit/01206c14d397df3e7901220257b826c721a66762
.. _`6e78851`: https://github.com/desihub/fiberassign/commit/6e78851160ebe10a172f5121391121c78242306f
.. _`#344`: https://github.com/desihub/fiberassign/pull/344
.. _`#347`: https://github.com/desihub/fiberassign/pull/347
.. _`#348`: https://github.com/desihub/fiberassign/pull/348
.. _`#349`: https://github.com/desihub/fiberassign/pull/349
.. _`#350`: https://github.com/desihub/fiberassign/pull/350

2.5.1 (2021-05-11)
------------------

* Adds ``bin/fba_launch_dc3r2_gama`` to support a special tile (PR `#345`_).

.. _`#345`: https://github.com/desihub/fiberassign/pull/345

2.5.0 (2021-05-11)
------------------

* Major refactor of ``bin/fba_launch`` into functions in
  ``fiberassign.fba_launch_io`` for reuse by other scripts (PR `#343`_).
* Headers record skybricks input version; support skybricks/v3 format
  (PR `#341`_, `#342`_).

.. _`#341`: https://github.com/desihub/fiberassign/pull/341
.. _`#342`: https://github.com/desihub/fiberassign/pull/342
.. _`#343`: https://github.com/desihub/fiberassign/pull/343

2.4.0 (2021-05-04)
------------------

* ``fba_launch --isodate`` option to set timestamp for MTL ledger reading
  (PR `#334`_).
* Assign stuck positioners to sky if possible, using skybricks/v2 lookup
  (PR `#337`_).
* Add per-slitblock sky fiber limits (PR `#338`_).
* Report counts of assigned fibers as fiberassign proceeds (PR `#339`_).
* Apply theta-phi offsets when computing locs of stuck positioners;
  fixes NaNs in outputs (PR `#340`_).
* Park unassigned positioners at phi=150 instead of 180. (PR `#340`_).

.. _`#334`: https://github.com/desihub/fiberassign/pull/334
.. _`#337`: https://github.com/desihub/fiberassign/pull/337
.. _`#338`: https://github.com/desihub/fiberassign/pull/338
.. _`#339`: https://github.com/desihub/fiberassign/pull/339
.. _`#340`: https://github.com/desihub/fiberassign/pull/340

2.3.0 (2021-04-22)
------------------

First used for tiles 98,179,198,209,231,287,315,375,423,438,441
on 2021-04-22 before making tag.

* Change assignment strategy of leftover fibers (PR `#321`_).
* Use UTC time everywhere (PR `#327`_, `#328`_).

.. _`#321`: https://github.com/desihub/fiberassign/pull/321
.. _`#327`: https://github.com/desihub/fiberassign/pull/327
.. _`#328`: https://github.com/desihub/fiberassign/pull/328

2.2.0 (2021-03-31)
------------------

* Support dedicated secondary programs (PR `#311`_).
* Migrate from Travis to GitHub workflows (PR `#313`_).
* Support sv1 tiles (PR `#314`_).
* Support new desimodel focal plane state format (PR `#315`_).
* Support sv2 tiles (PR `#318`_).
* new fba_launch wrapper script (PR `#319`_).
* Support matplotlib 3.3.4 (PR `#320`_).
* use desitarget write_skies instead of write_targets for skies
  (commit dd69bdd)

.. _`#311`: https://github.com/desihub/fiberassign/pull/311
.. _`#313`: https://github.com/desihub/fiberassign/pull/313
.. _`#314`: https://github.com/desihub/fiberassign/pull/314
.. _`#315`: https://github.com/desihub/fiberassign/pull/315
.. _`#318`: https://github.com/desihub/fiberassign/pull/318
.. _`#319`: https://github.com/desihub/fiberassign/pull/319
.. _`#320`: https://github.com/desihub/fiberassign/pull/320

2.1.1 (2021-02-11)
------------------

* Added bin/sv1-summary.py (PR `#301`_, `#308`_).
* Updates for secondary target support (PR `#303`_).
* Orion Rosette Praesepe support (PR `#306`_).
* Remove unnecessary (incorrect) -Wstrict-prototypes compile flag (PR `#309`_).

.. _`#301`: https://github.com/desihub/fiberassign/pull/301
.. _`#303`: https://github.com/desihub/fiberassign/pull/303
.. _`#306`: https://github.com/desihub/fiberassign/pull/306
.. _`#308`: https://github.com/desihub/fiberassign/pull/308
.. _`#309`: https://github.com/desihub/fiberassign/pull/309


2.1.0 (2020-12-23)
------------------

Major script and format updates for SV1 in December 2020.

* Add `SV1_*_TARGET` columns (PR `#287`_).
* fba_cmx gzip output (PR `#288`_).
* Add fba_sv1 script (PR `#289`_, `#291`_, `#293`_, `#294`_, `#299`_).
* Use read_targets_in_tiles quick=True option (PR `#290`_).
* Option for specifying proper motion epoch --pmtime (PR `#295`_).
* Update default fiberassign columns (PR `#297`_, `#298`_).

.. _`#287`: https://github.com/desihub/fiberassign/pull/287
.. _`#288`: https://github.com/desihub/fiberassign/pull/288
.. _`#289`: https://github.com/desihub/fiberassign/pull/289
.. _`#290`: https://github.com/desihub/fiberassign/pull/290
.. _`#291`: https://github.com/desihub/fiberassign/pull/291
.. _`#293`: https://github.com/desihub/fiberassign/pull/293
.. _`#294`: https://github.com/desihub/fiberassign/pull/294
.. _`#295`: https://github.com/desihub/fiberassign/pull/295
.. _`#297`: https://github.com/desihub/fiberassign/pull/297
.. _`#298`: https://github.com/desihub/fiberassign/pull/298
.. _`#299`: https://github.com/desihub/fiberassign/pull/299

2.0.0 (2020-12-11)
------------------

NOTE: New major version number due to fiberassign format changes.

* Added fba_cmx script for commissioning
  (PR `#277`_, `#280`_, `#281`_, `#283`_, `#286`_).
* Reduces the number of target columns propagated into the fiberassign
  file (PR `#279`_)
* Add SUPP_SKY targets to OBJTYPE=SKY (PR `#282`_).

.. _`#277`: https://github.com/desihub/fiberassign/pull/277
.. _`#279`: https://github.com/desihub/fiberassign/pull/279
.. _`#280`: https://github.com/desihub/fiberassign/pull/280
.. _`#281`: https://github.com/desihub/fiberassign/pull/281
.. _`#282`: https://github.com/desihub/fiberassign/pull/282
.. _`#283`: https://github.com/desihub/fiberassign/pull/283
.. _`#286`: https://github.com/desihub/fiberassign/pull/286

1.4.2 (2020-10-02)
------------------

* Support C++11, not requiring C++14 (PR `#273`_).

.. _`#273`: https://github.com/desihub/fiberassign/pull/273

1.4.1 (2020-08-04)
------------------

* Fix tests and qa-fiberassign (PR `#269`_).
* Simplify handling of MWS secondary bits in creating sv1_sciencemask (PR `#268`_).
* Fix bug in the range checking of positioner theta / phi angles (PR `#267`_).
* Move the checks for positioner reachability from the assignment code to the
  TargetsAvailable class (PR `#264`_).
* Use a specific rundate for unit tests, to ensure consistent focalplane
  model (PR `#262`_).

.. _`#262`: https://github.com/desihub/fiberassign/pull/262
.. _`#264`: https://github.com/desihub/fiberassign/pull/264
.. _`#267`: https://github.com/desihub/fiberassign/pull/267
.. _`#268`: https://github.com/desihub/fiberassign/pull/268
.. _`#269`: https://github.com/desihub/fiberassign/pull/269
  
1.4.0 (2020-03-19)
------------------

* Change assignment algorithm to be based on target order instead of
  fiber order (PR `#258`_).
* Fix radial platescale interpolation to work with latest desimodel (PR `#259`_).

.. _`#258`: https://github.com/desihub/fiberassign/pull/258
.. _`#259`: https://github.com/desihub/fiberassign/pull/259

1.3.1 (2020-03-13)
------------------

* Change targets to correctly look up desi and secondary mask (PR `#250`_).
* Add minisv2 bits (PR `#252`_).
* Extended QA (PR `#253`_).
* Avoid propagation of 2D target columns into FIBERASSIGN and TARGETS HDU (PR `#255`_)
* Increase target realism in unit tests (PR `#256`_)
* New SV0 science target bits from desitarget/0.37.0 (PR `#257`_)

.. _`#250`: https://github.com/desihub/fiberassign/pull/250
.. _`#252`: https://github.com/desihub/fiberassign/pull/252
.. _`#253`: https://github.com/desihub/fiberassign/pull/253
.. _`#255`: https://github.com/desihub/fiberassign/pull/255
.. _`#256`: https://github.com/desihub/fiberassign/pull/256
.. _`#257`: https://github.com/desihub/fiberassign/pull/257

1.3.0 (2019-12-20)
------------------

* Change output filenames to fba-*.fits and fiberassign-*.fits (PR `#235`_).
* Propagate run date/teim and depencency versions to outputs (PR `#240`_).
* Update documentation to more recent data releases (PR `#242`_).

.. _`#235`: https://github.com/desihub/fiberassign/pull/235
.. _`#240`: https://github.com/desihub/fiberassign/pull/240
.. _`#242`: https://github.com/desihub/fiberassign/pull/242

1.2.1 (2019-10-31)
------------------

* Implement GFA and petal boundary exclusion zones (PR `#233`_).
* Plot GFA and petal keepouts for all petals, not just petal zero (PR `#234`_).

.. _`#233`: https://github.com/desihub/fiberassign/pull/233
.. _`#234`: https://github.com/desihub/fiberassign/pull/234

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
