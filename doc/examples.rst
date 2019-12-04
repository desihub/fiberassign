.. _examples:

Examples
===============

It is useful to consider several concrete examples of running fiber assignment and also describe the techniques for exploring the impact of the assignment algorithm.

Debugging Example
----------------------

This example shows a recent debugging test done with a small dataset.  This was
using several separate target files.  I ran them with::

    %> export DESI_LOGLEVEL=DEBUG
    %> fba_run --targets mtl.fits standards-dark.fits sky.fits

And the outputs were written by default to date-stamped output directory.  Then
I made some plots::

    %> fba_plot_results \
       --dir "out_fiberassign_2018-11-24T11:09:48"

I opened one resulting image (fiberassign_001148.pdf) and zoomed into one
petal.  I noticed that device location 2233 was assigned to a standard.  Let's
look at the details of how this device was assigned::

    %> export DESI_DEBUG_LOCATION=2233
    %> fba_run --targets mtl.fits standards-dark.fits sky.fits

The log showed that there were several targets available to this fiber on tile
1148::

    DEBUG: targets avail:  tile 1148, loc 2233 append ID \
        288230398599630147 (type=1), total priority 3200.07
    DEBUG: targets avail:  tile 1148, loc 2233 append ID \
        288230398599629415 (type=1), total priority 3000.91
    DEBUG: targets avail:  tile 1148, loc 2233 append ID \
        288230398599629372 (type=1), total priority 3000.24
    DEBUG: targets avail:  tile 1148, loc 2233 append ID \
        288230398599627959 (type=1), total priority 2000.04
    DEBUG: targets avail:  tile 1148, loc 2233 append ID \
        288230398599631299 (type=3), total priority 1500.79
    DEBUG: targets avail:  tile 1148, loc 2233 append ID \
        288230398599631166 (type=1), total priority 1500.1
    DEBUG: targets avail:  tile 1148, loc 2233 append ID \
        288230398599632448 (type=4), total priority 0.293402
    DEBUG: targets avail:  tile 1148, loc 2233 append ID \
        288230398599632236 (type=4), total priority 0.179244

Looking at the target type (which are the 4 categories of target used
internally in fiberassign, defined in targets.py / targets.h), we see that
there are 2 sky targets, 5 science targets, and one target which is both a
science target and a standard.  Later, during the assignment of unused
locations to science targets, we see::

    DEBUG: assign unused science: tile 1148, petal 4 loc 2233 \
        available target 288230398599630147, subpriority 0.0692912
    DEBUG: assign unused science: tile 1148, petal 4 loc 2233 \
        available target 288230398599629415, subpriority 0.914937
    DEBUG: assign unused science: tile 1148, petal 4 loc 2233 \
        available target 288230398599629372, subpriority 0.241329
    DEBUG: assign unused science: tile 1148, petal 4 loc 2233 \
        available target 288230398599627959, subpriority 0.0415658
    DEBUG: assign unused science: tile 1148, petal 4 loc 2233 \
        available target 288230398599631299, subpriority 0.786657
    DEBUG: assign unused science: tile 1148, petal 4 loc 2233 \
        available target 288230398599631166, subpriority 0.0962114
    DEBUG: assign unused science: tile 1148, petal 4 loc 2233 \
        available target 288230398599632448 is wrong type (4)
    DEBUG: assign unused science: tile 1148, petal 4 loc 2233 \
        available target 288230398599632236 is wrong type (4)
    DEBUG: find_best: tile 1148, loc 2233, target \
        288230398599630147, type 1 accept with priority = 3200, \
        subpriority = 0.0692912, obs_remain = 2
    DEBUG: find_best: tile 1148, loc 2233, target \
        288230398599630147, type 1 SELECTED
    DEBUG: assign unused science: tile 1148, petal 4 loc 2233 \
        found best object 288230398599630147

So it skipped over the two sky targets and selected a high-priority science
target.  Later, during the redistribution step, we see::

    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 considering for swap...
    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 considering tile indices 0 to 16070
    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 avail T/F 1148,2233 already assigned
    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 avail T/F 6926,175 already assigned
    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 avail T/F 12688,968 already assigned
    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 avail T/F 12688,991 already assigned
    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 avail T/F 12689,2561 already assigned
    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 avail T/F 12689,2587 already assigned
    DEBUG: redist: tile 1148, loc 2233, target \
        288230398599630147 keeping assignment

So the science target initially assigned to this location had no other
unassigned available tile / locations for potential swapping.  Further in the
assignment process, we see that a standard was found available to this
location.  Here is the log snippet from this step::

    DEBUG: assign force standard: tile 1148, petal 4, \
        loc 2233, found object 288230398599631299 \
        with weight 1500.79
    DEBUG: assign force standard: tile 1148, petal 4, \
        class 1500, object 288230398599631299, subpriority \
        1500.79, available loc 2233 at target \
        288230398599630147 is wrong class (3200)
    DEBUG: assign force standard: tile 1148, petal 4, \
        class 1600, object 288230398599631299, subpriority \
        1500.79, available loc 2233 at target \
        288230398599630147 is wrong class (3200)
    DEBUG: assign force standard: tile 1148, petal 4, \
        class 2000, object 288230398599631299, subpriority \
        1500.79, available loc 2233 at target \
        288230398599630147 is wrong class (3200)
    DEBUG: assign force standard: tile 1148, petal 4, \
        class 2100, object 288230398599631299, subpriority \
        1500.79, available loc 2233 at target \
        288230398599630147 is wrong class (3200)
    DEBUG: assign force standard: tile 1148, petal 4, \
        class 3000, object 288230398599631299, subpriority \
        1500.79, available loc 2233 at target \
        288230398599630147 is wrong class (3200)
    DEBUG: assign force standard: tile 1148, petal 4, \
        class 3200, object 288230398599631299, subpriority \
        1500.79, available loc 2233 bumping science \
        target 288230398599630147
    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 considering for swap...
    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 considering tile indices 0 to 16070
    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 avail T/F 1148,2233 already assigned
    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 avail T/F 6926,175 already assigned
    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 avail T/F 12688,968 already assigned
    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 avail T/F 12688,991 already assigned
    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 avail T/F 12689,2561 already assigned
    DEBUG: reassign: tile 1148, loc 2233, target \
        288230398599630147 avail T/F 12689,2587 not OK to assign

What happened here is that a standard was found to replace the low-priority
science target assigned to location 2233.  The existing science target was tested
for other available tile / locations, but all but one of those locations were already
assigned, and that one remaining device would produce a collision.  During the
forced assignment of sky targets, this is what happens to this location::

    DEBUG: assign force sky: tile 1148, petal 4, loc 2233, \
        found object 288230398599632448 with weight 0.293402
    DEBUG: assign force sky: tile 1148, petal 4, loc 2233, \
        found object 288230398599632236 with weight 0.179244
    DEBUG: assign force sky: tile 1148, petal 4, class 1500, \
        object 288230398599632448, subpriority 0.293402, \
        available loc 2233 at science target \
        288230398599631299 is also a standard- skipping
    DEBUG: assign force sky: tile 1148, petal 4, class 1500, \
        object 288230398599632236, subpriority 0.179244, \
        available loc 2233 at science target \
        288230398599631299 is also a standard- skipping
    DEBUG: assign force sky: tile 1148, petal 4, class 1600, \
        object 288230398599632448, subpriority 0.293402, available \
        loc 2233 at target 288230398599631299 is wrong class (1500)
    DEBUG: assign force sky: tile 1148, petal 4, class 1600, \
        object 288230398599632236, subpriority 0.179244, available \
        loc 2233 at target 288230398599631299 is wrong class (1500)
    DEBUG: assign force sky: tile 1148, petal 4, class 2000, \
        object 288230398599632448, subpriority 0.293402, available \
        loc 2233 at target 288230398599631299 is wrong class (1500)
    DEBUG: assign force sky: tile 1148, petal 4, class 2000, \
        object 288230398599632236, subpriority 0.179244, available \
        loc 2233 at target 288230398599631299 is wrong class (1500)
    DEBUG: assign force sky: tile 1148, petal 4, class 2100, \
        object 288230398599632448, subpriority 0.293402, available \
        loc 2233 at target 288230398599631299 is wrong class (1500)
    DEBUG: assign force sky: tile 1148, petal 4, class 2100, \
        object 288230398599632236, subpriority 0.179244, available \
        loc 2233 at target 288230398599631299 is wrong class (1500)

So for this device, the existing assignment was recognized as both a science
target and a standard, and was therefore not considered for bumping to place a
sky target.


Small Reference Run
--------------------------

This example is run on cori.nersc.gov, using data files in the project space
here::

    /project/projectdirs/desi/datachallenge/reference_runs/19.10/targets

After building (and optionally installing) fiberassign you should get an
interactive session on a compute node for up to 4 hours::

    %> salloc -N 1 -C haswell -A desi --qos=interactive -t 04:00:0

Once that job launches and you are on the compute node, set up some environment
variables::

    %> export OMP_NUM_THREADS=32
    %> export DESI_LOGLEVEL=DEBUG
    %> export \
       targetdir=/project/projectdirs/desi/datachallenge/reference_runs/19.10/targets

Now run the assignment using the default footprint tiling from
desimodel::

    %> time fba_run \
        --targets ${targetdir}/mtl-dark.fits \
        ${targetdir}/sky.fits \
        --dir out_ref_19.10 | tee log_ref_19.10

Make a plot of all tiles (you can also plot only some tiles or petals- see
options for fba_plot_results)::

    %> time fba_plot_results --dir out_ref_19.10

Merge all columns of the original target files into a new set of fiberassign
outputs::

    %> time fba_merge_results \
    --targets ${targetdir}/mtl-dark.fits \
    ${targetdir}/sky.fits --dir out_ref_19.10

Another debugging example
----------------------------

Examining the output from the above, in tile 001149 we find that location 6525 has a green positioner, meaning that this is both a science target and a standards target.  How did this assignment come about?  We can learn this by running again but with::

  %>export DESI_DEBUG_LOCATION=6525

  %>time fba_run     --targets ${targetdir}/mtl-dark.fits\
  ${targetdir}/sky.fits --dir out_debug2 |tee log_debug

 Then we can examine log_debug for the occurrence of 1149 and 6525::

   DEBUG: targets avail:  tile 1149, location 6525 append ID 288230398385719030 (type=1), total priority 3000.42
   DEBUG: targets avail:  tile 1149, location 6525 append ID 288230398385719269 (type=1), total priority 3000.37
   DEBUG: targets avail:  tile 1149, location 6525 append ID 288230398385719156 (type=1), total priority 3000.29
   DEBUG: targets avail:  tile 1149, location 6525 append ID 288230398385719119 (type=1), total priority 3000.12
   DEBUG: targets avail:  tile 1149, location 6525 append ID 288230398385719090 (type=1), total priority 3000.02
   DEBUG: targets avail:  tile 1149, location 6525 append ID 288230398385721102 (type=3), total priority 1500.58
   DEBUG: targets avail:  tile 1149, location 6525 append ID 288230398385720982 (type=1), total priority 1400.95
   DEBUG: targets avail:  tile 1149, location 6525 append ID 288230398385720952 (type=1), total priority 1400.32

The eventual choice has ID  288230398385721102, but this was not the initial assignment for tile 1149, location 6525 since it didn't have the highest priority.::

  DEBUG: target 288230398385719030 has available tile / location 1149, 6525 (fiber 3006)
  DEBUG: target 288230398385719269 has available tile / location 1149, 6525 (fiber 3006)
  DEBUG: target 288230398385719156 has available tile / location 1149, 6525 (fiber 3006)
  DEBUG: target 288230398385719119 has available tile / location 1149, 6525 (fiber 3006)
  DEBUG: target 288230398385719090 has available tile / location 1149, 6525 (fiber 3006)
  DEBUG: target 288230398385721102 has available tile / location 1149, 6525 (fiber 3006)
  DEBUG: target 288230398385720982 has available tile / location 1149, 6525 (fiber 3006)
  DEBUG: target 288230398385720952 has available tile / location 1149, 6525 (fiber 3006)
  DEBUG: assign unused science: tile 1149, petal 6 location 6525 (fiber 3006) available target 288230398385719030, priority 3000, subpriority 0.415054
  DEBUG: assign unused science: tile 1149, petal 6 location 6525 (fiber 3006) available target 288230398385719269, priority 3000, subpriority 0.370536
  DEBUG: assign unused science: tile 1149, petal 6 location 6525 (fiber 3006) available target 288230398385719156, priority 3000, subpriority 0.294639
  DEBUG: assign unused science: tile 1149, petal 6 location 6525 (fiber 3006) available target 288230398385719119, priority 3000, subpriority 0.121289
  DEBUG: assign unused science: tile 1149, petal 6 location 6525 (fiber 3006) available target 288230398385719090, priority 3000, subpriority 0.0176536
  DEBUG: assign unused science: tile 1149, petal 6 location 6525 (fiber 3006) available target 288230398385721102, priority 1500, subpriority 0.583245
  DEBUG: assign unused science: tile 1149, petal 6 location 6525 (fiber 3006) available target 288230398385720982, priority 1400, subpriority 0.951584DEBUG: assign unused science: tile 1149, petal 6 location 6525 (fiber 3006) available target 288230398385720952, priority 1400, subpriority 0.317493
  DEBUG: find_best: tile 1149, loc 6525, target 288230398385719030 science target with no remaining obs
  DEBUG: find_best: tile 1149, loc 6525, target 288230398385719269, type 1 accept with priority = 3000, subpriority = 0.370536, obsremain = 1
  DEBUG: find_best: tile 1149, loc 6525, target 288230398385719269, type 1 SELECTED
  DEBUG: assign unused science: tile 1149, petal 6 location 6525 (fiber 3006) found best object 288230398385719269
  DEBUG: assign_tileloc: tile 1149, loc 6525, target 288230398385719269, type 1 N_tile now = 4045 N_petal now = 456

Target 288230398385719269 is selected::

  DEBUG: reassign: tile 1149, location 6525, target 288230398385719269 considering for swap...
  DEBUG: reassign: tile 1149, location 6525, target 288230398385719269 considering tile indices 0 to 16070
  DEBUG: reassign: tile 1149, loc 6525, target 288230398385719269 avail T/F 1148,443 already considered or swapped
  DEBUG: reassign: tile 1149, loc 6525, target 288230398385719269 avail T/F 1149,6525 already assigned
  DEBUG: reassign: tile 1149, loc 6525, target 288230398385719269 avail T/F 6911,5272 already assigned
  DEBUG: reassign: tile 1149, loc 6525, target 288230398385719269 avail T/F 12673,4365 already assigned
  DEBUG: reassign: tile 1149, loc 6525, target 288230398385719269 avail T/F 18451,9252 already considered or swapped
  DEBUG: reassign: tile 1149, loc 6525, target 288230398385719269 avail T/F 24213,1165 already assigned
  DEBUG: redist: tile 1149, location 6525 (fiber 3006), target 288230398385719269 keeping assignment
  DEBUG: reassign: tile 5330, loc 6525, target 288230398276670339 avail T/F 1149,421 already considered or swapped
  DEBUG: reassign: tile 5330, loc 6525, target 288230398276670339 avail T/F 1149,398 already considered or swapped
  DEBUG: reassign: tile 11108, loc 6525, target 288230398377331091 avail T/F 1149,3222 already considered or swapped
  DEBUG: reassign: tile 16870, loc 6525, target 288230398381528045 avail T/F 1149,5051 already considered or swapped
  DEBUG: reassign: tile 16870, loc 6525, target 288230398381528045 avail T/F 1149,5059 already considered or swapped
  DEBUG: reassign: tile 28408, loc 6525, target 288230398377330190 avail T/F 1149,2219 already considered or swapped
  DEBUG: reassign: tile 34170, loc 6525, target 288230398381523339 avail T/F 1149,1003 already considered or swapped
  DEBUG: reassign: tile 39932, loc 6525, target 288230398289252223 avail T/F 1149,8141 already considered or swapped
  DEBUG: assign unused standard: tile 1149, petal 6 location 6525 (fiber 3006) is already assigned
  DEBUG: assign force standard: tile 1149, petal 6, location 6525 (fiber 3006), found object 288230398385721102 with weight 1500.58
  DEBUG: assign force standard: tile 1149, petal 6, class 1400, object 288230398385721102, total priority 1500.58, available loc 6525 (fiber 3006) at target 288230398385719269 is wrong class (3000)
  DEBUG: assign force standard: tile 1149, petal 6, class 1500, object 288230398385721102, total priority 1500.58, available loc 6525 (fiber 3006) at target 288230398385719269 is wrong class (3000)
  DEBUG: assign force standard: tile 1149, petal 6, class 1600, object 288230398385721102, total priority 1500.58, available loc 6525 (fiber 3006) at target 288230398385719269 is wrong class (3000)
  DEBUG: assign force standard: tile 1149, petal 6, class 2000, object 288230398385721102, total priority 1500.58, available loc 6525 (fiber 3006) at target 288230398385719269 is wrong class (3000)
  DEBUG: assign force standard: tile 1149, petal 6, class 2998, object 288230398385721102, total priority 1500.58, available loc 6525 (fiber 3006) at target 288230398385719269 is wrong class (3000)
  DEBUG: assign force standard: tile 1149, petal 6, class 3000, object 288230398385721102, total priority 1500.58, available l
  DEBUG: assign force standard: tile 1149, petal 6, class 3000, object 288230398385721102, total priority 1500.58, available loc 6525 (fiber 3006) at target 288230398385719269 (subpriority 0.370536) could be bumped
  DEBUG: assign force standard: tile 1149, petal 6, class 3000, object 288230398385721102, total priority 1500.58, available loc 6525 (fiber 3006) bumping science target 288230398385719269
  DEBUG: assign force standard: tile 1149, petal 6, class 3000, object 288230398385721102, total priority 1500.58, available loc 6525 (fiber 3006) bumping science target 288230398385719269
  DEBUG: reassign: tile 1149, location 6525, target 288230398385719269 considering for swap...
  DEBUG: reassign: tile 1149, location 6525, target 288230398385719269 considering tile indices 0 to 16070
  DEBUG: reassign: tile 1149, loc 6525, target 288230398385719269 avail T/F 1148,443 already assigned
  DEBUG: reassign: tile 1149, loc 6525, target 288230398385719269 avail T/F 1149,6525 already assigned
  DEBUG: reassign: tile 1149, loc 6525, target 288230398385719269 avail T/F 6911,5272 already assigned
  DEBUG: reassign: tile 1149, loc 6525, target 288230398385719269 avail T/F 12673,4365 already assigned
  DEBUG: reassign: tile 1149, loc 6525, target 288230398385719269 avail T/F 18451,9252 already assigned
  DEBUG: reassign: tile 1149, loc 6525, target 288230398385719269 avail T/F 24213,1165 already assigned
  DEBUG: unassign_tileloc: tile 1149, loc 6525, target 288230398385719269, type 1 N_tile now = 3966 N_petal now = 433
  DEBUG: assign_tileloc: tile 1149, loc 6525, target 288230398385721102, type 1 N_tile now = 3967 N_petal now = 434
  DEBUG: assign_tileloc: tile 1149, loc 6525, target 288230398385721102, type 2 N_tile now = 71 N_petal now = 10
  DEBUG: assign unused science: tile 1149, petal 6 location 6525 (fiber 3006) is already assigned
  DEBUG: assign unused standard: tile 1149, petal 6 location 6525 (fiber 3006) is already assigned

Finally, in search of a standard, target  288230398385719269 is unassigned and 288230398385721102 is used in its place.



Large Run
-----------------

This large DR7 example is run on cori.nersc.gov, using data files in the
project space here::

    /project/projectdirs/desi/target/fiberassign/dr7.1/0.10.3-dark

After building (and optionally installing) fiberassign you should get an
interactive session on a compute node for up to 4 hours::

    %> salloc -N 1 -C haswell -A desi --qos=interactive -t 04:00:0

Once that job launches and you are on the compute node, set up some environment
variables::

    %> export OMP_NUM_THREADS=32
    %> export DESI_LOGLEVEL=DEBUG
    %> export \
       targetdir=/project/projectdirs/desi/target/fiberassign/dr7.1/0.10.3-dark

Now run the assignment.  This will use about half of the RAM on a cori
haswell compute node and take about an hour- but half of that time is
writing the output files (something to work on)::

    %> time fba_run \
    --footprint ${targetdir}/input_tiles.fits \
    --targets ${targetdir}/mtl_large.fits \
    ${targetdir}/std_large.fits \
    ${targetdir}/sky_large.fits \
    --dir out_dr7.1_dark | tee log_dr7.1_dark

To save time for this example, only plot one of the petals on each tile::

    %> time fba_plot_results \
       --dir out_dr7.1_dark \
       --petals 4

Merge results::

    %> time fba_merge_results \
    --targets ${targetdir}/mtl_large.fits \
    ${targetdir}/std_large.fits \
    ${targetdir}/sky_large.fits \
    --dir out_dr7.1_dark
