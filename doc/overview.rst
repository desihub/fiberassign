.. _overview:

Overview
========================

The purpose of this software is to determine the best configuration of fiber
positioners for each exposure.  The fiberassign code uses information about
the DESI instrument and applies that to lists of available targets with
pre-assigned prioritization.

Algorithm
---------------

The outermost command-line script basically goes through these steps:

#.  Load hardware information and list of tiles.

#.  Load all targets and categorize them as "science", "standards", or "sky".
    The "BAD_SKY" targets are treated as very low priority science targets.

#.  Compute available (reachable) targets for every fiber on every tile.

#.  Assign unused fibers to science targets (no maximum).

#.  Redistribute science targets.  If targets can be assigned to other tiles
    which have fewer science targets on the available petal, then move those.

#.  Assign unused fibers to standards up to some maximum per petal.

#.  Assign unused fibers to sky up to some maximum per petal.

#.  Forcibly assign standards to fibers up to some number per petal.  Bump
    low-priority science targets and attempt to reassign them.

#.  Forcibly assign sky to fibers up to some number per petal.  Bump
    low-priority science targets and attempt to reassign them.

#.  Assign unused fibers to sky (no maximum).

More details on the algorithms can be found in the docstrings for methods in
the Assignment class.
