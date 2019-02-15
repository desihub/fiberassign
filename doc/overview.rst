.. _overview:

Overview
==============

At DESI's focal plane there are robotically controlled positioners
capable of placing optical fibers on desired locations.
`fiberassign` implements the algorithms that decide which positioners
should be assigned to observe which objects in a series of telescope
pointings.

Constraints
-------------------

The focal plane has 5000 robotically controlled positioners with
a patrol radius of 6 mm and a positioning precision close to a few microns.
The positioners are arrayed in 10 petals shaped like pie slices.
With about 10000 pointings (also called tiles) during the life of the survey, the system can
reach about 50 million targets.
The instrumented area of the focal plane is 7.5 square degrees.
The anticipated density of targets is 25000 objets per squared degree.

Algorithm
--------------------

For each positioner in every tile choices must be made among the
potential targets.
Priorities are set in accordance with the scientific value of each
target class.
For the main (dark) survey program these classes include QSOs, LRGs, and
ELGs; for the secondary (bright) survey program these clases include
BGS (Bright Galaxies) and MWS (Milky Way Stars).
In addition, for both programs, fibers must be assigned to measure the
sky background (“sky fibers”) and calibration on standard stars.

The fiber assignment code needs a variety of inputs. The inputs
describe the hardware (i.e. locations of the positioners in the focal
plane), the observational strategy (a list of tiles), the targets
available for observation and calibration, together with their
relative observational priorities.

The outputs are a series of files describing the assignments.
There are individual files for each tile.
Each file contains the information about the target assigned to each
fiber, the list of of potential targets available to each fiber and
the list of targets for Guide/Focus/Alignment cameras.

History
--------

The first version of `fiberassign` was written by Martin White.
Later under the supervision of Bob Cahn it was expanded by Lile Wang,
Arthur Stril, Cyrille Doux, Aldo Riello, Louis Garrigue and Lucas
Pinol.  The most recent implementation was done by Ted Kisner.
Current development and maintenance is done by Stephen Bailey, Ted
Kisner and Jaime E. Forero-Romero.

A full list of contributors can be seen on the github repository:
https://github.com/desihub/fiberassign/graphs/contributors
