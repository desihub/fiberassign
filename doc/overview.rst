.. _overview:


Overview
==============

At the heart of DESI is the focal plane with robotically controlled
positioners capable of placing optical fibers on desired locations.
`fiberassignment` implements the algoritms that decide which fibers of
which should be assigned to observe which objects in series of
telescope pointings.   


The focal plane contains 5000 robotically controlled positioners with
a patrol radius of 6 mm and a positioning precision close to a few microns.
The positioners are arrayed in 10 petals shaped like pie slices. 
With about 10k pointings during the life of the survey, the system can
reach about 50 million targets. 
The instrumented area of the focal plane is 7.5 square degrees and
given the anticipated density of targets, about 25k targets fall on
the focal plane for each pointing (generally called a tile or
plate). With a total coverage of about 14,000 square degrees, the
average coverage of each point in the footprint is about 75/14=5.35. 

