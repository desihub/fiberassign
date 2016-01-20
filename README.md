# Fiber assignment code for DESI

## Contents

Several executables are available for fiber assignment.  The scripts for running the code are in the directory "scripts"

The two primary executables are fiberassign and fiberassign_surveysim.  fiberassign is intended to function as part of the data flow, receiving its input from the merged target list.  fiberassign_surveysim simulates the full data flow.  In particular, by looking at a file containing the truth for every galaxy target, it updates the target list, requiring additional observations where necessary and eliminating them where not required.  fiberassign_surveysim provides a summary of how well it did in observing each category of target.

In addition to the executables, configuration files (currently called features files) are provided, as well as sample scripts for running the code at NERSC on cori or edison. The features files specify the galaxy collections and the tiles to be used. The scripts for running fiber assign are as follows:

do_fa:  This runs fiberassign with the configuration file fa_features.txt.

do_fa_surveysim:  This runs fiberassign_surveysim with the configuration file fa_surveysim_features.txt.

do_short_fa_surveysim:  The same as do_fa_surveysim, but runs on a very reduced galaxy sample from 200 sq. deg.  It uses shortrun_features.txt.

do_bright_fa_surveysim:  Again, runs fiberassign_surveysim, but with a bright time galaxy collection as specified in bright_time_features.


In addition, there are executables for creating the files of targets and of standard stars and sky fibers.  In particular

do_split_targ_secret takes a collection of galaxies and creates the target and "secret" files.  The target file contains the targets without disclosing whether, for example, a QSO target is a Ly-a, a target QSO (z<2.1), or a fake QSO.  This information is in the secret file, which is used in updating the target list and in evaluating the performance of the fiber assignment algorithm.

do_split_targ_secret_sky_std and do_short_split_targ_secret_sky_std shouldn't be needed.'

##Building

By default, options for a generic Linux machine are used.  All 
dependencies (e.g. libcfitsio) are assumed to be installed in
the default compiler and linker search paths.  To build the 
software, go to the top source directory:

    $> cd fiberassign

To create all the executables,

    $> make clean (optional)
    $> make install


##The features files
fa_features: used by do_fa
    specifies all the input files for the targets, positioners, tiles, etc.
    specifies output files, which need to be modified to the user's area

fa_surveysims_features: used by do_fa_surveysim
    specifies all the input files for the targets, positioners, tiles, etc.
    specifies output files, which need to be modified to the user's area

shortrun_features: used by do_short_fa_surveysim

bright_time_features: used by do_bright_fa_surveysim

    
These are intended as templates.  You can make your own variations.

##From Ted Kisner:

You can customize the build options by using one of the 
configurations in the platforms directory.  This is done by 
exporting an environment variable and setting it to the name of 
the file in the platforms directory you want to use.  For example, 
when building at NERSC with HPCPorts on the Cray systems, you would do:

    $> hpcports shared_gnu
    $> module load cfitsio-hpcp
    $> export PLATFORM=hpcports
    $> make install

Similarly, on OSX with the Clang compiler, you would do:

    $> export PLATFORM=osx
    $> make install

To install the "assign" executable to a different location,
you can set the INSTALL environment variable to the path you
wish to use.

## Running

from fiberassign on cori or edison

sbatch ./script/do_fa

sbatch ./script/do_fa_surveysim

etc.
