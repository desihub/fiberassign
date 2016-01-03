# Fiber assignment code for DESI

## Contents

Several executables are available for fiber assignment.  The scripts for running the code are in the directory "scripts"

At present - January 2, 2016 - only the cori machine is available at NERSC.  The scripts for running the code on cori are as follows:

cori_fa:   This uses mocks for ELGs, LRGs, and QSOs, together with random mocks for standard stars and sky fibers, with densities chosen to be appropriate to the DESI experiment.  At the end, the yield of the various types is displayed.  The performance is controlled by an ascii file "fa_features."  In particular, you can choose have output giving the fiber assignments, and whether to have it as ascii or FITS. Unfortunately, this now takes a bit over 30 minutes. 

cori_shortrun_fa:  Same as cori_fa, except that there are only about one million targets (0<RA<10, -10<DEC<10).

cori_pipeline: This is a stripped down version of the full code.  It doesn't read a full galaxy collections but a reduced tareget file, which doesn't reveal the details, e.g. whether a QSO target is a  Lyman-alpha forest or not.  It also uses files for standard stars and sky fibers.

cori_run_mtl:  This generates a complete set of input files - Targ, SStars, SkyF, and Secret.  All these are derived from Martin White's mocks.  Targ = target files, SStars = standard stars, SkyF = skyfibers, Secret = information not available in pipeline, revealing true nature of each target

When edison is operative, the following scripts are appropriate.

run_fa:  This uses mocks for ELGs, LRGs, and QSOs, together with random mocks for standard stars and sky fibers, with densities chosen to be appropriate to the DESI experiment.  At the end, the yield of the various types is displayed.  The performance is controlled by an ascii file "fa_features."  In particular, you can choose have output giving the fiber assignments, and whether to have it as ascii or FITS.

shortrun_fa: Same as assign_fa, except that there are only about one million targets (0<RA<10, -10<DEC<10).

assign_pipeline This is a stripped down version of the full code.  It doesn't read a full galaxy collections but a reduced one, which doesn't reveal the details, e.g. whether a QSO target is a  Lyman-alpha forest or not.  The input is called the MTLfile and is specified in mtl_features.

run_mtl:  This generates a complete set of input files - Targ, SStars, SkyF, and Secret.  All these are derived from Martin White's mocks.  Targ = target files, SStars = standard stars, SkyF = skyfibers, Secret = information not available in pipeline, revealing true nature of each target

assign_mtl: This creates an MTLfile suitable for pipeline_fa

pipeline_fa: This is a stripped down version of the full code.  It doesn't read a full galaxy collections but a reduced tareget file, which doesn't reveal the details, e.g. whether a QSO target is a  Lyman-alpha forest or not.  It also uses files for standard stars and sky fibers.



##Building

By default, options for a generic Linux machine are used.  All 
dependencies (e.g. libcfitsio) are assumed to be installed in
the default compiler and linker search paths.  To build the 
software, go to the top source directory:

    $> cd fiberassign

To create the executable for run_fa,

    $> make clean_fa (optional)
    $> make install_fa

To create the executable for pipeline_fa

    $> make clean_pipeline (optional)
    $> make install_pipeline (optional)

To create the executable for run_mtl

    $> make clean_mtl
    $> make install_mtl

##The features files
fa_features: used by cori_fa, assign_fa
    specifies all the input files for the targets, positioners, tiles, etc.
    specifies output files, which need to be modified to the user's area

shortrun_features: used by cori_shortrun_fa, shortrun_fa
    specifies all the input files for the targets, positioners, tiles, etc.
    specifies output files, which need to be modified to the user's area

mtl_features: used to make target, standard star, skyfiber, secret files


pipeline_features: used by cori_pipeline
   

    
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

On cori


sbatch ./script/cori_fa 

or

sbatch ./script/cori_shortrun_fa

or 

sbatch ./script/cori_pipeline

