# Fiber assignment code for DESI

## Contents

Several executables are available for fiber assignment.  The scripts for running the code are in the directory "scripts"

At present - January 2, 2016 - only the cori machine is available at NERSC.  The scripts for running the code on cori are as follows:

cori_fa:   This uses mocks for ELGs, LRGs, and QSOs, together with random mocks for standard stars and sky fibers, with densities chosen to be appropriate to the DESI experiment.  At the end, the yield of the various types is displayed.  The performance is controlled by an ascii file "fa_features."  In particular, you can choose have output giving the fiber assignments, and whether to have it as ascii or FITS. Unfortunately, this now takes a bit over 30 minutes. 

cori_shortrun_fa:  Same as cori_fa, except that there are only about one million targets (0<RA<10, -10<DEC<10).

cori_pipeline: This is a stripped down version of the full code.  It doesn't read a full galaxy collections but a reduced one, which doesn't reveal the details, e.g. whether a QSO target is a  Lyman-alpha forest or not.  The input is called the MTLfile and is specified in mtl_features.

cori_run_mtl:  This generates a complete set of input files - Targ, SStars, SkyF, and Secret.  All these are derived from Martin White's mocks.  Targ = target files, SStars = standard stars, SkyF = skyfibers, Secret = information not available in pipeline, revealing true nature of each target

When edison is operative, the following scripts are appropriate.

run_fa:  This uses mocks for ELGs, LRGs, and QSOs, together with random mocks for standard stars and sky fibers, with densities chosen to be appropriate to the DESI experiment.  At the end, the yield of the various types is displayed.  The performance is controlled by an ascii file "fa_features."  In particular, you can choose have output giving the fiber assignments, and whether to have it as ascii or FITS.

shortrun_fa: Same as assign_fa, except that there are only about one million targets (0<RA<10, -10<DEC<10).

assign_pipeline This is a stripped down version of the full code.  It doesn't read a full galaxy collections but a reduced one, which doesn't reveal the details, e.g. whether a QSO target is a  Lyman-alpha forest or not.  The input is called the MTLfile and is specified in mtl_features.

cori_run_mtl:  This generates a complete set of input files - Targ, SStars, SkyF, and Secret.  All these are derived from Martin White's mocks.  Targ = target files, SStars = standard stars, SkyF = skyfibers, Secret = information not available in pipeline, revealing true nature of each target

run_mtl:  This generates a complete set of input files - Targ, SStars, SkyF, and Secret.  All these are derived from Martin White's mocks.  Targ = target files, SStars = standard stars, SkyF = skyfibers, Secret = information not available in pipeline, revealing true nature of each target

assign_mtl: This creates an MTLfile suitable for pipeline_fa

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




    pipeline_features
    mtl_features


    
are intended as templates.  You can make your own variations.


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


Some convenient scripts for running on NERSC are provided.  For example the script run_fa is

#PBS -S /bin/bash
#PBS -N Assign
#PBS -l mppwidth=24,walltime=1:29:59

#PBS -q premium
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=24
echo "# Running assign"

aprun -n 1 -N 1 -d 24 ./assign_fa fa_features.txt

Thus on edison, for example, from the fiberassign directory, you can submit a job as

qsub ./script/run_fa


# Features files
Fiber assignment is governed by features files.  For run_fa, the features file is fa_features. 
For the 200 sq deg version, shortrun_fa, there is shortrun_features.txt
For pipeline there is pipeline_features.

The intent is for the executables in master to work with their associated features files EXCEPT
that you should, if you want to write the fiber assignments, you must define an appropriate
outDir to which you can write.  Whether you do write is governed by the booleans PrintAscii and 
PrintFits in the features file.  If you are using your own input as MTLfile, you should set 
MaxSS and MaxSF to zero unless you have actually included appropriate numbers of standard stars and 
sky fibers.  Eventually, we will make regular collections of these available to augment any
target collections you develop.  


