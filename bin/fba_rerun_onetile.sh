#!/bin/bash

usage () {
    cat <<HELP_USAGE

    Bash wrapper to rerun the fiber assignment of SV3 or Main BRIGHT/DARK tiles.
    It assumes the desi environment is already loaded, e.g.:
        source /global/cfs/cdirs/desi/software/desi_environment.sh master
    This script can:
    - generate the fiberassign intermediate products TILEID-{tiles,sky,targ,scnd,too}.fits
    - generate the fba-TILEID.fits(.gz) and/org fiberassign-TILEID.fits(.gz) files
    - compare the rerun fiberassign to the original one.

    `basename $0` ORIGFN RERUNDIR INTERMEDIATE FBA FASSIGN RUNCHECK INTERMEDIATE_FINAL

    ORIGFN        original fiberassign-TILEID.fits.gz file we want to rerun
    RERUNDIR      folder where outputs are produced (will be in RERUNDIR/ABC/)
    INTERMEDIATE  True or False (generate the intermediate products?)
    FBA           none, unzip, or zip (generate the fba-TILEID.fits(.gz) files if unzip or zip)
    FASSIGN       none, unzip, or zip (generate the fiberassign-TILEID.fits(.gz) files if unzip or zip)
    RUNCHECK      True or False (generate a *diff file)
    INTERMEDIATE_FINAL   - or folder to move intermediate products to
    OVERWRITE     True or False
HELP_USAGE
}

[ -z "$1" ] && { usage; exit;}

ORIGFN=$1
RERUNDIR=$2
INTERMEDIATE=$3
FBA=$4
FASSIGN=$5
RUNCHECK=$6
INTERMEDIATE_FINAL=$7
OVERWRITE=$8

echo "Start at: " `date`
echo ORIGFN       = $ORIGFN
echo RERUNDIR     = $RERUNDIR
echo INTERMEDIATE = $INTERMEDIATE
echo FBA          = $FBA
echo FASSIGN      = $FASSIGN
echo RUNCHECK     = $RUNCHECK
echo INTERMEDIATE_FINAL = $INTERMEDIATE_FINAL
echo OVERWRITE    = $OVERWRITE


# HACK for development
export PYTHONPATH=/global/homes/r/raichoor/software_dev/fiberassign_fba_rerun4/py:$PYTHONPATH
# HACK

# AR generate the *sh file(s)
echo "Generate *sh file(s)"
python -c 'from fiberassign.fba_launch_io import fba_rerun_fbascript; fba_rerun_fbascript("'$ORIGFN'", "'$RERUNDIR'", '$INTERMEDIATE', fba="'$FBA'", fiberassign="'$FASSIGN'", run_check='$RUNCHECK', intermediate_dir_final="'$INTERMEDIATE_FINAL'", overwrite='$OVERWRITE')'


if [[ "$INTERMEDIATE" != "-" ]]
then
    FN=`python -c 'from fiberassign.fba_launch_io import get_fba_rerun_scriptname; a = get_fba_rerun_scriptname("'$ORIGFN'", "'$RERUNDIR'", "intermediate"); print(a)'`
    echo "Execute $FN"
    $FN
fi

FN=`python -c 'from fiberassign.fba_launch_io import get_fba_rerun_scriptname; a = get_fba_rerun_scriptname("'$ORIGFN'", "'$RERUNDIR'", "fa"); print(a)'`
# HACK for development
export PYTHONPATH=`echo $PYTHONPATH | tr ":" "\n" | grep -v "/global/homes/r/raichoor/software_dev/fiberassign_fba_rerun4/py" | tr "\n" ":"`
# HACK
echo "Execute $FN"
$FN

echo "Done at: " `date`
