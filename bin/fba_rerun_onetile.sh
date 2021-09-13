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

    `basename $0` ORIGFN RERUNDIR INTERMEDIATE FBA FASSIGN RUNCHECK INTERMEDIATE_FINAL OVERWRITE FAVER_NOSWAP

    ORIGFN        original fiberassign-TILEID.fits.gz file we want to rerun
    RERUNDIR      folder where outputs are produced (will be in RERUNDIR/ABC/)
    INTERMEDIATE  True or False (generate the intermediate products?)
    FBA           none, unzip, or zip (generate the fba-TILEID.fits(.gz) files if unzip or zip)
    FASSIGN       none, unzip, or zip (generate the fiberassign-TILEID.fits(.gz) files if unzip or zip)
    RUNCHECK      True or False (generate a *diff file)
    INTERMEDIATE_FINAL   - or folder to move intermediate products to
    OVERWRITE     True or False
    FAVER_NOSWAP  True or False
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
FAVER_NOSWAP=$9

echo "Start at: " `date`
echo ORIGFN       = $ORIGFN
echo RERUNDIR     = $RERUNDIR
echo INTERMEDIATE = $INTERMEDIATE
echo FBA          = $FBA
echo FASSIGN      = $FASSIGN
echo RUNCHECK     = $RUNCHECK
echo INTERMEDIATE_FINAL = $INTERMEDIATE_FINAL
echo OVERWRITE    = $OVERWRITE
echo FAVER_NOSWAP = $FAVER_NOSWAP


# AR generate the *sh file(s)
echo "Generate *sh file(s)"
python -c 'from fiberassign.fba_rerun_io import fba_rerun_fbascript; fba_rerun_fbascript("'$ORIGFN'", "'$RERUNDIR'", '$INTERMEDIATE', outfba_type="'$FBA'", outfiberassign_type="'$FASSIGN'", run_check='$RUNCHECK', intermediate_dir_final="'$INTERMEDIATE_FINAL'", overwrite='$OVERWRITE', faver_noswap='$FAVER_NOSWAP')'


# AR rerun intermediate files, if requested
if [[ "$INTERMEDIATE" != "-" ]]
then
    FN=`python -c 'from fiberassign.fba_rerun_io import get_fba_rerun_scriptname; a = get_fba_rerun_scriptname("'$ORIGFN'", "'$RERUNDIR'", "intermediate"); print(a)'`
    echo "Execute $FN"
    $FN
fi


# AR rerun fa files
FN=`python -c 'from fiberassign.fba_rerun_io import get_fba_rerun_scriptname; a = get_fba_rerun_scriptname("'$ORIGFN'", "'$RERUNDIR'", "fa"); print(a)'`
echo "Execute $FN"
$FN

echo "Done at: " `date`
