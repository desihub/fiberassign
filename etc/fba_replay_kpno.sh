#!/bin/bash

NJOB=8
SUBFOLDER=

usage () {
    cat <<HELP_USAGE

    Rerun fiberassign for a batch of tiles designed at KPNO.

    `basename $0` --kpno_copydir KPNO_COPYDIR --rerundir RERUNDIR --njob NJOB

    --kpno_copydir  KPNO_COPYDIR  the local folder where KPNO files are copied into
    --subfolder     SUBFOLDER     optional; if set, only rerun that subfolder (e.g. 001, for 1??? tileids)
    --rerundir      RERUNDIR      the local folder where outputs will be stored
    --njob          NJOB          nb of parallel jobs running (default=8)
    --rerun                       optional; toggles whether to do the rerun or not
    --compare                     optional; toggles comparing rerun fiberassign output to original KPNO output
    --concat                      optional; toggles concatenating all the diff files into one large diff file

    Should be run on one or multiple interactive (or debug) nodes.
    If the reproducibility is successful, then the \$RERUNDIR/fba_nersc_kpno-diff.asc
        file should be empty, i.e. only with the header (and a file size of 139 octets).

    The code will overwrite any existing file in \$RERUNDIR.

    This script was typically used to perform NERSC vs. KPNO reproducibility
        when developing fiberassign code (one would then load the local, in development
        version of the fiberassign code before running).

    If any of the toggles (--rerun, --compare and --concat) are toggled on, any toggles not explciitly
        set on will be set to false. If none of the toggles are passed, they are all assumed to be true
        and all three steps are run.

    Files from KPNO should be organized as follows:
        \$KPNO_COPYDIR/000/fiberassign-001000.fits.gz

    \$RERUNDIR must exist.
    Rerun fiberassign files will be in:
        \$RERUNDIR/000/fiberassign-001000.fits.gz

    The NERSC vs. KPNO comparison summary file will be:
        \$RERUNDIR/fba_nersc_kpno-diff.asc

    The overall environment must be loaded.
    For example to rerun with the main environment, one should run beforehand (from a fresh environment):

=====
source /global/cfs/cdirs/desi/software/desi_environment.sh main
export DESIMODEL=/global/common/software/desi/perlmutter/desiconda/current/code/desimodel/main
export SKYHEALPIXS_DIR=\$DESI_ROOT/target/skyhealpixs/v1
====

    If one wants to rerun with a fiberassign branch in \$BRANCH_DIR, one should in addition run:
====
export PATH=\$BRANCH_DIR/bin:\$PATH
export PYTHONPATH=\$BRANCH_DIR/py:\$PYTHONPATH
====

HELP_USAGE
}

[ -z "$1" ] && { usage; exit;}

# Switches for each of the three steps, in case we want to run any combination
# of less than all three of them
RERUN=false
COMPARE=false
CONCAT=false

# AR read the provided arguments
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    for i in "$@"
    do
        case $i in
            --kpno_copydir)
            KPNO_COPYDIR=$2
            shift
            shift
            ;;
            --subfolder)
            SUBFOLDER=$2
            shift
            shift
            ;;
            --rerundir)
            RERUNDIR=$2
            shift
            shift
            ;;
            --njob)
            NJOB=$2
            shift
            shift
            ;;
            --rerun)
            RERUN=true
            ;;
            --compare)
            COMPARE=true
            ;;
            --concat)
            CONCAT=true
            ;;
            *)
                  # unknown option
            ;;
        esac
    done
    set -- "${POSITIONAL[@]}" # restore positional parameters
done

# If all are false still then none of them were passed, so set them all to true
# to run all steps.
if ! $RERUN && ! $COMPARE && ! $CONCAT
then
    RERUN=true
    COMPARE=true
    CONCAT=true
fi

# AR print arguments
echo ""
echo "KPNO_COPYDIR="$KPNO_COPYDIR
echo "SUBFOLDER="$SUBFOLDER
echo "RERUNDIR="$RERUNDIR
echo "NJOB="$NJOB
echo "RERUN="$RERUN
echo "COMPARE="$COMPARE
echo "CONCAT="$CONCAT
echo ""

# AR to manage multiple jobs over one or several nodes, with GNU parallel
WRAP_PARALLEL_FN=`which fba_launch | xargs dirname | awk '{print $1"/../etc/wrap_parallel.sh"}'`
echo "WRAP_PARALLEL_FN="$WRAP_PARALLEL_FN
echo ""

# AR tiles file
TILESFN=$DESI_SURVEYOPS/ops/tiles-main.ecsv

# AR files with the command lines for the GNU parallel call
CMDSH=$RERUNDIR/fba_nersc_kpno-cmd.sh
CMDLOG=$RERUNDIR/fba_nersc_kpno-cmd.log

# AR files the NERSC vs. KPNO comparison
DIFFSH=$RERUNDIR/fba_nersc_kpno-diff.sh
DIFFLOG=$RERUNDIR/fba_nersc_kpno-diff.log
DIFFASC=$RERUNDIR/fba_nersc_kpno-diff.asc

# AR get the fiberassign files to rerun.
if [[ "$SUBFOLDER" == "" ]]
then
    ORIGFNS=`ls $KPNO_COPYDIR/???/fiberassign-??????.fits.gz`
else
    ORIGFNS=`ls $KPNO_COPYDIR/$SUBFOLDER/fiberassign-??????.fits.gz`
fi

if $RERUN
then
    # AR text file with the commands to rerun the KPNO fiberassign files
    rm $CMDSH
    echo "rerun fiberassign-TILEID.fits.gz"
    echo "    Start at: " `date`
    for ORIGFN in $ORIGFNS
    do
        echo 'fba_rerun --infiberassign '$ORIGFN' --outdir '$RERUNDIR' --nosteps qa' >> $CMDSH
    done

    # AR execute the rerun in parallel
    CMD="srun --ntasks $SLURM_NNODES --ntasks-per-node 1 --wait=0 $WRAP_PARALLEL_FN --input_fn $CMDSH --njob $NJOB > $CMDLOG 2>&1"
    echo "    "$CMD
    eval $CMD
echo "    Done at: " `date`
fi

if $COMPARE
then
    # AR text file with the_NERSC vs. KPNO comparison command for each fiberassign file
    rm $DIFFSH
    echo "run NERSC vs. KPNO comparison"
    echo "    Start at: " `date`
    for ORIGFN in $ORIGFNS
    do
        SUBDIR=`basename $ORIGFN | awk '{print substr($1, 13, 3)}'`
        RERUNFN=$RERUNDIR/$SUBDIR/`basename $ORIGFN`
        DIFFFN=`echo $RERUNFN | sed -e 's/.fits.gz/.diff/'`
        echo 'python -c '\''from fiberassign.fba_rerun_io import fba_rerun_check; fba_rerun_check("'$ORIGFN'", "'$RERUNFN'", "'$DIFFFN'")'\''' >> $DIFFSH
    done

    # AR execute the comparison in parallel
    CMD="srun --ntasks $SLURM_NNODES --ntasks-per-node 1 --wait=0 $WRAP_PARALLEL_FN --input_fn $DIFFSH --njob $NJOB > $DIFFLOG 2>&1"
    echo "    "$CMD
    eval $CMD
    echo "    Done at: " `date`
fi

if $CONCAT
then
    # AR make a single diff file for the NERSC vs. KPNO comparison
    echo "create a single diff file"
    echo "    Start at: " `date`
    COUNT=0
    for ORIGFN in $ORIGFNS
    do
        SUBDIR=`basename $ORIGFN | awk '{print substr($1, 13, 3)}'`
        RERUNFN=$RERUNDIR/$SUBDIR/`basename $ORIGFN`
        DIFFFN=`echo $RERUNFN | sed -e 's/.fits.gz/.diff/'`
        if [[ $COUNT == 0 ]]
        then
            head -n 1 $DIFFFN > $DIFFASC
        fi
        COUNT=`echo $COUNT | awk '{print $1+1}'`
        cat $DIFFFN | grep -v \# >> $DIFFASC
    done
    echo "    Done at: " `date`
fi
