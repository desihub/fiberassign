#!/bin/bash

# this script is a wrapper for a fba_launch call for one main tile
# there is no "svn up" here, no SVN checks
# the code will exit if a TILEID already existing in the official SVN folder exists

TILEID=$1 # e.g. 1000
QA=$2   # y or n

OUTDIR=$FA_HOLDING_PEN
TILESFN=$SURVEYOPS/ops/tiles-main.ecsv
DTCATVER=1.1.1

# setting the proper environment
# https://desi.lbl.gov/trac/wiki/SurveyOps/FiberAssignAtKPNO version #18
export DESI_PRODUCT_ROOT=/software/datasystems/desiconda/20200924
export DESI_ROOT=/data/datasystems
export DESI_TARGET=$DESI_ROOT/target
export DESI_SURVEYOPS=$DESI_ROOT/survey/ops/surveyops/trunk
module use $DESI_PRODUCT_ROOT/modulefiles
module load desiconda
module load desimodules/21.5
module swap desitarget/1.1.1
module swap fiberassign/5.0.0
module swap desimodel/master

# grabbing the RA, DEC, PROGRAM, STATUS, HA for TILEID
LINE=`awk '{if ($1 == '$TILEID') print $0}' $TILESFN`
RA=`echo $LINE | awk '{print $3}'`
DEC=`echo $LINE | awk '{print $4}'`
PROGRAM=`echo $LINE | awk '{print $5}'`
STATUS=`echo $LINE | awk '{print $8}'`
HA=`echo $LINE | awk '{print $10}'`

# small sanity check
if [[ "$STATUS" != "unobs" ]]
then
    echo "TILEID=$TILEID has already been observed; exiting"
    exit
fi

# calling fba_launch
if [[ "$QA" == "y" ]]
then
    CMD="fba_launch --outdir $OUTDIR --tileid $TILEID --tilera $RA --tiledec $DEC --survey main --program $PROGRAM --ha $HA --dtver $DTCATVER --steps qa --log-stdout --doclean y --svntiledir $FIBER_ASSIGN_DIR --worldreadable"
else
    CMD="fba_launch --outdir $OUTDIR --tileid $TILEID --tilera $RA --tiledec $DEC --survey main --program $PROGRAM --ha $HA --dtver $DTCATVER --nosteps qa --doclean n --svntiledir $FIBER_ASSIGN_DIR --log-stdout --worldreadable"
fi
echo $CMD
eval $CMD
