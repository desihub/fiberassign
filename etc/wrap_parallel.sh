#!/bin/bash

usage () {
    cat <<HELP_USAGE

    Run multiple tasks in parallel across several nodes with using GNU parallel.

    `basename $0` --input_fn INPUT_FN --njob NJOB

    --input_fn INPUT_FN  a text file with one bash command per line
    --njob     NJOB      nb of parallel jobs running

    This is from https://docs.nersc.gov/jobs/workflow/gnuparallel/#many-tasks-inside-a-multiple-node-allocation.
    The only minor adaptation is to include the two INPUT_FN and NJOB arguments.

    GNU parallel should be installed (e.g., if not already install, run "module load parallel").

    An example calling sequence is:
        srun --ntasks \$SLURM_NNODES --ntasks-per-node 1 --wait=0  `basename $0` --input_fn \$INPUT_FN --njob \$NJOB

    Should obviously not be run on a logging node.

HELP_USAGE
}

[ -z "$1" ] && { usage; exit;}

# AR read the provided arguments
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    for i in "$@"
    do
        case $i in
            --input_fn)
            INPUT_FN=$2
            shift
            shift
            ;;
            --njob)
            NJOB=$2
            shift
            shift
            ;;
            *)
                  # unknown option
            ;;
        esac
    done
    set -- "${POSITIONAL[@]}" # restore positional parameters
done

# AR print arguments
echo ""
echo "INPUT_FN="$INPUT_FN
echo "NJOB="$NJOB
echo ""


if [[ -z "${SLURM_NODEID}" ]]; then
    echo "need \$SLURM_NODEID set"
    exit
fi

if [[ -z "${SLURM_NNODES}" ]]; then
    echo "need \$SLURM_NNODES set"
    exit
fi

cat $INPUT_FN |                                            \
    awk -v NNODE="$SLURM_NNODES" -v NODEID="$SLURM_NODEID" \
    'NR % NNODE == NODEID' |                               \
    parallel -j $NJOB {}
