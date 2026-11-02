#!/usr/bin/env bash
# Regenerate the _static/*.inc files included by the .rst docs from each
# script's --help output.
#
# Requires a python environment with fiberassign's dependencies installed
# and the C++ extension built in place:
#
#   CC=g++ CXX=g++ CFLAGS="-fPIC" LDSHARED="g++ -shared" python setup.py build_ext --inplace
#   doc/update_help.sh
set -euo pipefail
cd "$(dirname "$0")"

export PYTHONPATH="$(cd ../py && pwd)${PYTHONPATH:+:$PYTHONPATH}"
export COLUMNS=80  # keep argparse line wrapping stable

scripts="fba_run fba_launch fba_rerun fba_fulltest fba_merge_results fba_run_qa fba_plot fba_plot_qa fiberassign"

for s in $scripts; do
    out="_static/${s}.inc"
    {
        echo ".. code-block:: console"
        echo
        ../bin/$s --help 2>/dev/null | awk '{ if (length($0) > 0) print "    " $0; else print "" }'
    } > "$out"
    echo "wrote $out"
done
