#!/bin/bash
#
# This script is intended to be run by desiInstall and only when installing
# branches. However, it can also be run by hand.
#
# The single command-line argument is the path to the Python executable
# that is calling desiInstall.
#
# This is meant to be run in ${INSTALL_DIR}, which is automatically
# defined by desiInstall.
#
# As of 2022, setup.py in generally deprecated, so whatever replaces it
# can go in this script.
#
py=$1
cd ${INSTALL_DIR}
${py} setup.py build_ext --inplace
