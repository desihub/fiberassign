#!/bin/bash
#
# This script is intended to be run by desiInstall and only when installing
# branches. However, it can also be run by hand.
#
# The single command-line argument is the path to the Python executable
# that is calling desiInstall.
#
# As of 2022, setup.py in generally deprecated, so whatever replaces it
# can go in this script.
#
py=$1
${py} setup.py build_ext --inplace
