# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.scripts.qa_plot
==============================

High-level functions for plotting QA output.

"""
from __future__ import absolute_import, division, print_function

import os
import argparse
import json

from ..vis import plot_qa


def parse_plot_qa(optlist=None):
    """Parse QA plotting options.

    This parses either sys.argv or a list of strings passed in.  If passing
    an option list, you can create that more easily using the
    :func:`option_list` function.

    Args:
        optlist (list, optional): Optional list of arguments to parse instead
            of using sys.argv.

    Returns:
        (namespace):  an ArgumentParser namespace.

    """
    parser = argparse.ArgumentParser()

    parser.add_argument("--qafile", type=str, required=True, default=None,
                        help="Input QA file.")

    parser.add_argument("--outroot", type=str, required=False, default=None,
                        help="Output root file name.  Default uses input.")

    parser.add_argument("--labels", required=False, default=False,
                        action="store_true",
                        help="Plot tile IDs at center of circles.")

    args = None
    if optlist is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(optlist)

    # Check directory
    if not os.path.isfile(args.qafile):
        raise RuntimeError("Input file {} does not exist".format(args.qafile))

    if args.outroot is None:
        args.outroot = os.path.splitext(args.qafile)[0]

    return args


def run_plot_qa(args):
    """Run QA plotting.

    This uses the previously parsed options to read input data and make a
    plot of the QA results.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        None

    """
    qadata = None
    with open(args.qafile, "r") as f:
        qadata = json.load(f)
    plot_qa(qadata, args.outroot, labels=args.labels)
    return
