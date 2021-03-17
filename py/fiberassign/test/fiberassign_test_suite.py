# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.test.fiberassign_test_suite
===================================================

Used to initialize the unit test framework via ``python setup.py test``.
"""
from __future__ import absolute_import, division, print_function

import sys

import unittest


def fiberassign_test_suite():
    """Returns unittest.TestSuite of desiutil tests.

    This is factored out separately from runtests() so that it can be used by
    ``python setup.py test``.
    """
    from os.path import dirname

    py_dir = dirname(dirname(__file__))
    return unittest.defaultTestLoader.discover(py_dir, top_level_dir=dirname(py_dir))


def runtests():
    """Run all tests in fiberassign.test.test_*."""
    # Load all TestCase classes from desispec/test/test_*.py
    tests = fiberassign_test_suite()
    # Run them and force exit with a non-zero process return value if they fail
    ret = unittest.TextTestRunner(verbosity=2).run(tests)
    if not ret.wasSuccessful():
        sys.exit(ret)


if __name__ == "__main__":
    runtests()
