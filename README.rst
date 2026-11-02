===========
fiberassign
===========

|Actions Status| |Coveralls Status| |Documentation Status|

.. |Actions Status| image:: https://github.com/desihub/fiberassign/workflows/CI/badge.svg
    :target: https://github.com/desihub/fiberassign/actions
    :alt: GitHub Actions CI Status

.. |Coveralls Status| image:: https://coveralls.io/repos/desihub/fiberassign/badge.svg
    :target: https://coveralls.io/github/desihub/fiberassign
    :alt: Test Coverage Status

.. |Documentation Status| image:: https://readthedocs.org/projects/fiberassign/badge/?version=latest
    :target: https://fiberassign.readthedocs.io/en/latest/
    :alt: Documentation Status

Introduction
------------

This repository contains code to plan the location of DESI fiber positioners.  Given information about the instrument, a list of pointings (tiles), and catalogs of prioritized targets, fiber assignment computes the best positioner locations for each tile.

Installation
------------

To install a fixed version::

    $> pip install git+https://github.com/desihub/fiberassign.git@main

For development, clone the repository, compile the C++ extension in place,
and add the ``py`` and ``bin`` directories to your search paths::

    $> CC=g++ CXX=g++ CFLAGS="-fPIC" LDSHARED="g++ -shared" python setup.py build_ext --inplace
    $> export PYTHONPATH=${PWD}/py
    $> export PATH=${PWD}/bin:$PATH

See the `installation instructions <https://fiberassign.readthedocs.io/en/latest/install.html>`_
for the required DESI dependencies (desiutil, desimodel, desitarget,
desimeter), the desimodel data, and further details.

Versioning
----------

If you have tagged a version and wish to set the package version based on your
current git location::

    $> desi_update_version [-t TAG] fiberassign

And then install as usual.

Full Documentation
------------------

Please visit `fiberassign on Read the Docs`_

.. _`fiberassign on Read the Docs`: https://fiberassign.readthedocs.io/en/latest/

License
-------

fiberassign is free software licensed under a 3-clause BSD-style license. For details see
the ``LICENSE.rst`` file.
