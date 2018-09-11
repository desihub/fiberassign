.. _install:


Installation
===============

The DESI spectroscopic pipeline requires and interfaces with other external software packages.  This document assumes that you are setting up a "development" software stack from scratch based on the latest versions of the DESI tools.


To compile the code, set ``$PLATFORM`` to one of the recipes in the
``platforms/`` and then run ``make install``;  *e.g.* at NERSC::

    make PLATFORM=nersc_desiconda install

The version of fiberassign can be updated with::

    make version

or::

    make TAG=1.2.3 version

to set the version when tagging.
