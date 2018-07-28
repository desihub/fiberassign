.. _install:


Installation
===============

The DESI fiber assignment code is a hybrid code written in C++ and Python.
Installation requires support for at least the C++11 standard.

External Dependencies
------------------------

In addition to the usual numpy / scipy software stack, the "fitsio" package
is also required.

DESI Affiliated Dependencies
---------------------------------

Fiberassign requires some additional DESI packages in order to function.
Currently these include desimodel and desitarget.

Installing fiberassign
-----------------------------

To install fiberassign to some location, run::

    %> python setup.py clean
    %> python setup.py install --prefix <somewhere>

To compile fiberassign and use it from the source tree::

    %> python setup.py clean; python3 setup.py build

In either case, you will need to modify your PATH and PYTHONPATH variables to
point to the copy you wish to use.  Alternatively, if you are installing to a
virtualenv or a conda environment, then you can omit the :code:`--prefix` option and
just install it to the default location.
