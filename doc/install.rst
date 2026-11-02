.. _install:

Installation
============

Requirements
------------

- Python >= 3.10 (see ``requires-python`` in ``pyproject.toml``).
- A C++ compiler with C++14 support.  ``setup.py`` probes for ``-std=c++14``
  and falls back to ``-std=c++11``.
- The usual scientific python stack: numpy, scipy, astropy, matplotlib,
  fitsio, pyyaml, healpy.
- The DESI packages desiutil, desimodel, desitarget and desimeter.  Note that
  some of these are imported directly by the fiberassign code even though they
  are not all listed in ``pyproject.toml`` dependencies.


DESI dependencies
-----------------

Install the DESI packages from github (the default branch is ``main``)::

    %> pip install git+https://github.com/desihub/desiutil.git@main
    %> pip install git+https://github.com/desihub/desimodel.git@main
    %> pip install git+https://github.com/desihub/desitarget.git@main
    %> pip install git+https://github.com/desihub/desimeter.git@main

fiberassign needs the desimodel data, pointed to by the ``$DESIMODEL``
environment variable.  For real survey work, export the full data::

    %> export DESIMODEL=${HOME}/software/desimodel
    %> svn export https://desi.lbl.gov/svn/code/desimodel/trunk/data ${DESIMODEL}/data

For running the test suite, the small test branch used by CI is sufficient::

    %> svn export https://desi.lbl.gov/svn/code/desimodel/branches/test-0.19/data ${DESIMODEL}/data


Installing fiberassign
----------------------

User
~~~~

Install a fixed version directly from github::

    %> pip install git+https://github.com/desihub/fiberassign.git@main

Developer / contributor
~~~~~~~~~~~~~~~~~~~~~~~

Clone the repository and compile the C++ extension in place::

    %> git clone https://github.com/desihub/fiberassign.git
    %> cd fiberassign
    %> CC=g++ CXX=g++ CFLAGS="-fPIC" LDSHARED="g++ -shared" python setup.py build_ext --inplace

Then use it from the source tree by adding the package and scripts to your
search paths::

    %> export PYTHONPATH=${PWD}/py
    %> export PATH=${PWD}/bin:$PATH

Alternatively install it into your environment with ``pip install .`` (or
``pip install -e .`` for an editable install).  The ``CC`` / ``CXX`` /
``CFLAGS`` / ``LDSHARED`` variables above are needed on some platforms and are
what CI uses; on others the default compiler invocation just works.

Note that ``python setup.py`` only supports the ``build_ext`` and ``clean``
commands.  Running ``test``, ``version``, ``build_docs`` or ``build_sphinx``
prints a message explaining the modern replacement and exits.


Verifying the installation
--------------------------

With ``$DESIMODEL`` set, run the test suite::

    %> pytest

and check that the package imports::

    %> python -c "import fiberassign; print(fiberassign.__version__)"


Building the documentation
--------------------------

The documentation is built with Sphinx::

    %> pip install -r doc/rtd_requirements.txt
    %> cd doc && make html

or equivalently ``sphinx-build -W --keep-going -b html doc doc/_build/html``.
The build imports fiberassign, so the C++ extension must be compiled first
(see above).


Use at NERSC
------------

At NERSC, you can use the pre-installed version of all DESI tools.  There is a
nightly snapshot of the main branch accessible with::

    %> source /global/common/software/desi/desi_environment.sh main

All tools (including fiberassign) should then be in your executable and python
search paths.  After loading this, you can separately check out and install
your own fiberassign version if desired.


Troubleshooting
---------------

- ``ModuleNotFoundError`` for a DESI package or for ``fiberassign._internal``:
  check that the DESI dependencies are installed, and that the C++ extension
  was built (``py/fiberassign/_internal*.so`` should exist after ``build_ext
  --inplace``) or that fiberassign was pip-installed.
- Tests failing with desimodel errors: make sure ``$DESIMODEL`` is set and
  points at a directory containing ``data/``.
- If a local build gets messed up, ``python setup.py clean`` removes the build
  artifacts and the compiled extension.
