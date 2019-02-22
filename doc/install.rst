.. _install:


Use at NERSC
==============================

At NERSC, you can use the pre-installed version of all DESI tools.  There is a
nightly snapshot of the master branch accessible with (assuming you use bash)::

    %> . /project/projectdirs/desi/software/desi_environment.sh master

All tools (including fiberassign) should then be in your executable and python
search paths.  After loading this, you can separately checkout and install your
own fiberassign version if desired.


Installation
===============

The DESI fiber assignment code is a hybrid code written in C++ and Python.
Installation requires support for at least the C++11 standard.  Compilation of
the internal fiberassign extension should "just work" using whatever compiler
is the default for the Python you install.  The compilation has been tested on
Linux with gcc 4.9 through 8.0 and on OS X with both the system clang and the
gcc wrapper.


External Dependencies
------------------------

In addition to the usual numpy / scipy software stack, the "fitsio" package is
also required.  There are multiple ways of installing a working python3 stack
on both Linux and OS X.  The solution you choose likely depends on what other
things you are using Python for- not just fiberassign.  In these examples,
we'll be creating a python stack in ${HOME}/software/desi, however if you
already have a python stack for use with DESI tools, just skip this section.

Use Anaconda...
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The officially supported DESI python environment is Anaconda.  After installing
the main Anaconda distribution, the "conda" command should be available.
Create a new conda environment::

  %> conda create --copy -m -p ${HOME}/software/desi
  %> conda activate ~/software/desi
  %> conda install numpy scipy astropy matplotlib
  %> pip install --no-binary :all: fitsio

... Or Use Virtualenv and Pip
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install Python3 and the virtualenv package with your OS package manager (Linux)
or using macports / homebrew (OS X).  Next create a virtualenv and activate
it::

  %> virtualenv -p python3 ${HOME}/software/desi
  %> source ${HOME}/software/desi/bin/activate

Now use pip to install the dependencies we need::

  %> pip install numpy scipy astropy matplotlib
  %> pip install --no-binary :all: fitsio


DESI Affiliated Dependencies
---------------------------------

Activate / load your python stack from the previous section.  Since you created
a conda environment or virtualenv directory specifically for DESI tools (or
perhaps even just for fiber assignment), you can always delete that directory
and make a new one as needed.

Fiberassign requires some additional DESI packages in order to function.
Currently these include desiutil, desimodel and desitarget.  The best way to
install and use these DESI packages depends on your use case.  Here we outline
some different scenarios.

DESI "User"
~~~~~~~~~~~~~~~~~~~~~

If you are using these DESI packages but have no intention of hacking on that
code, then you may not want or need a persistent git checkout of any code.  In
that case, you can pip install everything directly from github.::

    %> pip install git+https://github.com/desihub/desiutil.git@master#egg=desiutil
    %> pip install git+https://github.com/desihub/desimodel.git@master#egg=desimodel
    %> pip install git+https://github.com/desihub/desitarget.git@master#egg=desitarget

You also need to install the desimodel data::

    %> export DESIMODEL=${HOME}/software/desi
    %> svn export https://desi.lbl.gov/svn/code/desimodel/trunk/data \
     ${DESIMODEL}/data

You will need to ensure that $DESIMODEL is set in your environment before
running fiberassign.  I suggest making a shell function / alias that activates
your python stack and sets this environment variable.

DESI "Developer"
~~~~~~~~~~~~~~~~~~~~~

If you are routinely working on DESI code you likely already have git checkouts
of these core packages.  Make sure your git checkouts are up to date and
install each of the packages into your conda env / virtualenv with::

    %> python setup.py clean
    %> python setup.py install

And of course we still need a copy of the desimodel data.  If you already have
a local copy of that, then just make sure that the DESIMODEL environment
variable is set appropriately or export a copy into your conda env / virtualenv
as in the last section.


Installing fiberassign
-----------------------------

Here again we have a choice of how to install fiberassign depending on our use
case.

Fiberassign User
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are using fiberassign, but are not going to be working on the code or
submitting pull requests, etc, then you can just pip install directly from
github::

    %> pip install git+https://github.com/desihub/fiberassign.git@master#egg=fiberassign


Fiberassign Developer / Contributor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If there is any chance you might want to work on the fiberassign codebase, even
just to contribute small bug reports, it is useful to have a persistent git
checkout of the code.  Git clone the repository.  You can then install it into
your conda env / virtualenv with::

    %> python setup.py clean
    %> python setup.py install

Alternatively, you can install it to some other location::

    %> python setup.py clean
    %> python setup.py install --prefix <somewhere>

Or compile fiberassign and use it from the source tree::

    %> python setup.py clean; python setup.py build_ext --inplace

In the last two cases, you will need to modify your PATH and PYTHONPATH
variables to point to the copy you wish to use.

Something Went Wrong!
---------------------------

If something gets messed up, it's ok.  You are using a separate conda / virtualenv environment so you can just do::

    %> rm -rf ${HOME}/software/desi

and start over from the beginning.  If you encounter an installation problem, please open a github ticket and provide the following details:

- OS and version.

- Where you got your python (Anaconda, OS packages, homebrew, etc).

- Python version (:code:`python --version`).

- Any compiler error / output.
