.. _install:


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

In addition to the usual numpy / scipy software stack, the "fitsio" package
is also required.  There are multiple ways of installing a working python3 stack on both Linux and OS X.  The solution you choose likely depends on what other things you are using Python for- not just fiberassign.  In these examples, we'll be creating a python stack in ${HOME}/software/desi, however
if you already have a python stack for use with DESI tools, just skip this section.

Anaconda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The officially supported DESI python environment is Anaconda.  After installing the main Anaconda distribution, the "conda" command should be available.  Create a new conda environment::

  %> conda create --copy -m -p ${HOME}/software/desi
  %> conda activate ~/software/desi
  %> conda install numpy scipy astropy matplotlib
  %> pip install --no-binary :all: fitsio

Virtualenv and Pip
~~~~~~~~~~~~~~~~~~~~~~~

Install Python3 and the virtualenv package with your OS package manager (Linux) or using macports / homebrew (OS X).  Next create a virtualenv and activate it::

  %> virtualenv -p python3 ${HOME}/software/desi
  %> source ${HOME}/software/desi/bin/activate

Now use pip to install the dependencies we need::

  %> pip install numpy scipy astropy
  %> pip install --no-binary :all: fitsio




DESI Affiliated Dependencies
---------------------------------

Activate / load your python stack from the previous section.  Fiberassign
requires some additional DESI packages in order to function.  Currently these
include desimodel and desitarget.  If you already have these packages installed
for other purposes (for example you are developing these packages), then you
can skip this section.  Here we will be installing the required DESI packages
into our python stack using pip.

  %> pip install git+https://github.com/desihub/desiutil.git@master#egg=desiutil
  %> pip install git+https://github.com/desihub/desimodel.git@master#egg=desimodel
  %> pip install git+https://github.com/desihub/desitarget.git@master#egg=desitarget

We also need to install the desimodel data:

  %> export DESIMODEL=${HOME}/software/desi
  %> svn export https://desi.lbl.gov/svn/code/desimodel/trunk/data \
     ${DESIMODEL}/data

You will need to ensure that $DESIMODEL is set in your environment before
running fiberassign.  I suggest making a shell function / alias that activates
your python stack and sets this environment variable.


Installing fiberassign
-----------------------------

To install fiberassign to some arbitrary location, run::

    %> python setup.py clean
    %> python setup.py install --prefix <somewhere>

To compile fiberassign and use it from the source tree::

    %> python setup.py clean; python setup.py build_ext --inplace

In either case, you will need to modify your PATH and PYTHONPATH variables to
point to the copy you wish to use.  Alternatively, if you are installing to a
virtualenv or a conda environment, then you can omit the :code:`--prefix`
option and just install it to the default location.
