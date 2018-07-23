#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function
#
# Standard imports
#
import glob
import os
import sys
#
# setuptools' sdist command ignores MANIFEST.in
#
from distutils.command.sdist import sdist as DistutilsSdist
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
#
# DESI support code.
#
from desiutil.setup import DesiTest, DesiVersion, get_version
#
# Begin setup
#
setup_keywords = dict()
#
# THESE SETTINGS NEED TO BE CHANGED FOR EVERY PRODUCT.
#
setup_keywords['name'] = 'fiberassign'
setup_keywords['description'] = 'DESI Fiber Assignment Tools'
setup_keywords['author'] = 'DESI Collaboration'
setup_keywords['author_email'] = 'desi-data@desi.lbl.gov'
setup_keywords['license'] = 'BSD'
setup_keywords['url'] = 'https://github.com/desihub/fiberassign'
#
# END OF SETTINGS THAT NEED TO BE CHANGED.
#
setup_keywords['version'] = get_version(setup_keywords['name'])
#
# Use README.rst as long_description.
#
setup_keywords['long_description'] = ''
if os.path.exists('README.rst'):
    with open('README.rst') as readme:
        setup_keywords['long_description'] = readme.read()
#
# Set other keywords for the setup function.  These are automated, & should
# be left alone unless you are an expert.
#
# Treat everything in bin/ except *.rst as a script to be installed.
#
if os.path.isdir('bin'):
    setup_keywords['scripts'] = [fname for fname in glob.glob(os.path.join('bin', '*'))
        if not os.path.basename(fname).endswith('.rst')]
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['requires'] = ['Python (>2.7.0)']
# setup_keywords['install_requires'] = ['Python (>2.7.0)']
setup_keywords['zip_safe'] = False
setup_keywords['use_2to3'] = False
setup_keywords['packages'] = find_packages('py')
setup_keywords['package_dir'] = {'':'py'}
setup_keywords['cmdclass'] = {'version': DesiVersion, 'test': DesiTest, 'sdist': DistutilsSdist}
setup_keywords['test_suite']='{name}.test.{name}_test_suite.{name}_test_suite'.format(**setup_keywords)
#
# Autogenerate command-line scripts.
#
# setup_keywords['entry_points'] = {'console_scripts':['desiInstall = desiutil.install.main:main']}
#
# Add internal data directories.
#
#setup_keywords['package_data'] = {'fiberassign': ['data/*',]}
#

def find_cfitsio():
    """Helper function to find cfitsio include and library paths.

    If the user has specified the CFITSIO_DIR environment variable, then it
    uses that blindly.  Otherwise it tries to use pkg-config to get the needed
    flags.  If that fails, it just assumes that cfitsio is in the default CPATH
    and LD_LIBRARY_PATH / LIBRARY_PATH search locations.

    """
    inc = ""
    ldf = ""
    if "CFITSIO_DIR" in os.environ:
        inc = os.path.join(os.environ["CFITSIO_DIR"], "include")
        ldf = os.path.join(os.environ["CFITSIO_DIR"], "lib")
    else:
        import subprocess as sp
        import re
        try:
            sp.check_call(["pkg-config", "--version"])
            # Ok, we have pkg-config...
            inc = sp.check_output(
                ["pkg-config", "--cflags-only-I", "cfitsio"],
                universal_newlines=True).strip()
            ldf = sp.check_output(
                ["pkg-config", "--libs-only-L", "cfitsio"],
                universal_newlines=True).strip()
            inc = inc.replace("-I", "")
            ldf = ldf.replace("-L", "")
        except OSError as e:
            # Leave the variables blank and assume they are in the default
            # search path
            pass

    return inc, ldf


# These classes allow us to build a compiled extension that uses pybind11.
# For more details, see:
#
#  https://github.com/pybind/python_example
#

cfitsio_inc, cfitsio_ldf = find_cfitsio()
print("Using CFITSIO include directory: {}".format(cfitsio_inc))
print("Using CFITSIO library directory: {}".format(cfitsio_ldf))

ext_modules = [
    Extension(
        'fiberassign._internal',
        [
            'src/structs.cpp',
            'src/misc.cpp',
            'src/feat.cpp',
            'src/collision.cpp',
            'src/global.cpp',
            'src/fiberassign.cpp',
            'src/_pyfiberassign.cpp'
        ],
        include_dirs=[
            'src',
            cfitsio_inc,
        ],
        library_dirs=[cfitsio_ldf],
        libraries=['cfitsio'],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'],
        language='c++'
    ),
]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.

    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args.extend(opts)
        build_ext.build_extensions(self)


setup_keywords['ext_modules'] = ext_modules
setup_keywords['cmdclass']['build_ext'] = BuildExt

#
# Run setup command.
#
setup(**setup_keywords)
