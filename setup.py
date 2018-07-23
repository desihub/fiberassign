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
from setuptools import setup, find_packages
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

# These classes allow us to build a compiled extension that uses pybind11.
# For more details, see:
#
#  https://github.com/pybind/python_example
#

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


ext_modules = [
    Extension(
        '_fiberassign',
        [
            'src/structs.cpp',
            'src/misc.cpp',
            'src/feat.cpp',
            'src/collision.cpp',
            'src/global.cpp',
            'src/fiberassign.cpp',
        ],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True),
            'src',
        ],
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
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)


setup_keywords['ext_modules'] = ext_modules
setup_keywords['cmdclass']['build_ext'] = BuildExt

#
# Run setup command.
#
setup(**setup_keywords)
