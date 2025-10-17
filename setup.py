#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
#
# Standard imports
#
import glob
import os
import sys
import shutil
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.command.clean import clean
from distutils.errors import CompileError
#
# Set keywords for the setup function. Most metadata is in setup.cfg.
# The keywords set below need some amount of automation, & should
# be left alone unless you are an expert.
#
# Treat everything in bin/ except *.rst as a script to be installed.
#
setup_keywords = dict()
if os.path.isdir('bin'):
    setup_keywords['scripts'] = [fname for fname in
                                 glob.glob(os.path.join('bin', '*'))
                                 if not os.path.basename(fname).endswith('.rst')]

# Add a custom clean command that removes in-tree files like the
# compiled extension.

class RealClean(clean):
    def run(self):
        super().run()
        clean_files = [
            "./build",
            "./dist",
            "py/fiberassign/_internal*",
            "py/fiberassign/__pycache__",
            "py/fiberassign/test/__pycache__",
            "./*.egg-info",
            "py/*.egg-info"
        ]
        for cf in clean_files:
            # Make paths absolute and relative to this path
            apaths = glob.glob(os.path.abspath(cf))
            for path in apaths:
                if os.path.isdir(path):
                    shutil.rmtree(path)
                elif os.path.isfile(path):
                    os.remove(path)
        return


# These classes allow us to build a compiled extension that uses pybind11.
# For more details, see:
#
#  https://github.com/pybind/python_example
#

# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689


def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    devnull = None
    oldstderr = None
    try:
        with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
            f.write('int main (int argc, char **argv) { return 0; }')
            try:
                devnull = open('/dev/null', 'w')
                oldstderr = os.dup(sys.stderr.fileno())
                os.dup2(devnull.fileno(), sys.stderr.fileno())
                compiler.compile([f.name], extra_postargs=[flagname])
            except CompileError:
                return False
            return True
    finally:
        if oldstderr is not None:
            os.dup2(oldstderr, sys.stderr.fileno())
        if devnull is not None:
            devnull.close()


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

    if sys.platform.lower() == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        linkopts = []
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' %
                        self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
            if has_flag(self.compiler, '-fopenmp'):
                opts.append('-fopenmp')
                linkopts.append('-fopenmp')
            if sys.platform.lower() == 'darwin':
                linkopts.append('-stdlib=libc++')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' %
                        self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args.extend(opts)
            ext.extra_link_args.extend(linkopts)

        # remove -Wstrict-prototypes flag
        if '-Wstrict-prototypes' in self.compiler.compiler_so:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")

        build_ext.build_extensions(self)

ext_modules = [
    Extension(
        'fiberassign._internal',
        [
            'src/utils.cpp',
            'src/hardware.cpp',
            'src/tiles.cpp',
            'src/targets.cpp',
            'src/assign.cpp',
            'src/_pyfiberassign.cpp'
        ],
        include_dirs=[
            'src',
        ],
        language='c++'
    ),
]

setup_keywords['ext_modules'] = ext_modules
setup_keywords['cmdclass'] = {'build_ext': BuildExt, 'clean': RealClean}
#
# Warning about old features.
#
VERSION_HELP = """
Note: Generating version strings is no longer done using 'python setup.py version'. Instead
you will need to run:

    desi_update_version [-t TAG] desiutil

which is part of the desiutil package. If you don't already have desiutil installed, you can install it with:

    pip install desiutil
"""

TEST_HELP = """
Note: running tests is no longer done using 'python setup.py test'. Instead
you will need to run:

    pytest

If you don't already have pytest installed, you can install it with:

    pip install pytest
"""

DOCS_HELP = """
Note: building the documentation is no longer done using
'python setup.py {0}'. Instead you will need to run:

    sphinx-build -W --keep-going -b html doc doc/_build/html

If you don't already have Sphinx installed, you can install it with:

    pip install Sphinx
"""

message = {'test': TEST_HELP,
           'version': VERSION_HELP,
           'build_docs': DOCS_HELP.format('build_docs'),
           'build_sphinx': DOCS_HELP.format('build_sphinx'), }

for m in message:
    if m in sys.argv:
        print(message[m])
        sys.exit(1)

#
# Copy the version string into the C++ code.
#
py_version_file = os.path.join('py', 'fiberassign', '_version.py')
cpp_version_file = os.path.join("src", "_version.h")
with open(py_version_file) as f:
    data = f.read()
    pkg_version = data.split('=')[1].strip().strip("'").strip('"')
with open(cpp_version_file, "w") as f:
    f.write('// Generated by setup.py -- DO NOT EDIT THIS\n')
    f.write('const static std::string package_version("{}");\n\n'
            .format(pkg_version))

#
# Run setup command.
#
setup(**setup_keywords)
