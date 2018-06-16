#
# Template Makefile for use with desiInstall.  You can assume that
# desiInstall will set these environment variables:
#
# WORKING_DIR   : The directory containing the svn export
# INSTALL_DIR   : The directory the installed product will live in.
# (PRODUCT)     : Where (PRODUCT) is replaced with the name of the
#                 product in upper case, e.g. DESITEMPLATE.  This should
#                 be the same as WORKING_DIR for typical installs.
#
# Use this shell to interpret shell commands, & pass its value to sub-make
#
export SHELL = /bin/sh
#
# This is like doing 'make -w' on the command line.  This tells make to
# print the directory it is in.
#
MAKEFLAGS = w
#
# This is a list of subdirectories that make should descend into.  Makefiles
# in these subdirectories should also understand 'make all' & 'make clean'.
# This list can be empty, but should still be defined.
#
SUBDIRS = src
#
# Set INSTALL_DIR
#
ifndef INSTALL_DIR
  export INSTALL_DIR := $(shell pwd)
endif
#
# Check if we are using platform-specific options, otherwise use the generic
# configuration.
#
ifndef PLATFORM

ifdef NERSC_HOST

# At NERSC, we want to use the harpconfig tool from desiconda if
# available.  Otherwise, we will error.
ifeq "x$(DESICONDA_EXTRA)" "x"
define errmsg

You are building at NERSC, but do not have desiconda loaded.
Either load that or manually install CFITSIO with gcc and use
PLATFORM=nersc_gcc.

endef
$(error $(errmsg))
endif
PLATFORM := nersc_desiconda

else

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
PLATFORM := osx
else

define defaultmsg

You have not selected a build configuration.  Generic Linux
settings will be used.  To customize your build, make a new
file in the platforms directory and then do:

%> PLATFORM=<name> make <command>

When running make.

endef
$(warning $(defaultmsg))
PLATFORM := generic

endif

endif

endif

include platforms/$(PLATFORM)

# Make sure that these targets are 'actions' not files.
.PHONY : all clean install version

# This should compile all code prior to it being installed.

all :
	@ for f in $(SUBDIRS); do $(MAKE) -C $$f all ; done

install : all
	- /bin/mkdir -p $(INSTALL_DIR)
	@ for d in bin py script test; do if [ ! -d $(INSTALL_DIR)/$$d ]; then /bin/cp -a $$d $(INSTALL_DIR) ; fi ; done
	@ for f in $(SUBDIRS); do $(MAKE) -C $$f install ; done

clean :
	@ for f in $(SUBDIRS); do $(MAKE) -C $$f clean ; done

# Enable 'make version' to update the version string.
# Do make TAG=0.1.2 version to set the tag explicitly.
version :
	$(MAKE) -C src version
