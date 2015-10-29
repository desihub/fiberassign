
# check if we are using platform-specific options,
# otherwise use the 

ifndef PLATFORM
  PLATFORM := generic
endif

include platforms/$(PLATFORM)

# absolute path to top source directory

TOPDIR := $(shell pwd)
export TOPDIR

# where to place built executables

ifndef INSTALL
  INSTALL := $(TOPDIR)
endif

export INSTALL


all_mtl : 
	@cd src; $(MAKE) all_mtl

install_mtl :
	@cd src; $(MAKE) install_mtl

clean_mtl :
	@cd src; $(MAKE) clean_mtl

all_fa : 
	@cd src; $(MAKE) all_fa

install_fa :
	@cd src; $(MAKE) install_fa

clean_fa :
	@cd src; $(MAKE) clean_fa

all_mini : 
	@cd src; $(MAKE) all_mini

install_mini :
	@cd src; $(MAKE) install_mini

clean_mini :
	@cd src; $(MAKE) clean_mini
