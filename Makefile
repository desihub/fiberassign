
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


all : 
	@cd src; $(MAKE)

install :
	@cd src; $(MAKE) install

clean :
	@cd src; $(MAKE) clean

install_writeMTL:
	@cd src_writeMTL; $(MAKE) install