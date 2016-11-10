
# check if we are using platform-specific options,
# otherwise use the 

ifndef PLATFORM
  PLATFORM := generic
endif

include platforms/$(PLATFORM)

# absolute path to top source directory

export TOPDIR := $(shell pwd)

# where to place built executables
ifdef INSTALL_DIR
export INSTALL := $(INSTALL_DIR)/bin
endif
ifndef INSTALL
export INSTALL := $(TOPDIR)/bin
endif

all : 
	$(MAKE) -C src all

install :
	$(MAKE) -C src install

clean :
	$(MAKE) -C src clean
