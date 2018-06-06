# Makefile for building fiber assignment tools
# variables for compiler flags, etc are propagated
# from the top-level Makefile.  DO NOT use this
# Makefile directly (build from the top level).

HEADERS = \
collision.h \
feat.h \
global.h \
misc.h \
structs.h \
version.h \
_version.h

OBJS_SURVEYSIM = \
surveysim.o \
collision.o \
feat.o \
global.o \
misc.o \
structs.o \
version.o

OBJS_FIBERASSIGN = \
fiberassign.o \
collision.o \
feat.o \
global.o \
misc.o \
structs.o \
version.o

EXECUTABLES = fiberassign_exec
COMPILE = $(CXXFLAGS) $(OMPFLAGS) -I. $(CFITSIO_CPPFLAGS)
# valgrind
#COMPILE = -g -O0 -m64 -fPIC -DNDEBUG -std=c++11 $(OMPFLAGS) -I. $(CFITSIO_CPPFLAGS)
LINK = $(OMPFLAGS) $(CFITSIO)


all : $(EXECUTABLES)

install : all
	/bin/cp -a $(EXECUTABLES) $(INSTALL_DIR)/bin

fiberassign_exec : $(OBJS_FIBERASSIGN)
	$(CXX) -o $@ $^ $(LINK)

fiberassign_surveysim : $(OBJS_SURVEYSIM)
	$(CXX) -o $@ $^ $(LINK)


%.o : %.cpp $(HEADERS)
	$(CXX) $(COMPILE) -o $@ -c $<

#
# GNU make pre-defines $(RM).  The - in front of $(RM) causes make to
# ignore any errors produced by $(RM).
#
clean :
	- $(RM) $(EXECUTABLES)
	- $(RM) *.o *~
	@ for e in $(EXECUTABLES); do $(RM) $(INSTALL_DIR)/bin/$$e ; done
#
# Enable 'make version' to update the version string.
# Do make TAG=0.1.2 version to set the tag explicitly.
# LASTTAG and COUNT are defined inside the recipe so they
# only get defined when make version is invoked.
#
version :
	- $(RM) _version.h
	$(eval LASTTAG := $$(shell git describe --tags --dirty --always | cut -d- -f1))
	$(eval COUNT := $$(shell git rev-list --count HEAD))
	@ if test -n "$(TAG)"; then v=$(TAG); else v=$(LASTTAG).dev$(COUNT); fi; \
		echo "#define VERSION_STRING \"$$v\"" > _version.h
