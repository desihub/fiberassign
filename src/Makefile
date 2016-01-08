# Makefile for building fiber assignment tools
# variables for compiler flags, etc are propagated
# from the top-level Makefile.  DO NOT use this
# Makefile directly (build from the top level).

HEADERS = \
collision.h \
feat.h \
global.h \
misc.h \
structs.h

OBJS_SURVEYSIM = \
collision.o \
feat.o \
global.o \
fa.o \
misc.o \
structs.o

OBJS_SKYSTD = \
collision.o \
feat.o \
global.o \
mtl.o \
misc.o \
structs.o

OBJS_TARG_SECRET = \
collision.o \
feat.o \
global.o \
targ.o \
misc.o \
structs.o

OBJS_FIBERASSIGN = \
collision.o \
feat.o \
global.o \
pipeline.o \
misc.o \
structs.o

EXECUTABLES = fiberassign fiberassign_surveysim split_targ_secret_sky_std split_targ_secret
COMPILE = $(CXXFLAGS) $(OMPFLAGS) -I. $(CFITSIO_CPPFLAGS)
LINK = $(OMPFLAGS) $(CFITSIO)

all : $(EXECUTABLES)

install : all
	cp $(EXECUTABLES) $(INSTALL)

split_targ_secret_sky_std : $(OBJS_SKYSTD)
	$(CXX) -o $@ $^ $(LINK)

split_targ_secret : $(OBJS_TARG_SECRET)
	$(CXX) -o $@ $^ $(LINK)

fiberassign : $(OBJS_FIBERASSIGN)
	$(CXX) -o $@ $^ $(LINK)

fiberassign_surveysim : $(OBJS_SURVEYSIM)
	$(CXX) -o $@ $^ $(LINK)

%.o : %.cpp $(HEADERS)
	$(CXX) $(COMPILE) -o $@ -c $<

clean :
	rm -f $(EXECUTABLES)
	rm -f *.o *~
	@cd $(INSTALL); rm -f $(EXECUTABLES)
