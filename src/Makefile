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

OBJS_FA = \
collision.o \
feat.o \
global.o \
fa.o \
misc.o \
structs.o

OBJS_MTL = \
collision.o \
feat.o \
global.o \
mtl.o \
misc.o \
structs.o

COMPILE = $(CXXFLAGS) $(OMPFLAGS) -I. $(CFITSIO_CPPFLAGS)
LINK = $(OMPFLAGS) $(CFITSIO)


all_mtl : assign_mtl

install_mtl : all_mtl
	cp assign_mtl $(INSTALL)

all_fa : assign_fa

install_fa : all_fa
	cp assign_fa $(INSTALL)

assign_mtl : $(OBJS_MTL)
	$(CXX) -o $@ $^ $(LINK)

assign_fa : $(OBJS_FA)
	$(CXX) -o $@ $^ $(LINK)

%.o : %.cpp $(HEADERS)
	$(CXX) $(COMPILE) -o $@ -c $<

clean_mtl :
	rm -f assign_mtl  *.o *~

clean_fa :
	rm -f assign_fa  *.o *~
