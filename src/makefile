CC=CC
CFLAGS= -fopenmp -funroll-loops -O2 -lm -lpthread -dynamic
DEPS = structs.h macros.h misc.h structs.cpp main.cpp macros.cpp misc.cpp global.cpp global.h
OBJ = macros.o misc.o structs.o global.o main.o

all: $(OBJ)
	$(CC) -o ../assign $(OBJ) $(CFLAGS)

misc.o: misc.h misc.cpp
	$(CC) -c misc.cpp $(CFLAGS)

macros.o: misc.h macros.h macros.cpp
	$(CC) -c macros.cpp $(CFLAGS)

structs.o: misc.h macros.h structs.h structs.cpp
	$(CC) -c structs.cpp $(CFLAGS)

global.o: misc.h macros.h structs.h global.h global.cpp
	$(CC) -c global.cpp $(CFLAGS)

main.o: misc.h macros.h structs.h global.h main.cpp 
	$(CC) -c main.cpp $(CFLAGS)

clean:
	rm $(OBJ)
