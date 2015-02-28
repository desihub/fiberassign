CC=CC
CFLAGS= -fopenmp -funroll-loops -O2 -lm -lpthread -dynamic
DEPS = structs.h macros.h misc.h structs.cpp main.cpp macros.cpp misc.cpp
OBJ = macros.o misc.o structs.o main.o

all: $(OBJ)
	$(CC) -o ../assign $(OBJ) $(CFLAGS)

macros.o: macros.h macros.cpp
	$(CC) -c macros.cpp $(CFLAGS)

misc.o: macros.h misc.cpp misc.h
	$(CC) -c misc.cpp $(CFLAGS)

structs.o: macros.h misc.h structs.cpp structs.h
	$(CC) -c structs.cpp $(CFLAGS)

main.o: macros.h misc.h structs.h main.cpp 
	$(CC) -c main.cpp $(CFLAGS)

clean:
	rm $(OBJ)
