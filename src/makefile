CC=CC
CFLAGS= -fopenmp -funroll-loops -O2 -lm -lpthread -dynamic
DEPS = structs.h feat.h misc.h collision.h collision.cpp structs.cpp main.cpp feat.cpp misc.cpp global.cpp global.h
OBJ = misc.o collision.o feat.o structs.o global.o main.o

all: $(OBJ)
	$(CC) -o ../assign $(OBJ) $(CFLAGS)


misc.o: misc.h misc.cpp
	$(CC) -c misc.cpp $(CFLAGS)

collision.o: collision.h collision.cpp
	$(CC) -c collision.cpp $(CFLAGS)

feat.o: feat.h feat.cpp
	$(CC) -c feat.cpp $(CFLAGS)

structs.o: misc.h feat.h collision.h structs.h structs.cpp
	$(CC) -c structs.cpp $(CFLAGS)

global.o: misc.h feat.h collision.h structs.h global.h global.cpp
	$(CC) -c global.cpp $(CFLAGS)

main.o: misc.h collision.h feat.h structs.h global.h main.cpp 
	$(CC) -c main.cpp $(CFLAGS)

clean:
	rm $(OBJ)
