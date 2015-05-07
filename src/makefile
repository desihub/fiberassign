CC=CC
CFLAGS= -fopenmp -funroll-loops -O2 -lm -lpthread -dynamic
DEPS = structs.h feat.h misc.h collision.h collision.cpp structs.cpp main.cpp feat.cpp misc.cpp global.cpp global.h collision/collide.cpp collision/collide.h
OBJ = misc.o feat.o collision.o structs.o global.o main.o
OBJ2 = misc.o feat.o collision/collide.o collision/maincol.o

all: $(OBJ)
	$(CC) -o ../assign $(OBJ) $(CFLAGS)

allcol: $(OBJ2)
	$(CC) -o ../collide $(OBJ2) $(CFLAGS)




misc.o: misc.h misc.cpp
	$(CC) -c misc.cpp $(CFLAGS)

feat.o: feat.h feat.cpp
	$(CC) -c feat.cpp $(CFLAGS)

collision.o: collision.h collision.cpp
	$(CC) -c collision.cpp $(CFLAGS)

structs.o: misc.h feat.h collision.h structs.h structs.cpp
	$(CC) -c structs.cpp $(CFLAGS)

global.o: misc.h feat.h collision.h structs.h global.h global.cpp
	$(CC) -c global.cpp $(CFLAGS)

main.o: misc.h feat.h collision.h structs.h global.h main.cpp 
	$(CC) -c main.cpp $(CFLAGS)




collide.o: misc.h feat.h collision/collide.h collision/collide.cpp
	$(CC) -c collision/collide.cpp $(CFLAGS) collision/collide.o

maincol.o: misc.h feat.h collision/collide.h collision/maincol.cpp 
	$(CC) -c collision/maincol.cpp $(CFLAGS) collision/maincol.o

clean:
	rm $(OBJ)
	rm $(OBJ2)
