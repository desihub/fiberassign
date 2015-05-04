CC=CC
CFLAGS= -fopenmp -funroll-loops -O2 -lm -lpthread -dynamic
DEPS = structs.h macros.h misc.h collision.h collision.cpp structs.cpp main.cpp macros.cpp misc.cpp global.cpp global.h collision/collide.cpp collision/collide.h
OBJ = macros.o misc.o collision.o structs.o global.o main.o
OBJ2 = macros.o misc.o collision/collide.o collision/maincol.o

all: $(OBJ)
	$(CC) -o ../assign $(OBJ) $(CFLAGS)

allcol: $(OBJ2)
	$(CC) -o ../collide $(OBJ2) $(CFLAGS)




misc.o: misc.h misc.cpp
	$(CC) -c misc.cpp $(CFLAGS)

macros.o: macros.h macros.cpp
	$(CC) -c macros.cpp $(CFLAGS)

collision.o: collision.h collision.cpp
	$(CC) -c collision.cpp $(CFLAGS)

structs.o: misc.h macros.h collision.h structs.h structs.cpp
	$(CC) -c structs.cpp $(CFLAGS)

global.o: misc.h macros.h collision.h structs.h global.h global.cpp
	$(CC) -c global.cpp $(CFLAGS)

main.o: misc.h macros.h collision.h structs.h global.h main.cpp 
	$(CC) -c main.cpp $(CFLAGS)




collide.o: misc.h macros.h collision/collide.h collision/collide.cpp
	$(CC) -c collision/collide.cpp $(CFLAGS) collision/collide.o

maincol.o: misc.h macros.h collision/collide.h collision/maincol.cpp 
	$(CC) -c collision/maincol.cpp $(CFLAGS) collision/maincol.o

clean:
	rm $(OBJ)
	rm $(OBJ2)
