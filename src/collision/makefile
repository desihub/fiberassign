CC=CC
CFLAGS= -fopenmp -funroll-loops -O2 -lm -lpthread -dynamic
DEPS = ../macros.cpp ../misc.cpp ../macros.h ../misc.h main.cpp collide.cpp collide.h
OBJ = macros.o misc.o collide.o main.o

all: $(OBJ)
	$(CC) -o collision $(OBJ) $(CFLAGS)

misc.o: ../misc.h ../misc.cpp
	$(CC) -c ../misc.cpp $(CFLAGS)

macros.o: ../macros.h ../macros.cpp
	$(CC) -c ../macros.cpp $(CFLAGS)

collide.o: ../misc.h ../macros.h collide.h collide.cpp
	$(CC) -c collide.cpp $(CFLAGS)

main.o: ../misc.h ../macros.h collide.h main.cpp 
	$(CC) -c main.cpp $(CFLAGS)

clean:
	rm $(OBJ)
