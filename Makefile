CC = g++
CFLAGS = -O3 -mavx2 -Wall -fopenmp

all: solver.o utils.o data.o
	$(CC) -o solver solver.o utils.o data.o -fopenmp

solver.o: solver.cpp utils.h structs.h data.h simd.h
	$(CC) $(CFLAGS) -c solver.cpp

utils.o: utils.cpp utils.h structs.h simd.h
	$(CC) $(CFLAGS) -c utils.cpp

data.o: data.cpp data.h utils.h structs.h simd.h
	$(CC) $(CFLAGS) -c data.cpp

clean:
	rm -f solver solver.o utils.o data.o
