CC = g++
CFLAGS = -O3 -mavx2 -Wall

all: solver.o utils.o data.o structs.o
	$(CC) -o solver solver.o utils.o data.o structs.o

solver.o: solver.cpp utils.h structs.h data.h hyperparams.h simd.h
	$(CC) $(CFLAGS) -c solver.cpp

utils.o: utils.cpp utils.h structs.h simd.h
	$(CC) $(CFLAGS) -c utils.cpp

data.o: data.cpp data.h utils.h structs.h hyperparams.h data_values.h simd.h
	$(CC) $(CFLAGS) -c data.cpp

structs.o: structs.h
	$(CC) $(CFLAGS) -c structs.cpp

clean:
	rm -f solver solver.o utils.o data.o structs.o
