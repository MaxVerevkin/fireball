OBJECTS = solver.o utils.o data.o structs.o
CXXFLAGS = -O3 -mavx2 -Wall
LDFLAGS =

all: $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o solver

paralel:
	$(MAKE) CXXFLAGS="-O3 -mavx2 -Wall -fopenmp -DOP_PARALEL" LDFLAGS="-fopenmp"

solver.o: solver.cpp utils.h structs.h data.h hyperparams.h simd.h
utils.o: utils.cpp utils.h structs.h simd.h
data.o: data.cpp data.h utils.h structs.h hyperparams.h simd.h
structs.o: structs.h

clean:
	rm -f solver *.o
