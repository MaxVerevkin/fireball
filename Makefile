OBJECTS = solver.o utils.o data.o structs.o
CXXFLAGS = -O3 -mavx2 -Wall
LDFLAGS =

all: $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o solver

paralel:
	$(MAKE) CXXFLAGS="-O3 -mavx2 -Wall -fopenmp -DOP_PARALEL" LDFLAGS="-fopenmp"

profile:
	$(MAKE) CXXFLAGS="-O3 -mavx2 -Wall -pg" LDFLAGS="-pg"

debug:
	$(MAKE) CXXFLAGS="-O3 -mavx2 -Wall -g" LDFLAGS="-g"

solver.o: solver.cpp utils.h structs.h data.h hyperparams.h simd.h
utils.o: utils.cpp utils.h structs.h simd.h
data.o: data.cpp data.h utils.h structs.h hyperparams.h simd.h
structs.o: structs.h

clean:
	rm -f solver *.o
