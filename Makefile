OBJECTS = solver.o utils.o data.o structs.o
CXXFLAGS = -O3 -mavx2 -Wall
LDFLAGS =

all: $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o solver

profile:
	$(MAKE) CXXFLAGS="-O3 -mavx2 -Wall -pg" LDFLAGS="-pg"

debug:
	$(MAKE) CXXFLAGS="-O3 -mavx2 -Wall -g" LDFLAGS="-g"

solver.o: solver.cpp utils.h structs.h data.h hyperparams.h simd.h
utils.o: utils.cpp utils.h structs.h simd.h
data.o: data.cpp data.h utils.h structs.h hyperparams.h simd.h
structs.o: structs.cpp structs.h utils.h

clean:
	rm -f solver *.o
