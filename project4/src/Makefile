CC = g++
MPICC = mpic++
CFLAGS = -O3 -Wall -std=c++0x -fopenmp
LFLAGS = 

OBJS = 
LIBS = 
TARGETS = MPIising

all: $(TARGETS)

run: MPIising
	mpirun -n 4 MPIising.x out.dat

MPIising: MPIising.cpp
	$(MPICC) -O3 -o MPIising.x MPIising.cpp

clean :
	rm -f *~ \#*# *.o *.x

