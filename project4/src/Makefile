CC = g++
MPICC = mpic++
CFLAGS = -O3 -Wall -std=c++0x -fopenmp
LFLAGS = 

OBJS = functions.o
LIBS = 
TARGETS = MPIising MPIisingImproved

all: $(TARGETS)

MPIising: $(OBJS) MPIising.cpp 
	$(MPICC) -O3 $(OBJS) -o MPIising.x MPIising.cpp

MPIisingImproved: $(OBJS) MPIisingImproved.cpp 
	$(MPICC) -O3 $(OBJS) -o MPIisingImproved.x MPIisingImproved.cpp

functions.o: functions.cpp functions.h
	$(CC) $(CFLAGS) -c functions.cpp

clean :
	rm -f *~ \#*# *.o *.x

