CC = g++
CFLAGS = -O3 -Wall -std=c++0x -fopenmp
LFLAGS = 

OBJS = functions.o
LIBS = -lunittest++ 
TARGETS = gauss_legendre gauss_laguerre uniformmonteCarlo smartMonteCarlo

all: $(TARGETS)

smartMonteCarlo: smartMonteCarlo.o
	$(CC) $(CFLAGS) $(OBJS) smartMonteCarlo.o -o smartMonteCarlo.x $(LFLAGS)

smartMonteCarlo.o: $(OBJS) smartMonteCarlo.cpp
	$(CC) $(CFLAGS) -c smartMonteCarlo.cpp

uniformmonteCarlo: uniformmonteCarlo.o
	$(CC) $(CFLAGS) $(OBJS) uniformmonteCarlo.o -o uniformmonteCarlo.x $(LFLAGS)

uniformmonteCarlo.o: $(OBJS) uniformmonteCarlo.cpp
	$(CC) $(CFLAGS) -c uniformmonteCarlo.cpp


gauss_laguerre: gauss_laguerre.o
	$(CC) $(CFLAGS) $(OBJS) gauss_laguerre.o -o gauss_laguerre.x $(LFLAGS)

gauss_laguerre.o: $(OBJS) gauss_laguerre.cpp
	$(CC) $(CFLAGS) -c gauss_laguerre.cpp


gauss_legendre: gauss_legendre.o
	$(CC) $(CFLAGS) $(OBJS) gauss_legendre.o -o gauss_legendre.x $(LFLAGS)

gauss_legendre.o: $(OBJS) gauss_legendre.cpp
	$(CC) $(CFLAGS) -c gauss_legendre.cpp

functions.o: functions.cpp functions.h
	$(CC) $(CFLAGS) -c functions.cpp

clean :
	rm -f *~ \#*# *.o *.x

