CC = c++
CFLAGS = -O2 -Wall -std=c++11 -fopenmp
LFLAGS = 

OBJS = 
LIBS = -lunittest++ 
TARGETS = gauss_legendre gauss_laguerre uniformmonteCarlo

all: $(TARGETS)

uniformmonteCarlo: uniformmonteCarlo.o
	$(CC) $(CFLAGS) uniformmonteCarlo.o -o uniformmonteCarlo.x $(LFLAGS)

uniformmonteCarlo.o:  uniformmonteCarlo.cpp
	$(CC) $(CFLAGS) -c uniformmonteCarlo.cpp


gauss_laguerre: gauss_laguerre.o
	$(CC) $(CFLAGS) gauss_laguerre.o -o gauss_laguerre.x $(LFLAGS)

gauss_laguerre.o:  gauss_laguerre.cpp
	$(CC) $(CFLAGS) -c gauss_laguerre.cpp


gauss_legendre: gauss_legendre.o
	$(CC) $(CFLAGS) gauss_legendre.o -o gauss_legendre.x $(LFLAGS)

gauss_legendre.o:  gauss_legendre.cpp
	$(CC) $(CFLAGS) -c gauss_legendre.cpp

clean :
	rm -f *~ \#*# *.o *.x

