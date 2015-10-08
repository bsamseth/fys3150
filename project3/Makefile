CC = c++
CFLAGS = -O2 -Wall -std=c++11 -fopenmp
LFLAGS = 

OBJS = 
LIBS = -lunittest++ 
TARGETS = gauss_legendre

all: $(TARGETS)


gauss_legendre: gauss_legendre.o
	$(CC) $(CFLAGS) gauss_legendre.o -o gauss_legendre.x $(LFLAGS)

gauss_legendre.o:  gauss_legendre.cpp
	$(CC) $(CFLAGS) -c gauss_legendre.cpp

clean :
	rm -f *~ \#*# *.o *.x

