CC = g++
CFLAGS = -O3 -Wall -std=c++0x
LFLAGS = 

OBJS = Body.o Universe.o
LIBS = ../../lib/lib.o
TARGETS = main

all: $(TARGETS)

main: $(OBJS) main.cpp
	$(CC) $(CFLAGS) -Wno-unused-variable $(OBJS) $(LIBS) -o main.x main.cpp

Universe.o: Body.o Universe.cpp Universe.h
	$(CC) $(CFLAGS) -c Universe.cpp

Body.o: Body.cpp Body.h
	$(CC) $(CFLAGS) -c Body.cpp

clean :
	rm -f *~ \#*# *.o *.x