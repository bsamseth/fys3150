CC = g++
CFLAGS = -O3 -Wall -std=c++0x
LFLAGS = 

OBJS = Body.o Universe.o
LIBS = ../../lib/lib.o
TARGETS = main 2body probabilityDensity

all: $(TARGETS)

main: $(OBJS) main.cpp
	$(CC) $(CFLAGS) -Wno-unused-variable $(OBJS) $(LIBS) -o main.x main.cpp

2body: $(OBJS) 2body.cpp
	$(CC) $(CFLAGS) -Wno-unused-variable $(OBJS) $(LIBS) -o 2body.x 2body.cpp

Universe.o: Body.o Universe.cpp Universe.h
	$(CC) $(CFLAGS) -c Universe.cpp

Body.o: Body.cpp Body.h
	$(CC) $(CFLAGS) -c Body.cpp

2body: $(OBJS) 2body.cpp
	$(CC) $(CFLAGS) -Wno-unused-variable $(OBJS) $(LIBS) -o 2body.x 2body.cpp

probabilityDensity: $(OBJS) probabilityDensity.cpp
	$(CC) $(CFLAGS) -Wno-unused-variable $(OBJS) $(LIBS) -o probabilityDensity.x probabilityDensity.cpp

clean :
	rm -f *~ \#*# *.o *.x
