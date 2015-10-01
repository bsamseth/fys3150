CC = c++
CFLAGS = -g -Wall

OBJS = jacobiMethod.o 
LIBS = -lunittest++


interactingElectrons: interactingElectrons.o
	$(CC) $(CFLAGS) $(OBJS) interactingElectrons.o -o interactingElectrons.x

interactingElectrons.o: interactingElectrons.cpp
	$(CC) $(CFLAGS) -c interactingElectrons.cpp

findnstep: findnstep.o
	$(CC) $(CFLAGS) $(OBJS) findnstep.o -o findnstep.x 

findnstep.o: jacobiMethod.o
	$(CC) $(CFLAGS) -c findnstep.cpp

unitTest: unitTest.cpp $(OBJS) unitTest.o
	$(CC) $(CFLAGS) $(OBJS) unitTest.o -o unitTest.x ${LIBS}

unitTest.o: unitTest.cpp 
	$(CC) $(CFLAGS) -c unitTest.cpp 

jacobiMethod.o: jacobiMethod.cpp
	$(CC) $(CFLAGS) -c jacobiMethod.cpp

clean :
	rm -f *~ \#*# *.o
