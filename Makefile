CC = c++
CFLAGS = -g -Wall -std=c++11

OBJS = jacobiMethod.o 
LIBS = -lunittest++


interactingElectrons: interactingElectrons.o
	$(CC) $(CFLAGS) $(OBJS) interactingElectrons.o -o interactingElectrons.x

interactingElectrons.o: interactingElectrons.cpp
	$(CC) $(CFLAGS) -c interactingElectrons.cpp

solve_lower_3_states: solve_lower_3_states.o
	$(CC) $(CFLAGS) $(OBJS) solve_lower_3_states.o -o solve_lower_3_states.x

solve_lower_3_states.o: solve_lower_3_states.cpp
	$(CC) $(CFLAGS) -c solve_lower_3_states.cpp

unitTest: unitTest.cpp $(OBJS) unitTest.o
	$(CC) $(CFLAGS) $(OBJS) unitTest.o -o unitTest.x ${LIBS}

unitTest.o: unitTest.cpp 
	$(CC) $(CFLAGS) -c unitTest.cpp 

jacobiMethod.o: jacobiMethod.cpp
	$(CC) $(CFLAGS) -c jacobiMethod.cpp

clean :
	rm -f *~ \#*# *.o

