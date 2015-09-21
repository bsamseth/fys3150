CC = c++
CFLAGS = -g -Wall

OBJS = jacobiMethod.o


jacobiMethod: jacobiMethod.cpp
	$(CC) $(CFLAGS) jacobiMethod.cpp -o jacobiMethod.x 

unitTest: unitTest.cpp $(OBJS)
	$(CC) $(CFLAGS) unitTest.cpp -o unitTest.x 

clean :
	rm -f *~ \#*# $(OBJS)
