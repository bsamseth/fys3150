CC = c++
CFLAGS = -g -Wall

OBJS = jacobiMethod.o 
LIBS = /usr/lib/libgtest.a /usr/lib/libgtest_main.a -lpthread

unitTest: unitTest.cpp $(OBJS) unitTest.o
	$(CC) $(CFLAGS) $(OBJS) unitTest.o -o unitTest.x ${LIBS}

unitTest.o: unitTest.cpp 
	$(CC) $(CFLAGS) -c unitTest.cpp 

jacobiMethod.o: jacobiMethod.cpp
	$(CC) $(CFLAGS) -c jacobiMethod.cpp

clean :
	rm -f *~ \#*# *.o
