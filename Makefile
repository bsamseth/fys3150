CC = c++
CFLAGS = -g -Wall

OBJS = jacobiMethod.o


jacobiMethod: jacobiMethod.cpp
	$(CC) $(CFLAGS) jacobiMethod.cpp -o jacobiMethod.x 

clean :
	rm -f *~ \#*# $(OBJS)
