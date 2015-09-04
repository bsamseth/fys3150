CC = g++
CFLAGS = -std=c++11 -g -Wall

OBJS = 


main: main.cpp
	$(CC) $(CFLAGS) $(OBJS) main.cpp -o main.x

clean :
	rm -f *~ \#*# $(OBJS)
