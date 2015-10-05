CC = g++
CFLAGS = -std=c++11 -g -Wall

OBJS = 


main: main.cpp
	$(CC) $(CFLAGS) $(OBJS) main.cpp -o main.x 
d: d.cpp
	$(CC) $(CFLAGS) $(OBJS) d.cpp -o d.x -larmadillo -llapack

clean :
	rm -f *~ \#*# $(OBJS)
