FLAGS = -O3 -w
DEBUG = -g -p
CC = g++
SRCS = ljmd.cpp

OBJS = $(SRCS)

main: $(OBJS) 
	$(CC)  $(OBJS) $(FLAGS) -o run

debug: $(OBJS)
	$(CC)  $(OBJS) $(DEBUG) -o run

%.cpp.o: %.cpp
	$(CC) $(FLAGS) $(DEBUG) -c -o $@ $< 
clean:
	rm *.cpp.o
