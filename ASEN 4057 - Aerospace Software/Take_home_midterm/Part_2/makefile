OBJS	= Problem2dot1.o
SOURCE	= Problem2dot1.c
HEADER	=
OUT	= Problem2dot1
CC	 = gcc
FLAGS	 = -g -c -Wall
LFLAGS	 =

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)

Problem2dot1.o: Problem2dot1.c
	$(CC) $(FLAGS) Problem2dot1.c


clean:
	rm -f $(OBJS) $(OUT)
