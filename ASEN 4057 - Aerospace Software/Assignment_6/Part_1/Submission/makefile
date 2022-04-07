OBJS	= ThreeBody.o
SOURCE	= ThreeBody.c
HEADER	=
OUT	= Homework6
CC	 = gcc
FLAGS	 = -g -c -Wall
LFLAGS	 = -lm

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)

ThreeBody.o: ThreeBody.c
	$(CC) $(FLAGS) ThreeBody.c


clean:
	rm -f $(OBJS) $(OUT)

