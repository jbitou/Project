OBJS 	= distances.o DistancesTest.o
SOURCE	= distances.c DistancesTest.c
HEADER  = distances.h 
OUT  	= test
CC	    = gcc
FLAGS   = -c 	
# -c flag generates object code for separate files

all: $(OBJS)
	$(CC) -Wall -I$HOME/local/include  $(OBJS) -L$HOME/local/lib -lcunit -o $(OUT) -lm

#create/compile the individual files >>separately<< 
distances.o: distances.c distances.h
	$(CC) $(FLAGS) distances.c
	
DistancesTest.o: DistancesTest.c 
	$(CC) $(FLAGS) DistancesTest.c
	
#clean house
clean:
	rm -f $(OBJS) $(OUT)
