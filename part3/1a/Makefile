OBJS 	= main.o inputProcessing.o NNsearch.o
SOURCE	= main.c inputProcessing.c NNsearch.c
HEADER  = inputProcessing.h NNsearch.h
OUT  	= recommendation
CC	    = gcc
FLAGS   = -g -c 
# -g option enables debugging mode 
# -c flag generates object code for separate files

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) -lm

# create/compile the individual files >>separately<< 
main.o: main.c inputProcessing.h 
	$(CC) $(FLAGS) main.c

inputProcessing.o: inputProcessing.c inputProcessing.h
	$(CC) $(FLAGS) inputProcessing.c	

NNsearch.o: NNsearch.c NNsearch.h
	$(CC) $(FLAGS) NNsearch.c	
			
# clean house
clean:
	rm -f $(OBJS) $(OUT)

