OBJS 	= main.o inputProcessing.o NNsearch.o hash.o chain.o distances.o 
SOURCE	= main.c inputProcessing.c NNsearch.c hash.c chain.c distances.c 
HEADER  = inputProcessing.h NNsearch.h hash.h chain.h distances.h 
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

hash.o: hash.c hash.h
	$(CC) $(FLAGS) hash.c

chain.o: chain.c chain.h
	$(CC) $(FLAGS) chain.c	
	
distances.o: distances.c distances.h
	$(CC) $(FLAGS) distances.c	
						
			
# clean house
clean:
	rm -f $(OBJS) $(OUT)

