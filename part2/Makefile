OBJS 	= medoids.o chain.o hash.o distances.o nnrlist.o matrixMedoid.o inputProcessing.o initialization.o assignment.o
SOURCE	= medoids.c chain.c hash.c distances.c nnrlist.c matrixMedoid.c inputProcessing.c initialization.c assignment.c
HEADER  = medoids.h chain.h hash.h distances.h nnrlist.h matrixMedoid.h inputProcessing.h initialization.h assignment.h
OUT  	= medoids
CC	    = gcc
FLAGS   = -g -c 	
# -g option enables debugging mode 
# -c flag generates object code for separate files

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) -lm

# create/compile the individual files >>separately<< 
medoids.o: medoids.c matrixMedoid.h
	$(CC) $(FLAGS) medoids.c

chain.o: chain.c chain.h distances.h
	$(CC) $(FLAGS) chain.c
	
hash.o: hash.c hash.h
	$(CC) $(FLAGS) hash.c
	
distances.o: distances.c distances.h
	$(CC) $(FLAGS) distances.c
	
nnrlist.o: nnrlist.c nnrlist.h
	$(CC) $(FLAGS) nnrlist.c

inputProcessing.o: inputProcessing.c inputProcessing.h
	$(CC) $(FLAGS) inputProcessing.c
	
initialization.o: initialization.c initialization.h
	$(CC) $(FLAGS) initialization.c
	
assignment.o: assignment.c assignment.h
	$(CC) $(FLAGS) assignment.c	

matrixMedoid.o: matrixMedoid.c matrixMedoid.h hash.h
	$(CC) $(FLAGS) matrixMedoid.c		
	
# clean house
clean:
	rm -f $(OBJS) $(OUT)

