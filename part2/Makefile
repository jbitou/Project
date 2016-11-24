OBJS 	= medoids.o chain.o hash.o distances.o matrixMedoid.o inputProcessing.o initialization.o assignment.o update.o
SOURCE	= medoids.c chain.c hash.c distances.c matrixMedoid.c inputProcessing.c initialization.c assignment.c update.c
HEADER  = medoids.h chain.h hash.h distances.h matrixMedoid.h inputProcessing.h structure.h
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

inputProcessing.o: inputProcessing.c inputProcessing.h
	$(CC) $(FLAGS) inputProcessing.c
	
initialization.o: initialization.c structure.h
	$(CC) $(FLAGS) initialization.c
	
assignment.o: assignment.c structure.h
	$(CC) $(FLAGS) assignment.c	
	
update.o: update.c structure.h
	$(CC) $(FLAGS) update.c	

matrixMedoid.o: matrixMedoid.c matrixMedoid.h
	$(CC) $(FLAGS) matrixMedoid.c		
	
# clean house
clean:
	rm -f $(OBJS) $(OUT)

