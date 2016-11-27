OBJS 	= medoids.o chain.o hash.o distances.o matrixMedoid.o inputProcessing.o initialization.o assignment.o update.o silhouette.o algorithms.o vectorMedoid.o
SOURCE	= medoids.c chain.c hash.c distances.c matrixMedoid.c inputProcessing.c initialization.c assignment.c update.c silhouette.c algorithms.c vectorMedoid.c
HEADER  = medoids.h chain.h hash.h distances.h medoids.h inputProcessing.h structure.h silhouette.h algorithms.h
OUT  	= medoids
CC	    = gcc
FLAGS   = -g -c 	
# -g option enables debugging mode 
# -c flag generates object code for separate files

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) -lm

# create/compile the individual files >>separately<< 
medoids.o: medoids.c medoids.h
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
	
silhouette.o : silhouette.c silhouette.h
	$(CC) $(FLAGS) silhouette.c

algorithms.o : algorithms.c algorithms.h
	$(CC) $(FLAGS) algorithms.c

matrixMedoid.o: matrixMedoid.c medoids.h
	$(CC) $(FLAGS) matrixMedoid.c	
	
vectorMedoid.o: vectorMedoid.c medoids.h
	$(CC) $(FLAGS) vectorMedoid.c
	
# clean house
clean:
	rm -f $(OBJS) $(OUT)

