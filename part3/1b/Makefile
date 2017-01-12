OBJS 	= main.o inputProcessing.o clustering.o cross-validation.o silhouette.o points_list.o distances.o
SOURCE	= main.c inputProcessing.c clustering.c cross-validation.c silhouette.c points_list.c distances.c
HEADER  = inputProcessing.h clustering.h cross-validation.h silhouette.h points_list.h distances.h
OUT  	= recommendation
CC	    = gcc
FLAGS   = -g -c 
# -g option enables debugging mode 
# -c flag generates object code for separate files

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) -lm

# create/compile the individual files >>separately<< 
main.o: main.c 
	$(CC) $(FLAGS) main.c

inputProcessing.o: inputProcessing.c inputProcessing.h
	$(CC) $(FLAGS) inputProcessing.c
	
cross-validation.o: cross-validation.c cross-validation.h
	$(CC) $(FLAGS) cross-validation.c	
	
clustering.o: clustering.c clustering.h
	$(CC) $(FLAGS) clustering.c	
	
silhouette.o: silhouette.c silhouette.h
	$(CC) $(FLAGS) silhouette.c	
	
distances.o: distances.c distances.h
	$(CC) $(FLAGS) distances.c	
	
points_list.o: points_list.c points_list.h
	$(CC) $(FLAGS) points_list.c	

# clean house
clean:
	rm -f $(OBJS) $(OUT)

