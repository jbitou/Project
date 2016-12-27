OBJS 	= main.o inputProcessing.o distanceCRMSD.o clustering.o points_list.o silhouette.o
SOURCE	= main.c inputProcessing.c distanceCRMSD.c clustering.c points_list.c silhouette.c
HEADER  = inputProcessing.h distanceCRMSD.h lapacke.h cblas.h clustering.h points_list.h silhouette.h
OUT  	= crmsd
CC	    = gcc
FLAGS   = -g -c 
# -g option enables debugging mode 
# -c flag generates object code for separate files

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) -lm -lblas -llapacke

# create/compile the individual files >>separately<< 
main.o: main.c inputProcessing.h silhouette.h
	$(CC) $(FLAGS) main.c

inputProcessing.o: inputProcessing.c inputProcessing.h
	$(CC) $(FLAGS) inputProcessing.c
	
distanceCRMSD.o: distanceCRMSD.c distanceCRMSD.h lapacke.h cblas.h
	$(CC) $(FLAGS) distanceCRMSD.c

clustering.o: clustering.c clustering.h
	$(CC) $(FLAGS) clustering.c	
		
points_list.o: points_list.c points_list.h
	$(CC) $(FLAGS) points_list.c

silhouette.o: silhouette.c silhouette.h clustering.h
	$(CC) $(FLAGS) silhouette.c		
			
# clean house
clean:
	rm -f $(OBJS) $(OUT)

