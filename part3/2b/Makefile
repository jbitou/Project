OBJS 	= main.o inputProcessing.o distanceDRMSD.o
SOURCE	= main.c inputProcessing.c  distanceDRMSD.c
HEADER  = inputProcessing.h  distanceDRMSD.h
OUT  	= drmsd
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
		
distanceDRMSD.o: distanceDRMSD.c distanceDRMSD.h
	$(CC) $(FLAGS) distanceDRMSD.c	
			
# clean house
clean:
	rm -f $(OBJS) $(OUT)

