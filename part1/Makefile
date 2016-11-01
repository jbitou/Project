OBJS 	= lsh.o chain.o hash.o distances.o nnrlist.o
SOURCE	= lsh.c chain.c hash.c distances.c nnrlist.c
HEADER  = lsh.h chain.h hash.h distances.h nnrlist.h
OUT  	= lsh
CC	    = gcc
FLAGS   = -g -c 	
# -g option enables debugging mode 
# -c flag generates object code for separate files

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) -lm

# create/compile the individual files >>separately<< 
lsh.o: lsh.c hash.h
	$(CC) $(FLAGS) lsh.c

chain.o: chain.c chain.h distances.h
	$(CC) $(FLAGS) chain.c
	
hash.o: hash.c hash.h
	$(CC) $(FLAGS) hash.c
	
distances.o: distances.c distances.h
	$(CC) $(FLAGS) distances.c
	
nnrlist.o: nnrlist.c nnrlist.h
	$(CC) $(FLAGS) nnrlist.c
	
# clean house
clean:
	rm -f $(OBJS) $(OUT)

