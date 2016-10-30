OBJS 	= chain.o hash.o distances.o nnrlist.o HashTest.o
SOURCE	= chain.c hash.c distances.c nnrlist.c HashTest.c
HEADER  = chain.h hash.h distances.h nnrlist.h
OUT  	= test
CC	    = gcc
FLAGS   = -c 	
# -c flag generates object code for separate files

all: $(OBJS)
	$(CC) -Wall -I$HOME/local/include  $(OBJS) -L$HOME/local/lib -lcunit -o $(OUT) -lm

#create/compile the individual files >>separately<< 
chain.o: chain.c chain.h distances.h
	$(CC) $(FLAGS) chain.c
	
hash.o: hash.c hash.h
	$(CC) $(FLAGS) hash.c
	
distances.o: distances.c distances.h
	$(CC) $(FLAGS) distances.c
	
nnrlist.o: nnrlist.c nnrlist.h
	$(CC) $(FLAGS) nnrlist.c
	
HashTest.o: HashTest.c hash.h
	$(CC) $(FLAGS) HashTest.c
	
#clean house
clean:
	rm -f $(OBJS) $(OUT)
