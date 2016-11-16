#include <stdint.h>
#include <inttypes.h>
#include "nnrlist.h"

typedef struct chain_node *chainp;
typedef struct chain_node
{
	char *key;
	uint64_t *value;
	double *p;
	int id;
	int *distances;
	chainp next;
}chain;

typedef struct nn_node
{
	char *key;
	double distance;
}nn;

int make_item(char *item);
void insert_chain(char *, void *, chainp *, int, int, int);
void search_chain_NNR(chainp, void *, double, nnrp *, int, int, int);
void search_chain_NN(chainp, void *, int, int, int, int, int *, int, nnrp*, double *);
void destroy_chain(chainp *, int);
