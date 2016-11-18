#include <stdint.h>
#include <inttypes.h>

typedef struct chain_node *chainp;
typedef struct chain_node
{
	char *key;
	uint64_t *value;
	double *p;
	int id;
	chainp next;
}chain;

typedef struct nnr_node *nnrp;
typedef struct nnr_node
{
	chain *neighbor;
	nnrp next;
}nnr;

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
void insert_nnrlist(chain *, nnrp *);
void destroy_nnrlist(nnrp *);
void print_nnrlist(nnrp *,FILE *);
void combine_nnrlist(nnrp *, nnrp *);
void display_nnrlist(nnrp);
