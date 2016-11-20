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
	chain neighbor;
	nnrp next;
}nnr;

void insert_chain(char *, void *, chainp *, int, int, int);
void search_chain_NNR(chainp, void *, double, nnrp *, int, int, int);
void destroy_chain(chainp *, int);
void insert_nnrlist(chainp, nnrp *);
//void combine_nnrlist(chainp *, nnrp *);
void print_nnrlist(nnrp *,FILE *);
void display_nnrlist(nnrp);
int make_item(char *);
