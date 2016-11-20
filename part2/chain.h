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

void insert_chain(char *, void *, chainp *, int, int, int);
void search_chain_NNR(chainp, void *, double, chainp *, int, int, int);
void print_chain(chainp);
void destroy_chain(chainp *, int);
void insert_list(char *, chainp *, int);
//void print_nnrlist(nnrp *,FILE *);
int make_item(char *);
