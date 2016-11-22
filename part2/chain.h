#include <stdint.h>
#include <inttypes.h>

typedef struct chain_node *chainp;
typedef struct chain_node
{
	char *key;
	uint64_t *value;
	/**distance: for part2, LSH assignment**/
	double *p, distance;
	int id;
	chainp next;
}chain;

void insert_chain(char *, void *, chainp *, double, int, int, int);
int search_chain_NNR(chainp *, void *, double, chainp *, chainp *, int, int, int, int *);
void move_chain_nodes(chainp *, chainp);
void print_chain(chainp);
void destroy_chain(chainp *, int);
void delete_from_chain(chainp *, char *);
//void print_nnrlist(nnrp *,FILE *);
int make_item(char *);
