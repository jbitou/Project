#include "distances.h"

/**Bucket of hash tables**/
typedef struct chain_node *chainp;
typedef struct chain_node {
	int euclid;
	user usr;
	chainp next;
}chain;

/**Nearest neighbours**/
typedef struct nnr_node *nnrp;
typedef struct nnr_node {
	double distance;
	user usr;
	nnrp next;
}nnr;

typedef struct nnr_info {
	int length;
	nnrp list;
}nnrlist;


void vector_insert_chain(chainp *, user, int, int); 
nnrlist search_chain_NNR(chainp, nnrlist, user, double, double, int, int); 
void destroy_chain(chainp *);
void print_chain(chainp);

int insert_nnrlist(nnrp *, user, int, double);
void combine_nnrlist(nnrlist *, nnrlist *, int, int);
int nnrlist_length(nnrp);
void destroy_nnrlist(nnrp *);
void print_nnrlist(nnrp);
