#include "chain.h"

typedef struct hash_table {
	chainp *table;
	int size;
}hash_table;

typedef struct g_node * ghashp;
typedef struct g_node {
	int t, r, t1;
	double *v;
}ghash;

int mod(int, long long);
void init_table(hash_table *, int);
//void init_hash_Ham(ghashp *, int, int, char *);
void init_hash_Eucl(ghashp *, int);
void init_hash_Cos(ghashp *, int);
//int hash_func_Ham(ghashp, char *, int);
int hash_func_Eucl(ghashp, user, int);
int hash_func_Cos(ghashp, user, int);
hash_table *vector_hash_insert(ghashp *, user *, int, int, int);
nnrlist search_table_NNR(hash_table *, ghashp *, user, double, double, int, int, int);
void destroy_table(hash_table *);
