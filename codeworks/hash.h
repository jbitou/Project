#include "chain.h"

typedef struct hash_table 
{
	chainp *table;
	int size;
}hash_table;

typedef struct g_node * ghashp;
typedef struct g_node
{
	int t, r, t1;
	double *v;
}ghash;

int mod(int, long long);
void init_table(int, hash_table *, int);
void init_hash_Ham(ghashp *, int, int, char*);
void init_hash_Eucl(ghashp *, int, int, int);
void init_hash_Cos(ghashp *, int, int, int);
void init_hash_matrix(ghashp *, int **, int, int, int);
int hash_func_Ham(ghashp, char *, int);
int hash_func_Eucl(ghashp, double *, int, int);
int hash_func_Cos(ghashp, double *, int, int);
int hash_func_Matrix(ghashp, int, int **, int, int);
int hash_func_MSearch(ghashp, int *, int **, int, int);
void search_table_NNR(int, hash_table, void *, double, nnrp *, int, int, int);
nn search_table_NN(ghashp *, hash_table *, void *, int **, int, int, int, int, int);
nn brute_force_table(hash_table, void *, int, int, int, int, int);
void destroy_table(hash_table *, int);
