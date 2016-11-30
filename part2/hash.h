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
int *sortMatrix(int *, int);
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
int **create_small_distance_matrix(int **, int, int, int *);
int **create_ham_distance_table(hash_table, int);
double **create_vector_distance_table(hash_table, int, int, int);
int search_table_NNR(int, hash_table *, void *, double, pointp *, chainp *, int, int, int, int *);
double *find_vector_info(hash_table, int);
char *find_ham_info(hash_table, int);
void destroy_table(hash_table *, int);
