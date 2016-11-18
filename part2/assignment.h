#include "hash.h"

typedef struct cluster_node *pcluster;
typedef struct cluster_node {
	chainp items;
	void *centroid;
}cluster;

hash_table *matrix_insert_hash(hash_table *, ghashp *, int **, int, int, int);
pcluster matrix_simplest_assignment(pcluster, int **, hash_table, int *, int);
int compute_objective_function(pcluster, int **, int);
