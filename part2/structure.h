#include "inputProcessing.h" 
#include "hash.h"

/**Input info**/
typedef struct park_jun_info {
	int index;
	double v;
}pj_info;

/**Centroids info**/
typedef struct centroid_node {
	void *center;
	void *info;
}centroid;

/**Clusters info**/
typedef struct cluster_node *pcluster;
typedef struct cluster_node {
	chainp items;
	centroid center;
}cluster;

/**Initialization functions**/
int binarySearch(int, int, int *);
pj_info *sortArray(pj_info *, int);
centroid *matrix_init_kmedoids(int **, pinfo, int);
centroid *matrix_init_concentrate(int **, pinfo, int);

/**Assignment functions**/
hash_table *matrix_insert_hash(hash_table *, ghashp *, int **, int, int, int);
pcluster matrix_simplest_assignment(pcluster, int **, hash_table, centroid *, int);
pcluster matrix_reverse_approach(pcluster, int **, hash_table *, ghashp *, centroid *, int, int, int, int);
int matrix_compute_start_radius(int **, centroid *, int);
int matrix_compute_objective_function(pcluster, int **, int);
