#include "inputProcessing.h" 
#include "hash.h"

/**Input info**/
typedef struct park_jun_info {
	int index;
	double v;
}pj_info;

/**Clusters info**/
typedef struct cluster_node *pcluster;
typedef struct cluster_node {
	pointp items;
	centroid center;
}cluster;

/**Initialization functions**/
int binarySearch(int, int, int *);
pj_info *sortArray(pj_info *, int);
centroid *matrix_init_kmedoids(int **, pinfo, int);
centroid *matrix_init_concentrate(int **, pinfo, int);

/**Assignment functions**/
hash_table *matrix_insert_hash(hash_table *, ghashp *, int **, pinfo);
pcluster matrix_simplest_assignment(pcluster, int **, hash_table, centroid *, int);
pcluster matrix_reverse_approach(pcluster, int **, hash_table *, ghashp *, centroid *, pinfo);
pcluster lsh_second_cluster(pcluster, int **, int);
pcluster matrix_remove_clusters_duplicates(pcluster, int);
pcluster matrix_assign_rest(pcluster, int **, hash_table *, ghashp *, chainp **, pinfo);
int matrix_compute_start_radius(int **, centroid *, int);
int matrix_compute_objective_function(pcluster, int **, int);

/**Update functions**/
centroid *matrix_update_alaloyds(pcluster, centroid *, int, int **, pinfo);
pointp matrix_calculate_medoid(pointp, int **);
