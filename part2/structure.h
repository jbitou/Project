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
int doublebinarySearch(int, double, double *);
pj_info *sortArray(pj_info *, int);
centroid *matrix_init_krandom(pinfo, int **, int);
centroid *vector_init_krandom(pinfo, double **, int);
centroid *matrix_init_kmedoids(int **, pinfo, int);
centroid *vector_init_kmedoids(double **, pinfo, int);
centroid *matrix_init_concentrate(int **, pinfo, int);
centroid *vector_init_concentrate(double **, pinfo, int);

/**Insert functions**/
hash_table *vector_insert_hash(hash_table *, ghashp *, pinfo, FILE *, int);
hash_table *matrix_insert_hash(hash_table *, ghashp *, int **, pinfo);
hash_table *hamming_insert_hash(hash_table *, ghashp *, FILE *, pinfo);
hash_table hamming_one_hash(hash_table, ghashp *, FILE *, pinfo, int *, int);
hash_table matrix_one_hash(hash_table, ghashp *, int **, int *, pinfo, int);
hash_table vector_one_hash(hash_table, ghashp *, pinfo, FILE *, int , int *, int);

/**Assignment functions**/
pcluster matrix_simplest_assignment(pcluster, int **, hash_table, centroid *, int);
pcluster vector_simplest_assignment(pcluster, double **, hash_table, centroid *, int);
pcluster matrix_reverse_approach(pcluster, int **, hash_table *, ghashp *, centroid *, pinfo, int);
pcluster vector_reverse_approach(pcluster, double **, hash_table *, ghashp *, centroid *, pinfo, int);
pcluster matrix_lsh_second_cluster(pcluster, int **, int);
pcluster vector_lsh_second_cluster(pcluster, double **, int);
pcluster remove_clusters_duplicates(pcluster, int);
pcluster matrix_assign_rest(pcluster, int **, hash_table *, ghashp *, chainp **, pinfo, int);
pcluster vector_assign_rest(pcluster, double **, hash_table *, ghashp *, chainp **, pinfo, int);
int matrix_compute_start_radius(int **, centroid *, int);
double vector_compute_start_radius(double **, centroid *, int);
double matrix_compute_objective_function(pcluster, int **, int);
double vector_compute_objective_function(pcluster, double **, int);

/**Update functions**/
centroid *vector_pam_update(pcluster, centroid *, double, double **, pinfo, int);
centroid *matrix_pam_update(pcluster, centroid *, double, int **, pinfo, int);
centroid *matrix_update_alaloyds(pcluster, centroid *, double, int **, pinfo);
centroid *vector_update_alaloyds(pcluster, centroid *, double, double **, pinfo);
centroid *matrix_update_clarans(pcluster, centroid *, cpair, int **, int, pinfo);
centroid *vector_update_clarans(pcluster, centroid *, cpair, double **, double, pinfo);
pointp matrix_calculate_medoid(pointp, int **);
pointp vector_calculate_medoid(pointp, double **);
int compare_centroids(centroid *, centroid *, int);
cpair matrix_select_pairs(centroid *, int **, pinfo);
cpair vector_select_pairs(centroid *, double **, pinfo);
cpair sortPairs(cpair, int);
