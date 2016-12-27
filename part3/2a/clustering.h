#include "points_list.h"

/**Clusters info**/
typedef struct cluster_node *pcluster;
typedef struct cluster_node {
	pointp items;
	int center;
}cluster;

pcluster clustering(double **, int, int, int);
int *vector_init_kmedoids(double **, int, int, int);
int doublebinarySearch(int, double, double *);
pcluster vector_simplest_assignment(pcluster, double **, int *, int, int, int);
double vector_compute_objective_function(pcluster, double **, int, int);
int *vector_update_alaloyds(pcluster, int *, double, double **, int, int);
pointp vector_calculate_medoid(pointp, double **, int);
