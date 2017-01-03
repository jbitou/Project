#include "points_list.h"
#include "distanceDRMSD.h"


/**Clusters info**/
typedef struct cluster_node *pcluster;
typedef struct cluster_node {
	pointp items;
	centroid center;
}cluster;


void experiment(double **, dinfo **, int, int, int, int, int, FILE *);
pcluster k_clustering(double **, int, int, int, int *, double *);
pcluster clustering(double **, int, int, int);
centroid *vector_init_kmeans(double **, int, int, int);
int doublebinarySearch(int, double, double *);
pcluster vector_simplest_assignment(pcluster, double **, centroid *, int, int, int);
double vector_compute_objective_function(pcluster, double **, int, int);
centroid *vector_update_loyds(pcluster, centroid *, double, double **, int, int);
double *calculate_mean(pointp, double **, int);
int compare_centroids(centroid *, centroid *, int, int);
