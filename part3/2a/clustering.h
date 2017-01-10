#include "points_list.h"

/**Clusters info**/
typedef struct cluster_node *pcluster;
typedef struct cluster_node {
	pointp items;
	int center;
}cluster;

pcluster k_clustering(double **, int, int, int, int *, double *);
pcluster clustering(double **, int, int, int);
int *init_krandom(int, int);
int *init_kmedoids(double **, int, int, int);
int binarySearch(int, double, double *);
pcluster simplest_assignment(pcluster, double **, int *, int, int, int);
double compute_objective_function(pcluster, double **, int, int);
int *update_alaloyds(pcluster, int *, double, double **, int, int);
pointp calculate_medoid(pointp, double **, int);
int compare_centroids(int *, int *, int);
