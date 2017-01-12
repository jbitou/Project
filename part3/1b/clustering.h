#include "points_list.h"
#include "distances.h"

/**Clusters info**/
typedef struct cluster_node *pcluster;
typedef struct cluster_node {
	pointp items;
	int center;
}cluster;

pcluster k_clustering(user *, int, int, int, int *, double *);
pcluster clustering(user *, int, int, int);
int *init_krandom(int, int);
pcluster simplest_assignment(pcluster, user *, int *, int, int, int);
double compute_objective_function(pcluster, user *, int, int);
int *update_alaloyds(pcluster, int *, double, user *, int, int);
pointp calculate_medoid(pointp, user *, int);
int compare_centroids(int *, int *, int);
