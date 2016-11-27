#include <stdint.h>
#include <inttypes.h>
#include "distances.h"

/**Chains info**/
typedef struct chain_node *chainp;
typedef struct chain_node
{
	char *key;
	uint64_t *value;
	double *p;
	int id;
	int position;
	chainp next;
}chain;

/**Centroids info**/
typedef struct centroid_node {
	void *center;
	void *info;
}centroid;

/**Points info**/
typedef struct point_node *pointp;
typedef struct point_node {
	char *key;
	int position;
	int duplicate;
	double mindistance, secdistance;
	centroid second;
	pointp next;
}point;

/**Clarans pairs info**/
typedef struct clarans_pair *cpair;
typedef struct clarans_pair {
	centroid m;
	centroid t;
}pair;

void insert_chain(char *, void *, chainp *, int, int, int, int);
int search_chain_NNR(chainp *, void *, double, pointp *, chainp *, int, int, int, int *);
void move_chain_nodes(chainp *, chainp);
void print_chain(chainp);
void destroy_chain(chainp *, int);
int insert_points(pointp *, char *, double, double, centroid, int);
int chain_length(pointp);
void print_points(pointp);
void delete_from_chain(pointp *, char *);
void destroy_points(pointp *);
pointp clone(pointp);
int make_item(char *);
