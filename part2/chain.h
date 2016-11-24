#include <stdint.h>
#include <inttypes.h>

typedef struct chain_node *chainp;
typedef struct chain_node
{
	char *key;
	uint64_t *value;
	double *p;
	int id;
	chainp next;
}chain;

/**Centroids info**/
typedef struct centroid_node {
	void *center;
	void *info;
}centroid;

typedef struct point_node *pointp;
typedef struct point_node {
	char *key;
	int duplicate;
	double mindistance, secdistance;
	centroid second;
	pointp next;
}point;

void insert_chain(char *, void *, chainp *, int, int, int);
int search_chain_NNR(chainp *, void *, double, pointp *, chainp *, int, int, int, int *);
void move_chain_nodes(chainp *, chainp);
void print_chain(chainp);
void destroy_chain(chainp *, int);
int insert_points(pointp *, char *, double, double, centroid);
int chain_length(pointp);
void print_points(pointp);
void delete_from_chain(pointp *, char *);
void destroy_points(pointp *);
int make_item(char *);
