
/**Centroids info**/
typedef struct centroid_node {
	int center;
	double *vector;
}centroid;

/**Points info**/
typedef struct point_node *pointp;
typedef struct point_node {
	int position;
	centroid second;
	double mindistance;
	double secdistance;
	pointp next;
}point;

int insert_points(pointp *, double, double, centroid, int);
pointp clone_points(pointp);
int points_length(pointp);
void destroy_points(pointp *);
void printndestroy_points(pointp *, FILE *);
void print_points(pointp);
