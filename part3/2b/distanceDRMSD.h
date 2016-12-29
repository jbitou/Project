
typedef struct distance_info {
	int point1;
	int point2;
	double distance;
}dinfo;

dinfo *create_distances(double **, int, int, int, int);
dinfo **create_vectors(double **, int, int);
double distance_Euclidean(double *, double *, int);
