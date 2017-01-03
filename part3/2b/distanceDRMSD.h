
typedef struct distance_info {
	int point1;
	int point2;
	double distance;
}dinfo;

double **create_vectors(dinfo **, dinfo *, int, int, int);
dinfo *create_distances(dinfo *, double **, int, int, int, int);
dinfo **get_all_distances(double **, int, int);
double distance_Euclidean(double *, double *, int);
dinfo *sortdistances(dinfo *, int);
