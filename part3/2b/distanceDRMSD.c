#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "distanceDRMSD.h"

double **create_vectors(dinfo **alldistances, dinfo *distances, int numConform, int N, int r) {
	int i, j, z;
	double **vectors;
	vectors = malloc(numConform*sizeof(double *));
	vectors[0] = malloc(r*sizeof(double));
	for (i=0; i < r; i++) vectors[0][i] = distances[i].distance;
	for (i=1; i < numConform; i++) {
		vectors[i] = malloc(r*sizeof(double));
		for (j=0; j < r; j++) {
			for (z=0; z < N*(N-1)/2; z++) {
				if ((alldistances[i][z].point1 == distances[j].point1) && (alldistances[i][z].point2 == distances[j].point2)) {
					vectors[i][j] = alldistances[i][z].distance;
					break;
				}
			}
		}
	}
	return vectors;
}	

dinfo *create_distances(dinfo *alldistances, double **data, int numConform, int N, int r, int T) {
	int i, j, flag, index, size;
	dinfo *rdis, *temp;
	/**Choose r distances**/
	size = N*(N-1)/2;
	if (r == size) {
		temp =  malloc(size*sizeof(dinfo));
		for (i=0; i < size; i++)	temp[i] = alldistances[i];
		return temp;
	}
	rdis = malloc(r*sizeof(dinfo));
	if (T == 3) {
		for (i=0; i < r; i++) {
			flag = 0;
			index = (rand() / (RAND_MAX + 1.0)) * size;	
			/**Choose unique distances**/
			for (j=0; j < i; j++) {
				if ((alldistances[index].point1 == rdis[j].point1) && (alldistances[index].point2 == rdis[j].point2)) {
					i--;
					flag = 1;
				}
				if (flag) break;
			}
			if (!flag)	rdis[i] = alldistances[index];
		}
	}
	else {
		alldistances = sortdistances(alldistances,size);
		if (T == 1) {
			for (i=0; i < r; i++) rdis[i] = alldistances[i];
		}
		if (T == 2) {
			j = 0;
			for (i=size-1; i >= size-r; i--) {
				rdis[j] = alldistances[i];
				j++;
			}
		}
	}
	return rdis;
}

dinfo **get_all_distances(double **data, int numConform, int N) {
	int i, j, z, y, start;
	double *v1, *v2;
	dinfo **alldistances;
	alldistances = malloc(numConform*sizeof(dinfo *));
	/**For each conformation**/
	for (i=0; i < numConform; i++) {
		alldistances[i] = malloc((N*(N-1)/2)*sizeof(dinfo));
		v1 = malloc(3*sizeof(double));
		start = i*N;
		y = 0;
		/**Calculate all distances from the others**/
		for (j=start; j < start+N; j++) {
			v1[0] = data[j][0];
			v1[1] = data[j][1];
			v1[2] = data[j][2];
			v2 = malloc(3*sizeof(double));
			for (z=j+1; z < start+N; z++) {
				v2[0] = data[z][0];
				v2[1] = data[z][1];
				v2[2] = data[z][2];
				alldistances[i][y].point1 = abs(start-j);
				alldistances[i][y].point2 = abs(start-z);
				alldistances[i][y].distance = distance_Euclidean(v1,v2,3);
				y++;
			}
			free(v2);
		}
		free(v1);
	}
	return alldistances;
}

/**Returns the euclidean difference between two vectors**/
double distance_Euclidean(double *v1, double *v2, int d) {
	int i;
	double distance = 0.0, diff = 0.0 ;
	for (i=0; i < d; i++) {
		diff = v1[i] - v2[i];
		distance += diff*diff;
	}
	return sqrt(distance);
}


dinfo *sortdistances(dinfo *array, int N) {
	int i, j;
	dinfo a;
	for (i = 0; i < N; i++) {
		for (j = i + 1; j < N; j++) {
			if (array[i].distance > array[j].distance) {
				a =  array[i];
                array[i] = array[j];
                array[j] = a;
            }
        }
	}
	return array;
}
