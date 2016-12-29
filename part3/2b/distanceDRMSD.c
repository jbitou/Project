#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "distanceDRMSD.h"

dinfo *create_distances(double **data, int numConform, int N, int r, int T) {
	int i, j, flag, index;
	dinfo **alldistances, *rdis;
	/**Get all distances**/
	alldistances = create_vectors(data,numConform,N);
	/**Choose r distances**/
	if (T == 3) {
		rdis = malloc(r*sizeof(dinfo));
		for (i=0; i < r; i++) {
			flag = 0;
			index = (rand() / (RAND_MAX + 1.0)) * (N*(N-1)/2);
			printf("index=%d\n",index);	
			for (j=0; j < i; j++) {
				if ((alldistances[0][index].point1 == rdis[j].point1) && (alldistances[0][index].point2 == rdis[j].point2)) {
					i--;
					flag = 1;
				}
				if (flag) break;
			}
			if (!flag)	{
				rdis[i] = alldistances[0][index];
				printf("rdis[%d] -> p1 = %d, p2 = %d, distance = %lf\n",i,rdis[i].point1,rdis[i].point2,rdis[i].distance);
			}
		}
	}
	else {
	}
	for (i=0; i < numConform; i++)	free(alldistances[i]);
	free(alldistances);
	return rdis;
}

dinfo **create_vectors(double **data, int numConform, int N) {
	int i, j, z, y, start;
	double *v1, *v2;
	dinfo **alldistances;
	alldistances = malloc(numConform*sizeof(dinfo *));
	/**For each conformation**/
	for (i=0; i < numConform; i++) {
		printf("Conform %d\n",i+1);
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
				printf("item %d distance with item %d\n",j,z);
				v2[0] = data[z][0];
				v2[1] = data[z][1];
				v2[2] = data[z][2];
				alldistances[i][y].point1 = j;
				alldistances[i][y].point2 = z;
				alldistances[i][y].distance = distance_Euclidean(v1,v2,3);
				printf("distance = %lf\n",alldistances[i][y].distance);
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
