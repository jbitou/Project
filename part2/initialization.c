#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "initialization.h"

int binarySearch(int Nk, int search, int *array) {
	int first = 0;
	int last = Nk;
	int middle = (first+last)/2;
	if (search == 0)  return search;
	while (first <= last) {
		if (array[middle] < search)
			first = middle + 1;    
		else if (array[middle-1] < search && search <= array[middle])
			break;
		else
			last = middle - 1;
		middle = (first + last)/2;
	}
	return middle; 
}

pj_info *sortArray(pj_info *array, int N) {
	int i, j;
	pj_info a;
	for (i = 0; i < N; i++) {
		for (j = i + 1; j < N; j++) {
			if (array[i].v > array[j].v) {
				a =  array[i];
                array[i] = array[j];
                array[j] = a;
            }
        }
	}
	return array;
}

int *matrix_init_kmedoids(int **distances, pinfo info, int N) {
	int k = 1, i, j, min, max, x, z, flag;
	int *centroids = malloc((info->k)*sizeof(int)); 
	centroids[0] = (rand() / (RAND_MAX + 1.0)) * N;
	while (k < info->k) {
		int *D = malloc((N-k)*sizeof(int));
		int *P = malloc((N-k+1)*sizeof(int));
		z = 0;
		for (i=0; i < N; i++) {
			flag = 0;
			for (j=0; j < k; j++) {
				if (i == centroids[j]) flag = 1;
			}
			if (flag == 1) continue;
			if (i < centroids[0]) D[z] = distances[i][centroids[0]-i-1];
			else if (i > centroids[0]) D[z] = distances[centroids[0]][i-centroids[0]-1];
			for (j=1; j < k; j++) {
				if (i < centroids[j])	min = distances[i][centroids[j]-i-1];
				else if (i > centroids[j]) min = distances[centroids[j]][i-centroids[j]-1];
				if (min < D[z])	D[z] = min;
			}
			z++;
		}
		max = D[0];
		for (i=1; i < N-k; i++) {
			if (D[i] > max)  max = D[i];
		}
		P[0] = 0;
		for (i=1; i < N-k+1; i++) {
			P[i] = 0;
			for (j=0; j < i; j++) {
				P[i] += pow(D[j],2);
			}
			P[i] /= max;
		}
		x = (rand() / (RAND_MAX + 1.0)) * (P[N-k]+1);
		centroids[k] = binarySearch(N - k, x, P); //r
		free(D);
		free(P);
		k++;
	}
	return centroids;
}

int *matrix_init_concentrate(int **distances, pinfo info, int N) {
	int i, j, t, denominator;
	pj_info *v = malloc(N*sizeof(pj_info));
	/**For each object**/
	for (i=0; i < N; i++) {
		/**For outer Σ**/
		v[i].v = 0;
		v[i].index = i;
		for (j=0; j < N; j++) {
			/**For inner Σ**/
			denominator = 0;
			for (t=0; t < N; t++) {
				if (j < t)  denominator += distances[j][t-j-1];
				else if (j > t)  denominator += distances[t][j-t-1];
			}
			if (i < j)  {
				v[i].v += (double) distances[i][j-i-1] / denominator;
			}
			else if (i > j)  v[i].v += (double) distances[j][i-j-1] / denominator;
		}
	}
	v = sortArray(v, N);
	int *centroids = malloc((info->k)*sizeof(int));
	for (i=0; i < info->k; i++)  centroids[i] = v[i].index;	
	free(v);
	return centroids;
}
