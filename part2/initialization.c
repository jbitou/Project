#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "structure.h"

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

centroid *matrix_init_kmedoids(int **distances, pinfo info, int N) {
	int k = 1, i, j, x, z, min, max, flag, *temp;
	centroid *centroids = malloc((info->k)*sizeof(centroid)); 
	for (i=0; i < info->k; i++) 	centroids[i].info = malloc(N*sizeof(int));
	centroids[0].center = (void *)(intptr_t)((rand() / (RAND_MAX + 1.0)) * N);
	/**Create k centroids**/
	while (k < info->k) {
		int *D = malloc((N-k)*sizeof(int));
		int *P = malloc((N-k+1)*sizeof(int));
		z = 0;
		/**For each item**/
		for (i=0; i < N; i++) {
			flag = 0;
			for (j=0; j < k; j++) {
				if (i == (int)(intptr_t)centroids[j].center) flag = 1;
			}
			if (flag == 1) continue;
			if (i < (int)(intptr_t)centroids[0].center) 		D[z] = distances[i][(int)(intptr_t)centroids[0].center-i-1];
			else if (i > (int)(intptr_t)centroids[0].center) 	D[z] = distances[(int)(intptr_t)centroids[0].center][i-(int)(intptr_t)centroids[0].center-1];
			for (j=1; j < k; j++) {
				if (i < (int)(intptr_t)centroids[j].center)			min = distances[i][(int)(intptr_t)centroids[j].center-i-1];
				else if (i > (int)(intptr_t)centroids[j].center) 	min = distances[(int)(intptr_t)centroids[j].center][i-(int)(intptr_t)centroids[j].center-1];
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
			for (j=0; j < i; j++)	P[i] += pow(D[j],2);
			P[i] /= max;
		}
		x = (rand() / (RAND_MAX + 1.0)) * (P[N-k]+1);
		centroids[k].center = (void *)(intptr_t)binarySearch(N - k, x, P); //r
		for (i=0; i < k; i++) {
			if (centroids[k].center == centroids[i].center) {
				k--;
				break;
			}
		}		
		free(D);
		free(P);
		k++;
	}
	/**For each centroid**/
	for (i=0; i < info->k; i++) {
		temp = (int *)centroids[i].info;
		/**Create info line with distances**/
		for (j=0; j < (int)(intptr_t)centroids[i].center; j++) 
			temp[j] = distances[j][(int)(intptr_t)centroids[i].center-j-1];
		temp[j] = 0;
		for (j=((int)(intptr_t)centroids[i].center + 1); j < N; j++) 
			temp[j] = distances[(int)(intptr_t)centroids[i].center][j-(int)(intptr_t)centroids[i].center-1];
	}
	return centroids;
}

centroid *matrix_init_concentrate(int **distances, pinfo info, int N) {
	int i, j, t, denominator, *temp;
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
	centroid *centroids = malloc((info->k)*sizeof(centroid));
	for (i=0; i < info->k; i++) {
		centroids[i].info = malloc(N*sizeof(int));
		centroids[i].center = (void *)(intptr_t)v[i].index;
	}
	/**For each centroid**/
	for (i=0; i < info->k; i++) {
		temp = (int *)centroids[i].info;
		/**Create info line with distances**/
		for (j=0; j < (int)(intptr_t)centroids[i].center; j++) 
			temp[j] = distances[j][(int)(intptr_t)centroids[i].center-j-1];
		temp[j] = 0;
		for (j=((int)(intptr_t)centroids[i].center + 1); j < N; j++) 
			temp[j] = distances[(int)(intptr_t)centroids[i].center][j-(int)(intptr_t)centroids[i].center-1];
	}	
	free(v);
	return centroids;
}
