#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "structure.h"

int intbinarySearch(int Nk, int search, int *array) {
	int first = 0, last = Nk, middle;
	middle = (first+last)/2;
	if (search == 0)  return search;
	while (first <= last) {
		if (array[middle] < search) first = middle + 1;    
		else if (array[middle-1] < search && search <= array[middle]) break;
		else last = middle - 1;
		middle = (first + last)/2;
	}
	return middle; 
}

int doublebinarySearch(int Nk, double search, double *array) {
	int first = 0, last = Nk, middle;
	middle = (first+last)/2;
	if (search == 0)  return search;
	while (first <= last) {
		if (array[middle] < search) first = middle + 1;    
		else if (array[middle-1] < search && search <= array[middle]) break;
		else last = middle - 1;
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

centroid *matrix_init_krandom(pinfo info, int **distances, int n) {
	int i, j, id, *qdata, flag;
	centroid *centroids = malloc((info->k)*sizeof(centroid)); 
	for (i=0; i < info->k; i++) {
		flag = 1;
		while (flag == 1) {
			flag = 0;
			centroids[i].center = (void *)(intptr_t)((rand() / (RAND_MAX + 1.0)) * n);
			for (j=0; j < i; j++) {
				if (centroids[i].center == centroids[j].center) flag = 1;
			}
			if (flag == 0)  break;
		}
		centroids[i].info = malloc(info->N*sizeof(int));
		id = (int)(intptr_t)centroids[i].center;
		qdata = (int *)centroids[i].info;
		for (j=0; j < id; j++) qdata[j] = distances[j][id-j-1];
		qdata[id] = 0;
		for (j=(id+1); j < n; j++) qdata[j] = distances[id][j-id-1];
	}
	return centroids;
}

centroid *vector_init_krandom(pinfo info, double **distances, int n) {
	int i, j, id, flag;
	double *qdata;
	centroid *centroids = malloc((info->k)*sizeof(centroid)); 
	for (i=0; i < info->k; i++) {
		flag = 1;
		while (flag == 1) {
			flag = 0;
			centroids[i].center = (void *)(intptr_t)((rand() / (RAND_MAX + 1.0)) * n);
			for (j=0; j < i; j++) {
				if (centroids[i].center == centroids[j].center) flag = 1;
			}
			if (flag == 0)  break;
		}
		centroids[i].info = malloc(info->N*sizeof(double));
		id = (int)(intptr_t)centroids[i].center;
		qdata = (double *)centroids[i].info;
		for (j=0; j < id; j++) qdata[j] = distances[j][id-j-1];
		qdata[id] = 0;
		for (j=(id+1); j < n; j++) qdata[j] = distances[id][j-id-1];
	}
	return centroids;
}

centroid *matrix_init_kmedoids(int **distances, pinfo info, int N) {
	int k = 1, i, j, x, z, min, max, flag, *temp;
	centroid *centroids = malloc((info->k)*sizeof(centroid)); 
	for (i=0; i < info->k; i++) 	centroids[i].info = malloc(N*sizeof(int));
	centroids[0].center = (void *)(intptr_t)((rand() / (RAND_MAX + 1.0)) * N);
	/**Create k centroids**/
	while (k < info->k) {
		/**For each (k < number of centroids) allocate new Distances and Probabilities array**/
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
		/**Compute maximum of D**/
		max = D[0];
		for (i=1; i < N-k; i++) {
			if (D[i] > max)  max = D[i];
		}
		P[0] = 0;
		/**Compute each P as summary of D^2**/
		for (i=1; i < N-k+1; i++) {
			P[i] = 0;
			for (j=0; j < i; j++)	P[i] += D[j]*D[j];
			/**Make probability < 1**/
			P[i] /= max;
		}
		x = (rand() / (RAND_MAX + 1.0)) * (P[N-k]+1);
		/**r is the new centroid center**/
		centroids[k].center = (void *)(intptr_t)intbinarySearch(N - k, x, P);
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

centroid *vector_init_kmedoids(double **distances, pinfo info, int N) {
	int k = 1, i, j, x, z, flag;
	double min, max, *D, *P, *temp;
	centroid *centroids = malloc((info->k)*sizeof(centroid)); 
	for (i=0; i < info->k; i++) 	centroids[i].info = malloc(N*sizeof(double));
	centroids[0].center = (void *)(intptr_t)((rand() / (RAND_MAX + 1.0)) * N);
	/**Create k centroids**/
	while (k < info->k) {
		/**For each (k < number of centroids) allocate new Distances and Probabilities array**/
		D = malloc((N-k)*sizeof(double));
		P = malloc((N-k+1)*sizeof(double));
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
		/**Compute maximum of D**/
		max = D[0];
		for (i=1; i < N-k; i++) {
			if (D[i] > max)  max = D[i];
		}
		P[0] = 0;
		/**Compute each P as summary of D^2**/
		for (i=1; i < N-k+1; i++) {
			P[i] = 0;
			for (j=0; j < i; j++)	P[i] += D[j]*D[j];
			/**Make probability < 1**/
			P[i] /= max;
		}
		x = (rand() / (RAND_MAX + 1.0)) * (P[N-k]+1);
		/**r is the new centroid center**/
		centroids[k].center = (void *)(intptr_t)doublebinarySearch(N - k, x, P);
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
		temp = (double *)centroids[i].info;
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
	int i, j, t, denominator, *temp, *sums;
	pj_info *v = malloc(N*sizeof(pj_info));
	/**Compute distances of each item with every other**/
	sums = malloc(info->N*sizeof(int));
	for (i=0; i < N; i++) {
		sums[i] = 0;
		for (j=0; j < N; j++) {
			if (j < i)  sums[i] += distances[j][i-j-1];
			else if (j > i)  sums[i] += distances[i][j-i-1];
		}
	}
	/**For each object**/
	for (i=0; i < N; i++) {
		v[i].v = 0;
		v[i].index = i;
		/**For outer Σ**/
		for (j=0; j < N; j++) {
			if (i < j) v[i].v += (double) distances[i][j-i-1] / sums[j];
			else if (i > j)  v[i].v += (double) distances[j][i-j-1] / sums[j];
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
	free(sums);
	return centroids;
}


centroid *vector_init_concentrate(double **distances, pinfo info, int N) {
	int i, j, t;
	double  denominator, *temp, *sums;
	pj_info *v = malloc(N*sizeof(pj_info));
	/**Compute distances of each item with every other**/
	sums = malloc(info->N*sizeof(double));
	for (i=0; i < N; i++) {
		sums[i] = 0;
		for (j=0; j < N; j++) {
			if (j < i)  sums[i] += distances[j][i-j-1];
			else if (j > i)  sums[i] += distances[i][j-i-1];
		}
	}
	/**For each object**/
	for (i=0; i < N; i++) {
		v[i].v = 0;
		v[i].index = i;
		/**For outer Σ**/
		for (j=0; j < N; j++) {
			if (i < j) v[i].v += distances[i][j-i-1] / sums[j];
			else if (i > j)  v[i].v +=  distances[j][i-j-1] / sums[j];
		}
	}
	v = sortArray(v, N);
	centroid *centroids = malloc((info->k)*sizeof(centroid));
	for (i=0; i < info->k; i++) {
		centroids[i].info = malloc(N*sizeof(double));
		centroids[i].center = (void *)(intptr_t)v[i].index;
	}
	/**For each centroid**/
	for (i=0; i < info->k; i++) {
		temp = (double *)centroids[i].info;
		/**Create info line with distances**/
		for (j=0; j < (int)(intptr_t)centroids[i].center; j++) 
			temp[j] = distances[j][(int)(intptr_t)centroids[i].center-j-1];
		temp[j] = 0;
		for (j=((int)(intptr_t)centroids[i].center + 1); j < N; j++) 
			temp[j] = distances[(int)(intptr_t)centroids[i].center][j-(int)(intptr_t)centroids[i].center-1];
	}	
	free(v);
	free(sums);
	return centroids;
}
