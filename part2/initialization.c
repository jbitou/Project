#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "initialization.h"

int binarySearch(int N, int search, int *array) {
	int first = 0;
	int last = N - 1;
	int middle = (first+last)/2;
	while (first <= last) {
		if (array[middle] < search)
			first = middle + 1;    
		else if (array[middle-1] < search && search <= array[middle]) {
			printf("%d found at location %d.\n", search, middle);
			break;
		}
		else
			last = middle - 1;
		middle = (first + last)/2;
	}
	if (first > last)
		printf("Not found! %d is not present in the list.\n", search);
	return middle; 
}

int *matrix_init_kmedoids(int **distances, pinfo info, int N) {
	int k = 1, i, j, min, max, x;
	int *centroids = malloc((info->k)*sizeof(int));
	srand(time(NULL));
	centroids[0] = (rand() / (RAND_MAX + 1.0)) * N;
	while (k < info->k) {
		int *D = malloc((N-k)*sizeof(int));
		int *P = malloc((N-k)*sizeof(int));
		for (i=0; i < N-k; i++) {
			if (i < centroids[0])	D[i] = distances[i][centroids[0]-i-1];
			for (j=1; j < k; j++) {
				if (i < centroids[j])	min = distances[i][centroids[j]-i-1];
				else if (i > centroids[j]) min = distances[centroids[j]][i-centroids[j]-1];
				if (min < D[i])	D[i] = min;
			}
		}
		max = D[0];
		for (i=1; i < N-k; i++) {
			if (D[i] > max) D[i] = max;
		}
		for (i=0; i < N-k; i++) {
			P[i] = 0;
			for (j=0; j < i; j++) {
				P[i] += pow(D[j],2);
			}
			P[i] /= max;
		}
		x = (rand() / (RAND_MAX + 1.0)) * (P[N-k-1]+1);
		centroids[k] = binarySearch(N, x, P); //r
		free(D);
		free(P);
		k++;
	}
	return centroids;
}
