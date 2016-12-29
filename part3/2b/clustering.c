#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "clustering.h"

pcluster clustering(double **vectors, int numConform, int r, int k) {
	int i, j, times = 0;
	double J = 0.0, prevJ, *mean;
	pcluster clusters;
	centroid *centroids;
	/**k-medoids++ initialization**/
	centroids = vector_init_kmeans(vectors,numConform,k,r);
	for (i=0; i < k; i++) {
		printf("centroid = %d: ",centroids[i].center);
		for (j=0; j < r; j++) printf("%lf\t",centroids[i].vector[j]);
		printf("\n");
	}
	/*do {
		if (times > 0) {
			prevJ = J;
			for (i=0; i < k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
		}*/
		/**Allocate memory for clusters**/
		clusters = malloc(k*sizeof(cluster));
		for (i=0; i < k; i++) {
			clusters[i].center.vector = malloc(r*sizeof(double));
			clusters[i].items = NULL;
		}
		/**Simplest assignemt**/
		clusters = vector_simplest_assignment(clusters,vectors,centroids,numConform,r,k);
		J = vector_compute_objective_function(clusters,vectors,r,k);
		printf("J = %lf\n",J);	
		mean = calculate_mean(clusters[0].items,vectors,r);
		for (i=0; i < r; i++)	printf("%lf\t",mean[i]);
		mean = calculate_mean(clusters[1].items,vectors,r);
		for (i=0; i < r; i++)	printf("%lf\t",mean[i]);
		printf("\n");
		/**Update à la Lloyd’s**/
		/*centroids = vector_update_alaloyds(clusters,centroids,J,data,N,k);
		if (times == 0)	prevJ = J;	
		times++;
	}while ((J < prevJ) || (times == 1));*/
	free(centroids);
	for (i=0; i < k; i++) {
		printf("Cluster's centroid %d: ",clusters[i].center.center);
		for (j=0; j < r; j++) printf("%lf\t",clusters[i].center.vector[j]);
		printf("\n");
		print_points(clusters[i].items);
	}
	return clusters;
}


centroid *vector_init_kmeans(double **vectors, int numConform, int totalk, int r) {
	int k = 1, i, j, x, z, flag;
	double min, max, *D, *P, *temp;
	centroid *centroids = malloc(totalk*sizeof(centroid)); 
	centroids[0].center = (rand() / (RAND_MAX + 1.0)) * numConform;
	centroids[0].vector = malloc(r*sizeof(double));
	for (i=0; i < r; i++) centroids[0].vector[i] = vectors[centroids[0].center][i];
	/**Create k centroids**/
	while (k < totalk) {
		/**For each (k < number of centroids) allocate new Distances and Probabilities array**/
		D = malloc((numConform-k)*sizeof(double));
		P = malloc((numConform-k+1)*sizeof(double));
		z = 0;
		/**For each item**/
		for (i=0; i < numConform; i++) {
			flag = 0;
			for (j=0; j < k; j++) {
				if (i == centroids[j].center) flag = 1;
			}
			if (flag == 1) continue;
			D[z] = distance_Euclidean(vectors[i],centroids[0].vector,r);
			for (j=1; j < k; j++) {
				min = distance_Euclidean(vectors[i],centroids[j].vector,r);
				if (min < D[z])	D[z] = min;
			}
			z++;
		}
		/**Compute maximum of D**/
		max = D[0];
		for (i=1; i < numConform-k; i++) {
			if (D[i] > max)  max = D[i];
		}
		P[0] = 0;
		/**Compute each P as summary of D^2**/
		for (i=1; i < numConform-k+1; i++) {
			P[i] = 0;
			for (j=0; j < i; j++)	P[i] += D[j]*D[j];
			/**Make probability < 1**/
			P[i] /= max;
		}
		x = (rand() / (RAND_MAX + 1.0)) * (P[numConform-k]+1);
		/**r is the new centroid center**/
		centroids[k].center = doublebinarySearch(numConform - k, x, P);
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
	for (i=0; i < totalk; i++) {
		centroids[i].vector = malloc(r*sizeof(double));
		for (j=0; j < r; j++)	centroids[i].vector[j] = vectors[centroids[i].center][j];
	}
	return centroids;
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

pcluster vector_simplest_assignment(pcluster clusters, double **vectors, centroid *centroids, int numConform, int r, int k) {
	int i, j, z, secondcentroid, mincentroid;
	double distance, mindistance, seconddistance;
	/**For each comform**/
	for (i=0; i < numConform; i++) {
		mindistance = 0;
		mincentroid = 0;
		mindistance = distance_Euclidean(vectors[i],centroids[0].vector,r);
		secondcentroid = -1;			
		for (j=1; j < k; j++) {
			distance = distance_Euclidean(vectors[i],centroids[j].vector,r);
			if (distance < mindistance) {
				seconddistance = mindistance;
				secondcentroid = mincentroid;
				mindistance = distance;
				mincentroid = j;
			}
		}
		/**If previous for loop didn't find second centroid**/
		if (secondcentroid == -1) {
			/**Typically store a second distance**/
			for (j=0; j < k; j++) {
				if (j != mincentroid) {
					seconddistance = distance_Euclidean(vectors[i],centroids[j].vector,r);
					secondcentroid = j;
					break;
				}
			}
			/**Find minimum but exclude the minimum that has been found**/
			for (z=0; z < k; z++) {
				if (z != mincentroid) {
					distance = distance_Euclidean(vectors[i],centroids[z].vector,r);
					if (distance < seconddistance) {
						seconddistance = distance;
						secondcentroid = z;
					}
				}
			}
		}
		insert_points(&(clusters[mincentroid].items),mindistance,seconddistance,centroids[secondcentroid],i);
	}
	for (i=0; i < k; i++) {
		clusters[i].center.center = centroids[i].center;
		for (j=0; j < r; j++)	clusters[i].center.vector[j] = centroids[i].vector[j];
	}
	return clusters;
}

double vector_compute_objective_function(pcluster clusters, double **vectors, int r, int k) {
	int i; 
	double distance, J = 0.0;
	pointp temp;
	for (i=0; i < k; i++) {
		temp = clusters[i].items;
		while (temp != NULL) {
			distance = distance_Euclidean(clusters[i].center.vector,vectors[temp->position],r);
			J += distance*distance;
			temp = temp->next;
		}
	}
	return J;
}

double *calculate_mean(pointp items, double **vectors, int r) {
	int i, length;
	pointp temp;
	double *mean, sum;
	mean = malloc(r*sizeof(double));
	length =  points_length(items);
	for (i=0; i < r; i++) {
		sum = 0.0;
		temp = items;
		while (temp != NULL) {
			sum += vectors[temp->position][i];
			temp = temp->next;
		}
		mean[i] = sum/length;
	}
	return mean;
}

