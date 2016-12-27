#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "clustering.h"
#include "distanceCRMSD.h"

pcluster clustering(double **data, int numConform, int N, int k) {
	printf("start\n");
	int i, j, *centroids, times = 0;
	double J = 0.0, prevJ;
	pcluster clusters;
	/**k-medoids++ initialization**/
	centroids = vector_init_kmedoids(data,numConform,k,N);
	printf("1st centroid = %d , 2nd centroid = %d\n",centroids[0],centroids[1]);
	do {
		if (times > 0) {
			prevJ = J;
			for (i=0; i < k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
		}
		/**Allocate memory for clusters**/
		clusters = malloc(k*sizeof(cluster));
		for (i=0; i < k; i++)	 clusters[i].items = NULL;
		/**Simplest assignemt**/
		clusters = vector_simplest_assignment(clusters,data,centroids,numConform,N,k);
		J = vector_compute_objective_function(clusters,data,N,k);
		printf("J = %lf\n",J);	
		/**Update à la Lloyd’s**/
		centroids = vector_update_alaloyds(clusters,centroids,J,data,N,k);
		if (times == 0)	prevJ = J;	
		times++;
	}while ((J < prevJ) || (times == 1));
	free(centroids);
	for (i=0; i < k; i++) {
		printf("Cluster's centroid: conform %d\n",clusters[i].center);
		print_points(clusters[i].items);
	}
	return clusters;
}

int *vector_init_kmedoids(double **data, int numConform, int totalk, int N) {
	int k = 1, i, j, x, z, flag;
	double min, max, *D, *P, *temp;
	int *centroids = malloc(totalk*sizeof(int)); 
	centroids[0] = (rand() / (RAND_MAX + 1.0)) * numConform;
	/**Create k centroids**/
	while (k < totalk) {
		printf("k = %d\n",k);
		/**For each (k < number of centroids) allocate new Distances and Probabilities array**/
		D = malloc((numConform-k)*sizeof(double));
		P = malloc((numConform-k+1)*sizeof(double));
		z = 0;
		/**For each item**/
		for (i=0; i < numConform; i++) {
			flag = 0;
			for (j=0; j < k; j++) {
				if (i == centroids[j]) flag = 1;
			}
			if (flag == 1) continue;
			D[z] = distanceCRMSD(data,N,i,centroids[0]);
			for (j=1; j < k; j++) {
				min = distanceCRMSD(data,N,i,centroids[j]);
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
		centroids[k] = doublebinarySearch(numConform - k, x, P);
		for (i=0; i < k; i++) {
			if (centroids[k] == centroids[i]) {
				k--;
				break;
			}
		}		
		free(D);
		free(P);
		k++;
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

pcluster vector_simplest_assignment(pcluster clusters, double **data, int *centroids, int numConform, int N, int k) {
	int i, j, z, secondcentroid, mincentroid;
	double distance, mindistance, seconddistance;
	/**For each comform**/
	for (i=0; i < numConform; i++) {
		mindistance = 0;
		mincentroid = 0;
		mindistance = distanceCRMSD(data,N,i,centroids[0]);
		secondcentroid = -1;			
		for (j=1; j < k; j++) {
			distance = distanceCRMSD(data,N,i,centroids[j]);
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
					seconddistance = distanceCRMSD(data,N,i,centroids[j]);
					secondcentroid = j;
					break;
				}
			}
			/**Find minimum but exclude the minimum that has been found**/
			for (z=0; z < k; z++) {
				if (z != mincentroid) {
					distance = distanceCRMSD(data,N,i,centroids[z]);
					if (distance < seconddistance) {
						seconddistance = distance;
						secondcentroid = z;
					}
				}
			}
		}
		insert_points(&(clusters[mincentroid].items),mindistance,seconddistance,centroids[secondcentroid],i);
	}
	for (i=0; i < k; i++) 	clusters[i].center = centroids[i];
	return clusters;
}

double vector_compute_objective_function(pcluster clusters, double **data, int N, int k) {
	int i, center; 
	double J = 0.0;
	pointp temp;
	/**For each cluster**/
	for (i=0; i < k; i++) {
		temp = clusters[i].items;
		center = clusters[i].center;
		/**For each item in cluster**/
		while (temp != NULL) {
			J += distanceCRMSD(data,N,temp->position,center);
			temp = temp->next;
		}
	}
	return J;
}

int *vector_update_alaloyds(pcluster clusters, int *centroids, double J, double **data, int N, int k) {
	pointp medoid, temp, delete;
	int i, j, z, s, ci, m, id1, id2, newcenter;
	double tdistance, SDj, J1;
	/**For each cluster**/
	for (i=0; i < k; i++) {
		SDj = 0;
		medoid = vector_calculate_medoid(clusters[i].items,data,N);
		if (medoid == NULL) continue;
		id1 = medoid->position;
		/**Insert medoid's info to newcenter**/
		newcenter = id1;
		m = clusters[i].center;
		/**Check all items using clusters**/
		for (j=0; j < k; j++) {
			ci = clusters[j].center;		
			temp = clusters[j].items;
			/**For each item in the cluster j**/
			while (temp != NULL) {
				id2 = temp->position;
				if (id1 == id2) {
					temp = temp->next;
					continue;
				}
				tdistance = distanceCRMSD(data,N,id1,id2);
				if (ci == m)  {
					/**If dist(i,t) > dist(i,c')**/
					if (tdistance > temp->secdistance)	SDj += temp->secdistance - temp->mindistance;
					else  SDj += tdistance - temp->mindistance;
				}
				else {
					/**if dist(i,t) < dist(i,c(i))**/
					if (tdistance < temp->mindistance)	SDj += tdistance - temp->mindistance;
				}
				temp = temp->next;
			}
		}
		J1 = J + SDj;
		/**If J' < J, then swap centroid with medoid**/
		if (J1 < J)	centroids[i] = newcenter;
	}
	return centroids;
}

pointp vector_calculate_medoid(pointp items, double **data, int N) {
	int id1, id2;
	double min, distance, sum;
	pointp temp, curr, first, medoid;
	temp = first = items;
	if (temp == NULL) return NULL;
	min = 0;
	/**For each item calculate total distance from first item**/
	id1 = first->position;
	medoid = temp;
	while (temp != NULL) {
		id2 = temp->position;
		distance = distanceCRMSD(data,N,id1,id2);
		min += distance;
		temp = temp->next;
	}
	/**For each item, beginning from second**/
	temp = first->next;
	while (temp != NULL) {
		id1 = temp->position;
		sum = 0;
		curr = items;
		/**For each item calculate sum**/
		while (curr != NULL) {
			id2 = curr->position;
			distance = distanceCRMSD(data,N,id1,id2);
			sum += distance;
			curr = curr->next;
		}
		if (sum < min) {
			min = sum;
			medoid = temp;
		}
		temp = temp->next;
	}
	return medoid;
}
