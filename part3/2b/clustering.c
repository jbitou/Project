#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "inputProcessing.h"
#define UPPER 3
#define LOWER 2

void experiment(double **data, dinfo **alldistances, int numConform, int N, int r, int T, int k, FILE *fe) {
	int i, j, bestk;
	double **vectors, bestS, total_t;
	pcluster bestclusters;
	dinfo *distances;
	clock_t start_t, end_t;
	printf("start distances\n");
	distances = create_distances(alldistances[0],data,numConform,N,r,T);
	printf("start vectors\n");
	vectors = create_vectors(alldistances,distances,numConform,N,r);
	printf("done vectors\n");
	start_t = clock();
	/*bestclusters = clustering(vectors,numConform,r,100);
	bestS = compute_silhouette(bestclusters,vectors,r,numConform,100);*/
	bestclusters = k_clustering(vectors,k,r,numConform,&bestk,&bestS);
	printf("done k-clustering for r=%d and T=%d\n",r,T);
	end_t = clock();
	total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
	if (T == 0)	fprintf(fe,"%d\t-\t%d\t%.6lf\t%.5lf\n",r,bestk,bestS,total_t);
	else 	fprintf(fe,"%d\t%d\t%d\t%.6lf\t%.5lf\n",r,T,bestk,bestS,total_t);
	if ((clean_memory(vectors,distances,bestclusters,numConform,bestk)) == -1) printf("Failed to free memory\n");
	return;
}

pcluster k_clustering(double **vectors, int k, int r, int numConform, int *bestk, double *bestS) {
	int i, j, counter, betterS = -1, smallerS = 0, previousBest = 0;
	pcluster clusters, bestclusters;
	double silhouette, previousS;
	previousS = -1.1;
	/**Check LOWER times for smaller silhouette and UPPER times for larger silhouette**/
	while ((smallerS < LOWER) && (betterS < UPPER)) {
		printf("smallerS = %d and betterS = %d\n",smallerS,betterS);
		if (k > (int)(numConform * 0.25))	break;
		clusters = clustering(vectors,numConform,r,k);
		silhouette = compute_silhouette(clusters,vectors,r,numConform,k);
		if (silhouette > previousS && silhouette > previousBest) {
			if (betterS != -1) {
				for (i=0; i < *bestk; i++) {
					free(bestclusters[i].center.vector);
					destroy_points(&(bestclusters[i].items));
				}
				free(bestclusters);
			}
			betterS++;
			smallerS = 0;
			*bestk = k;
			*bestS = silhouette;
			bestclusters = malloc(k*sizeof(cluster));
			for (i=0; i < k; i++)	{
				bestclusters[i].center.vector = malloc(r*sizeof(double));
				bestclusters[i].items = NULL;	
				bestclusters[i].items = clone_points(clusters[i].items,r);
			}
			previousBest = silhouette;
		}
		else {
			smallerS++;
			betterS = 0;
		}
		/*for (i=0; i < k; i++) {
			printf("Cluster %d: ",i);
			for (j=0; j < r; j++) printf("%lf\t",clusters[i].center.vector[j]);
			printf("\n");
			print_points(clusters[i].items);
		}*/
		for (i=0; i < k; i++) {
			free(clusters[i].center.vector);
			destroy_points(&(clusters[i].items));
		}
		free(clusters);
		previousS = silhouette;
		k += 5;
	}
	return bestclusters;
}

pcluster clustering(double **vectors, int numConform, int r, int k) {
	int i, j, diff, times = 0;
	double J = 0.0;
	pcluster clusters;
	centroid *centroids, *previous;
	/**k-medoids++ initialization**/
	centroids = vector_init_kmeans(vectors,numConform,k,r);
	previous = malloc(k*sizeof(centroid));
	for (i=0; i < k; i++)	previous[i].vector = malloc(r*sizeof(double));
	do {
		if (times > 0) {
			for (i=0; i < k; i++) {
				free(clusters[i].center.vector);
				destroy_points(&(clusters[i].items));
			}
			free(clusters);
		}
		/**Allocate memory for clusters**/
		clusters = malloc(k*sizeof(cluster));
		for (i=0; i < k; i++) {
			clusters[i].center.vector = malloc(r*sizeof(double));
			clusters[i].items = NULL;
		}	
		/**Simplest assignemt**/
		clusters = vector_simplest_assignment(clusters,vectors,centroids,numConform,r,k);
		J = vector_compute_objective_function(clusters,vectors,r,k);
		for (i=0; i < k; i++) {
			for (j=0; j < r; j++)	previous[i].vector[j] = centroids[i].vector[j];	
		}
		/**Lloyd’s update**/
		centroids = vector_update_loyds(clusters,centroids,J,vectors,r,k);
		diff = compare_centroids(centroids,previous,r,k);		
		times++;
	}while (diff > 0);
	for (i=0; i < k; i++) {
		free(centroids[i].vector);
		free(previous[i].vector);
	}
	free(centroids);
	free(previous);
	return clusters;
}


centroid *vector_init_kmeans(double **vectors, int numConform, int totalk, int r) {
	int k = 1, i, j, x, z, flag;
	double min, max, *D, *P, *temp;
	centroid *centroids = malloc(totalk*sizeof(centroid)); 
	centroids[0].center = (rand() / (RAND_MAX + 1.0)) * numConform;
	for (i=0; i < totalk; i++) centroids[i].vector = malloc(r*sizeof(double));
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
		for (j=0; j < r; j++)	centroids[k].vector[j] = vectors[centroids[k].center][j];	
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
		insert_points(&(clusters[mincentroid].items),mindistance,seconddistance,centroids[secondcentroid],i,r);
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

centroid *vector_update_loyds(pcluster clusters, centroid *centroids, double J, double **vectors, int r, int k) {
	int i, j, z, flag;
	double tdistance, SDj, J1, *m, *ci, *mean;
	centroid newcenter;
	pointp temp;
	/**For each cluster**/
	for (i=0; i < k; i++) {
		SDj = 0;
		/**Calculate mean for this cluster**/
		mean = calculate_mean(clusters[i].items,vectors,r);
		/**Insert mean to newcenter**/
		newcenter.vector = malloc(r*sizeof(double));
		for (j=0; j < r; j++)	newcenter.vector[j] = mean[j];	
		m = malloc(r*sizeof(double));
		for (j=0; j < r; j++)	m[j] = clusters[i].center.vector[j];
		/**Check all items using clusters**/
		for (j=0; j < k; j++) {
			flag = 0;
			ci = malloc(r*sizeof(double));	
			/**Check if ci == m**/
			for (z=0; z < r; z++) {
				ci[z] = clusters[j].center.vector[z];	
				if (m[z] != ci[z]) {
					flag = 1;
					break;
				}
			}
			/**For each item in the cluster j**/
			temp = clusters[j].items;
			while (temp != NULL) {
				tdistance = distance_Euclidean(newcenter.vector,vectors[temp->position],r);
				if (!flag)  {
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
			free(ci);
		}
		J1 = J + SDj;
		/**If J' < J, then swap centroid with mean**/
		if (J1 < J)	{
			for (j=0; j < r; j++)	centroids[i].vector[j] = newcenter.vector[j];
		}
		free(m);
		free(mean);
		free(newcenter.vector);
	}
	return centroids;
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


int compare_centroids(centroid *centroids, centroid *previous, int r, int k) {
	int i, j, diff = 0;
	for (i=0; i < k; i++) {
		for (j=0; j < r; j++) {
			if (centroids[i].vector[j] != previous[i].vector[j]) {	
				diff++;
				break;
			}
		}
	}
	return diff;
}

