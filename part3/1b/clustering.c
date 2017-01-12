#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "silhouette.h"
#define UPPER 10
#define LOWER 5

pcluster k_clustering(user *users, int numofusers, int numofitems, int k, int *bestk, double *bestS) {
	int i, betterS = 0, smallerS = 0;
	double silhouette, previousS;
	pcluster clusters, bestclusters;
	previousS = -1.1;
	*bestS = -1.1;
	*bestk = k;
	/**Check LOWER times for smaller silhouette and UPPER times for larger silhouette**/
	/*while ((smallerS <= LOWER) && (betterS <= UPPER)) {
		printf("smallerS = %d and betterS = %d\n",smallerS,betterS);
		if (k == numConform / 2)	break;*/
		clusters = clustering(users,numofusers,numofitems,k);
		silhouette = compute_silhouette(clusters,users,numofitems,numofusers,k);
		printf("Silhouette = %lf for k = %d\n",silhouette,k);
		/*if (silhouette > previousS)	{
			if (betterS != 0) {
				for (i=0; i < *bestk; i++) destroy_points(&(bestclusters[i].items));
				free(bestclusters);
			}
			betterS++;
			*bestk = k;
			*bestS = silhouette;
			/**Clone clusters to keep the best**/
			/*bestclusters = malloc(k*sizeof(cluster));
			for (i=0; i < k; i++)	{
				bestclusters[i].items = NULL;	
				bestclusters[i].items = clone_points(clusters[i].items);
			}
		}
		else smallerS++;
		for (i=0; i < k; i++) destroy_points(&(clusters[i].items));
		free(clusters);
		previousS = silhouette;
		k++;
	}
	return bestclusters;*/
	return clusters;
}

pcluster clustering(user *users, int numofusers, int numofitems, int k) {
	int i, j, *centroids, *previous, times = 0, diff;
	double J = 0.0;
	pcluster clusters;
	/**Centroids' initialization**/
	centroids = init_krandom(numofusers,k);
	previous = malloc(k*sizeof(int));
	do {
		if (times > 0) {
			for (i=0; i < k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
		}
		/**Allocate memory for clusters**/
		clusters = malloc(k*sizeof(cluster));
		for (i=0; i < k; i++)	 clusters[i].items = NULL;
		/**Simplest assignemt**/
		clusters = simplest_assignment(clusters,users,centroids,numofusers,numofitems,k);
		for (i=0; i < k; i++) {
			printf("Cluster's centroid: user %d\n",users[clusters[i].center].userId);
			print_points(clusters[i].items);
		}
		J = compute_objective_function(clusters,users,numofitems,k);
		printf("J = %lf\n",J);
		for (i=0; i < k; i++)	previous[i] = centroids[i];	
		/**Update à la Lloyd’s**/
		centroids = update_alaloyds(clusters,centroids,J,users,numofitems,k);
		diff = compare_centroids(centroids,previous,k);	
		printf("diff = %d\n",diff);
		times++;
	}while (diff > 0);
	free(centroids);
	/*for (i=0; i < k; i++) {
		printf("Cluster's centroid: user %d\n",users[clusters[i].center].userId);
		print_points(clusters[i].items);
	}*/
	return clusters;
}

int *init_krandom(int numofusers, int totalk) {
	int i, j, flag;
	int *centroids = malloc(totalk*sizeof(int)); 
	/**Create k centroids**/
	for (i=0; i < totalk; i++) {
		flag = 1;
		while (flag == 1) {
			flag = 0;
			centroids[i] = (rand() / (RAND_MAX + 1.0)) * numofusers;
			/**Check if centroid already exists**/
			for (j=0; j < i; j++) {
				if (centroids[i] == centroids[j]) flag = 1;
			}
			if (flag == 0)  break;
		}
	}
	return centroids;
}

pcluster simplest_assignment(pcluster clusters, user *users, int *centroids, int numofusers, int numofitems, int k) {
	int i, j, secondcentroid, mincentroid;
	double distance, tempdis, mindistance, seconddistance;
	/**For each user**/
	for (i=0; i < numofusers; i++) {
		secondcentroid = centroids[1];
		mincentroid = 0;
		if (metric == 1) {
			mindistance = distance_Euclidean(users[i].ratings,users[centroids[0]].ratings,numofitems);	
			seconddistance = distance_Euclidean(users[i].ratings,users[centroids[1]].ratings,numofitems);
		}
		else if (metric == 2) {
			mindistance = distance_Cosine(users[i].ratings,users[centroids[0]].ratings,numofitems);	
			seconddistance = distance_Cosine(users[i].ratings,users[centroids[1]].ratings,numofitems);
		}
		if (mindistance > seconddistance) {
			tempdis = mindistance;
			mindistance = seconddistance;
			seconddistance = tempdis;
			secondcentroid = centroids[0];
			mincentroid = 1;
		}	
		for (j=2; j < k; j++) {
			if (metric == 1)		distance = distance_Euclidean(users[i].ratings,users[centroids[j]].ratings,numofitems);
			else if (metric == 2)	distance = distance_Cosine(users[i].ratings,users[centroids[j]].ratings,numofitems);
			if (distance < mindistance) {
				seconddistance = mindistance;
				secondcentroid = centroids[mincentroid];
				mindistance = distance;
				mincentroid = j;
			}
			else if ((distance > mindistance) && (distance < seconddistance)) {
				seconddistance = distance;
				secondcentroid = centroids[j];
			}
		}
		/**Put each centroid in its cluster**/
		for (j=0; j < k; j++) {
			if (centroids[j] == i) {
				mincentroid = j;
				mindistance = 0.0;
				break;
			}
		}
		insert_points(&(clusters[mincentroid].items),mindistance,seconddistance,secondcentroid,i);
	}
	for (i=0; i < k; i++) 	clusters[i].center = centroids[i];
	return clusters;
}

double compute_objective_function(pcluster clusters, user *users, int numofitems, int k) {
	int i, center; 
	double J = 0.0;
	pointp temp;
	/**For each cluster**/
	for (i=0; i < k; i++) {
		temp = clusters[i].items;
		center = clusters[i].center;
		/**For each item in cluster**/
		while (temp != NULL) {
			if (metric == 1)		J += distance_Euclidean(users[temp->position].ratings,users[center].ratings,numofitems);
			else if (metric == 2)	J += distance_Cosine(users[temp->position].ratings,users[center].ratings,numofitems);
			temp = temp->next;
		}
	}
	return J;
}

int *update_alaloyds(pcluster clusters, int *centroids, double J, user *users, int numofitems, int k) {
	pointp medoid, temp, delete;
	int i, j, z, s, ci, m, id1, id2, newcenter;
	double tdistance, SDj, J1;
	/**For each cluster**/
	for (i=0; i < k; i++) {
		SDj = 0;
		medoid = calculate_medoid(clusters[i].items,users,numofitems);
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
				if (metric == 1)		tdistance = distance_Euclidean(users[id1].ratings,users[id2].ratings,numofitems);
				else if (metric == 2)	tdistance = distance_Cosine(users[id1].ratings,users[id2].ratings,numofitems);
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

pointp calculate_medoid(pointp items, user *users, int numofitems) {
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
		if (metric == 1)		distance = distance_Euclidean(users[id1].ratings,users[id2].ratings,numofitems);
		else if (metric == 2)	distance = distance_Cosine(users[id1].ratings,users[id2].ratings,numofitems);
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
			if (metric == 1)		distance = distance_Euclidean(users[id1].ratings,users[id2].ratings,numofitems);
			else if (metric == 2)	distance = distance_Cosine(users[id1].ratings,users[id2].ratings,numofitems);
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

int compare_centroids(int *centroids, int *previous, int k) {
	int i, diff = 0;
	for (i=0; i < k; i++) {
		if (centroids[i] != previous[i]) diff++;
	}
	return diff;
}
