#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "assignment.h"
#define ITEM_ID 15

hash_table *matrix_insert_hash(hash_table *htable, ghashp *g, int **distances, int L, int num_of_hash, int N) {
	int i, j, pos;
	char itemID[ITEM_ID];
	/**For each hash table**/
	for(i = 0; i < L; i++) { 
		/**For each item**/
		for(j = 0; j < N; j++)	{
			sprintf(itemID,"item%d",j+1);
			pos = hash_func_Matrix(g[i],j,distances,num_of_hash,N);
			//pos = mod(pos,htable[i].size); 	IF TABLESIZE = POINTS/8
			insert_chain(itemID,NULL,&(htable[i].table[pos]),3,0,0);
		}
	}
	return htable;
}

pcluster matrix_simplest_assignment(pcluster clusters, int **distances, hash_table htable, int *centroids, int k) {
	int i, j, distance, mindistance, mincentroid, id;
	chainp temp;
	/**For each Bucket**/
	for (i=0; i < htable.size; i++) {
		temp = htable.table[i];
		/**For each Item**/
		while (temp != NULL) {
			id = make_item(temp->key) - 1;
			/**For each Centroid**/
			if (id < centroids[0])  mindistance = distances[id][centroids[0]-id-1];
			else if (id > centroids[0])  mindistance = distances[centroids[0]][id-centroids[0]-1];
			else {
				/**Exclude centroids**/
				temp = temp->next;
				continue;
			}
			mincentroid = 0;
			for (j=1; j < k; j++) {
				if (id < centroids[j])  distance = distances[id][centroids[j]-id-1];
				else if (id > centroids[j])  distance = distances[centroids[j]][id-centroids[j]-1];
				else {
					/**Exclude centroids**/
					mindistance = -1;
					break;
				}
				if (distance < mindistance) {
					 mindistance = distance;
					 mincentroid = j;
				}
			}
			/**Exclude centroids**/
			if (mindistance != -1)	insert_chain(temp->key,NULL,&(clusters[mincentroid].items),3,0,0);
			temp = temp->next;
		}
		
	}
	for (i=0; i < k; i++) 	clusters[i].centroid = (void *)(intptr_t)centroids[i];
	return clusters;
}

pcluster matrix_reverse_approach(pcluster clusters, int **distances, hash_table *htable, int *centroids, int k) {
	int radii;
	radii = compute_start_radius(distances,centroids,k);
	printf("radius: %d\n",radii);
	return clusters;
}

int compute_start_radius(int **distances, int *centroids, int k) {
	int i, j, distance1, distance2, mindistance;
	/**Suppose that at least two centroids exist! First, mindistance is distance of first two centroids 0,1**/
	if (centroids[0] < centroids[1])  mindistance = distances[centroids[0]][centroids[1]-centroids[0]-1];
	else if (centroids[0] > centroids[1])  mindistance = distances[centroids[1]][centroids[0]-centroids[1]-1];
	for (i=0; i < (k - 1); i++) {
		if (centroids[i+1] < centroids[i])  distance1 = distances[centroids[i+1]][centroids[i]-centroids[i+1]-1];
		else if (centroids[i+1] > centroids[i])  distance1 = distances[centroids[i]][centroids[i+1]-centroids[i]-1];
		if (distance1 < mindistance)   mindistance = distance1;
		for (j=(i + 2); j < k; j++) {
			if (centroids[i] < centroids[j])  	distance2 = distances[centroids[i]][centroids[j]-centroids[i]-1];
			else if (centroids[i] > centroids[j])  distance2 = distances[centroids[j]][centroids[i]-centroids[j]-1];
			if (distance2 < mindistance)  mindistance = distance2;
		}
	}
	return (mindistance / 2);
}

int compute_objective_function(pcluster clusters, int **distances, int k) {
	int J = 0, i, id, centroid;
	chainp temp;
	for (i=0; i < k; i++) {
		temp = clusters[i].items;
		centroid = (int)(intptr_t)clusters[i].centroid;
		while (temp != NULL) {
			id = make_item(temp->key) - 1;
			if (id < centroid)  J += distances[id][centroid-id-1];
			else if ( id > centroid)  J += distances[centroid][id-centroid-1];
			temp = temp->next;
		}
	}
	return J;
}
