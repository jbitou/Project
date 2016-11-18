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
			insert_chain(itemID,NULL,&(htable[i].table[pos]),3,0,0);
		}
	}
	return htable;
}

pcluster matrix_simplest_assignment(pcluster clusters, int **distances, hash_table htable, int *centroids ,int k) {
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
