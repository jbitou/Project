#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structure.h"
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

pcluster matrix_simplest_assignment(pcluster clusters, int **distances, hash_table htable, centroid *centroids, int k) {
	int i, j, distance, mindistance, mincentroid, id;
	chainp temp;
	/**For each bucket**/
	for (i=0; i < htable.size; i++) {
		temp = htable.table[i];
		/**For each item**/
		while (temp != NULL) {
			id = make_item(temp->key) - 1;
			/**For each centroid**/
			if (id < (int)(intptr_t)centroids[0].center)  mindistance = distances[id][(int)(intptr_t)centroids[0].center-id-1];
			else if (id > (int)(intptr_t)centroids[0].center)  mindistance = distances[(int)(intptr_t)centroids[0].center][id-(int)(intptr_t)centroids[0].center-1];
			else {
				/**Exclude centroids**/
				temp = temp->next;
				continue;
			}
			mincentroid = 0;
			for (j=1; j < k; j++) {
				if (id < (int)(intptr_t)centroids[j].center)  distance = distances[id][(int)(intptr_t)centroids[j].center-id-1];
				else if (id > (int)(intptr_t)centroids[j].center)  distance = distances[(int)(intptr_t)centroids[j].center][id-(int)(intptr_t)centroids[j].center-1];
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
	for (i=0; i < k; i++) 	clusters[i].center = centroids[i];
	return clusters;
}

pcluster matrix_reverse_approach(pcluster clusters, int **distances, hash_table *htable, ghashp *g, centroid *centroids, int k, int num_of_hash, int N, int L) {
	int i, j, z, radii, pos, done, all, previous;
	chainp **barriers, temp1, temp2;
	/**L tables of pointers**/
	barriers = malloc(L*sizeof(chainp *));
	for (i=0; i < L; i++) 
		barriers[i] = malloc(pow(2,num_of_hash)*sizeof(chainp));
	radii = matrix_compute_start_radius(distances,centroids,k);
	done = all = 0;
	do {
		/**For each cluster**/
		previous = done;
		for (i=0; i < k; i++) {
			/**For each table**/
			for (j=0; j < L; j++) {
				pos = hash_func_MSearch(g[j],(int *)centroids[i].info,distances,num_of_hash,N);
				done += search_table_NNR(pos,&(htable[j]),(int *)centroids[i].info,radii,&clusters[i].items,barriers[j],3,0,0,&all);			
			}
			clusters[i].center = centroids[i];
		}
		radii *= 2;
	}while ((done - previous > 1) && (all < L));
	/**If an item is in more than one clusters**/
	int distance1,distance2;
	/**For each cluster**/
	for (i=0; i < k; i++) {
		temp1 = clusters[i].items;
		/**For each item inside the cluster**/
		while (temp1 != NULL) {
			int id = make_item(temp1->key);
			/**For each item inside every other cluster**/
			for (j=i+1; j < k; j++) {
				temp2 = clusters[j].items;
				while (temp2 != NULL) {
					if (strcmp(temp1->key,temp2->key) == 0) {
						if (id < (int)(intptr_t)centroids[i].center) 
							distance1 =  distances[id][(int)(intptr_t)centroids[i].center-id-1];
						if (id > (int)(intptr_t)centroids[i].center) 
							distance1 =  distances[(int)(intptr_t)centroids[i].center][id-(int)(intptr_t)centroids[i].center-1];	
						if (id < (int)(intptr_t)centroids[j].center) 
							distance2 =  distances[id][(int)(intptr_t)centroids[j].center-id-1];
						if (id > (int)(intptr_t)centroids[j].center) 
							distance2 =  distances[(int)(intptr_t)centroids[j].center][id-(int)(intptr_t)centroids[j].center-1];
					}
				}
			}
		}
	}
	return clusters;
}

int matrix_compute_start_radius(int **distances, centroid *centroids, int k) {
	int i, j, distance1, distance2, mindistance;
	/**Suppose that at least two centroids exist! First, mindistance is distance of first two centroids 0,1**/
	if (centroids[0].center < centroids[1].center)  
		mindistance = distances[(int)(intptr_t)centroids[0].center][(int)(intptr_t)centroids[1].center-(int)(intptr_t)centroids[0].center-1];
	else if (centroids[0].center > centroids[1].center)  
		mindistance = distances[(int)(intptr_t)centroids[1].center][(int)(intptr_t)centroids[0].center-(int)(intptr_t)centroids[1].center-1];
	for (i=0; i < (k - 1); i++) {
		if (centroids[i+1].center < centroids[i].center)  
			distance1 = distances[(int)(intptr_t)centroids[i+1].center][(int)(intptr_t)centroids[i].center-(int)(intptr_t)centroids[i+1].center-1];
		else if (centroids[i+1].center > centroids[i].center)  
			distance1 = distances[(int)(intptr_t)centroids[i].center][(int)(intptr_t)centroids[i+1].center-(int)(intptr_t)centroids[i].center-1];
		if (distance1 < mindistance)   mindistance = distance1;
		for (j=(i + 2); j < k; j++) {
			if (centroids[i].center < centroids[j].center)  	
				distance2 = distances[(int)(intptr_t)centroids[i].center][(int)(intptr_t)centroids[j].center-(int)(intptr_t)centroids[i].center-1];
			else if (centroids[i].center > centroids[j].center)  
				distance2 = distances[(int)(intptr_t)centroids[j].center][(int)(intptr_t)centroids[i].center-(int)(intptr_t)centroids[j].center-1];
			if (distance2 < mindistance)  mindistance = distance2;
		}
	}
	return (mindistance / 2);
}

int matrix_compute_objective_function(pcluster clusters, int **distances, int k) {
	int J = 0, i, id; 
	centroid center;
	chainp temp;
	for (i=0; i < k; i++) {
		temp = clusters[i].items;
		center = clusters[i].center;
		while (temp != NULL) {
			id = make_item(temp->key) - 1;
			if (id < (int)(intptr_t)center.center)  J += distances[id][(int)(intptr_t)center.center-id-1];
			else if (id > (int)(intptr_t)center.center)  J += distances[(int)(intptr_t)center.center][id-(int)(intptr_t)center.center-1];
			temp = temp->next;
		}
	}
	return J;
}
