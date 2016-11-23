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
			insert_chain(itemID,NULL,&(htable[i].table[pos]),0,3,0,0);
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
			if (mindistance != -1)	insert_chain(temp->key,NULL,&(clusters[mincentroid].items),0,3,0,0);
			temp = temp->next;
		}
		
	}
	for (i=0; i < k; i++) 	clusters[i].center = centroids[i];
	return clusters;
}

pcluster matrix_reverse_approach(pcluster clusters, int **distances, hash_table *htable, ghashp *g, centroid *centroids, int k, int num_of_hash, int N, int L) {
	int i, j, radii, done, all, previous, pos;
	chainp **barriers, temp;
	/**L tables of pointers**/
	barriers = malloc(L*sizeof(chainp *));
	for (i=0; i < L; i++)	{
		barriers[i] = malloc(pow(2,num_of_hash)*sizeof(chainp));
		for (j=0; j < pow(2,num_of_hash); j++) barriers[i][j] = NULL;
	}
	radii = matrix_compute_start_radius(distances,centroids,k);
	done = all = 0;
	/**Range Search**/
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
	clusters = matrix_remove_clusters_duplicates(clusters,k);
	/**Assign unassigned items**/
	clusters = matrix_assign_rest(clusters,distances,htable,g,barriers,k,num_of_hash,N,L);
	for (i=0; i < L; i++)	free(barriers[i]);
	free(barriers);
	return clusters;
}

pcluster matrix_remove_clusters_duplicates(pcluster clusters, int k) {
	int i, j, jump;
	chainp temp1, temp2;
	/**For each cluster**/
	for (i=0; i < k; i++) {
		temp1 = clusters[i].items;
		/**For each item inside the cluster**/
		while (temp1 != NULL) {
			int jump = 0;
			/**For each item inside every other cluster**/
			for (j=i+1; j < k; j++) {
				temp2 = clusters[j].items;
				while (temp2 != NULL) {
					if (strcmp(temp1->key,temp2->key) == 0) {
						char *tkey = temp2->key;
						if (temp1->distance < temp2->distance)	delete_from_chain(&(clusters[j].items),tkey);
						else  {	
							temp1 = temp1->next;
							delete_from_chain(&(clusters[i].items),tkey);
							jump = 1;
						}
						break;	
					}
					temp2 = temp2->next;
				}
				if (jump == 1) break;
			}
			if (jump != 1)	temp1 = temp1->next;
		}
	}
	return clusters;
}

pcluster matrix_assign_rest(pcluster clusters, int **distances, hash_table *htable, ghashp *g, chainp **barriers, int k, int num_of_hash, int N, int L) {
	int i, j, z, id, pos, assigned, *qdata, distance, mincentroid, center;
	double mindistance;
	chainp temp, check;
	for (i=0; i < htable->size; i++) {
		temp = htable[0].table[i];
		/**For each item**/
		while (temp != NULL) {
			/**Stop when the barrier is found (All assigned points are placed after the barrier)**/
			if ((barriers[0][i] != NULL) && (strcmp(temp->key,barriers[0][i]->key) == 0))  	break;
			id = make_item(temp->key) - 1;
			/**Create info line with distances**/
			qdata = malloc(N*sizeof(int));
			for (z=0; z < id; z++) qdata[z] = distances[z][id-z-1];
			qdata[id] = 0;
			for (z=(id+1); z < N; z++) qdata[z] = distances[id][z-id-1];
			assigned = 0;
			/**Hash in all L tables to check if item is assigned**/
			for (j=1; j < L; j++) {
				pos = hash_func_MSearch(g[j],qdata,distances,num_of_hash,N);
				check = barriers[j][pos];
				while (check != NULL) {
					if (strcmp(check->key,temp->key) == 0) {
						assigned = 1;
						break;
					}
					check = check->next;
				}
				if (assigned == 1) break;
			}
			free(qdata);
			if (assigned == 1) {
				temp = temp->next;
				continue;
			}
			/**For each centroid**/
			center = (int)(intptr_t)clusters[0].center.center;
			if (id < center)  mindistance = distances[id][center-id-1];
			else if (id > center)  mindistance = distances[center][id-center-1];
			mincentroid = 0;
			for (j=1; j < k; j++) {
				center = (int)(intptr_t)clusters[j].center.center;
				if (id < center)  distance = distances[id][center-id-1];
				else if (id > center)  distance = distances[center][id-center-1];
				if (distance < mindistance) {
					 mindistance = distance;
					 mincentroid = j;
				}
			}
			insert_chain(temp->key,NULL,&(clusters[mincentroid].items),mindistance,3,0,0);
			temp = temp->next;
		}
	}
	return clusters;
}

int matrix_compute_start_radius(int **distances, centroid *centroids, int k) {
	int i, j, distance1, distance2, mindistance;
	/**Suppose that at least two centroids/clusters exist! First, mindistance is distance of first two centroids [0],[1]**/
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
	/**For each cluster**/
	for (i=0; i < k; i++) {
		temp = clusters[i].items;
		center = clusters[i].center;
		/**For each item in cluster**/
		while (temp != NULL) {
			id = make_item(temp->key) - 1;
			if (id < (int)(intptr_t)center.center)  J += distances[id][(int)(intptr_t)center.center-id-1];
			else if (id > (int)(intptr_t)center.center)  J += distances[(int)(intptr_t)center.center][id-(int)(intptr_t)center.center-1];
			temp = temp->next;
		}
	}
	return J;
}
