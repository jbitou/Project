#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structure.h"
#define ITEM_ID 15

hash_table *matrix_insert_hash(hash_table *htable, ghashp *g, int **distances, pinfo info) {
	int i, j, pos;
	char itemID[ITEM_ID];
	/**For each hash table**/
	for(i = 0; i < info->L; i++) { 
		/**For each item**/
		for(j = 0; j < info->N; j++)	{
			sprintf(itemID,"item%d",j+1);
			pos = hash_func_Matrix(g[i],j,distances,info->num_of_hash,info->N);
			//pos = mod(pos,htable[i].size); 	IF TABLESIZE = POINTS/8
			insert_chain(itemID,NULL,&(htable[i].table[pos]),3,0,0);
		}
	}
	return htable;
}

pcluster matrix_simplest_assignment(pcluster clusters, int **distances, hash_table htable, centroid *centroids, int k) {
	int i, j, z, distance, mindistance, seconddistance, secondcentroid, mincentroid, id, center;
	chainp temp;
	/**For each bucket**/
	for (i=0; i < htable.size; i++) {
		temp = htable.table[i];
		/**For each item**/
		while (temp != NULL) {
			mindistance = 0;
			id = make_item(temp->key) - 1;
			center = (int)(intptr_t)centroids[0].center;
			if (id < center)  mindistance = distances[id][center-id-1];
			else if (id > center)  mindistance = distances[center][id-center-1];
			else {
				/**Exclude centroids**/
				temp = temp->next;
				continue;
			}
			mincentroid = 0;
			secondcentroid = -1;
			for (j=1; j < k; j++) {
				center = (int)(intptr_t)centroids[j].center;
				if (id < center)  distance = distances[id][center-id-1];
				else if (id > center)  distance = distances[center][id-center-1];
				else {
					/**Exclude centroids**/
					mindistance = -1;
					break;
				}
				if (distance < mindistance) {
					seconddistance = mindistance;
					secondcentroid = mincentroid;
					mindistance = distance;
					mincentroid = j;
				}
			}
			/**If previous for loop didn't find second centroid**/
			if ((mindistance != -1) && (secondcentroid == -1)) {
				/**Typically store a second distance**/
				for (j=0; j < k; j++) {
					if (j != mincentroid) {
						center = (int)(intptr_t)centroids[j].center;
						if (id < center)  seconddistance = distances[id][center-id-1];
						else if (id > center)  seconddistance = distances[center][id-center-1];
						secondcentroid = j;
						break;
					}
				}
				/**Find minimum but exclude the minimum that has been found**/
				for (z=0; z < k; z++) {
					if (z != mincentroid) {
						center = (int)(intptr_t)centroids[z].center;
						if (id < center)  distance = distances[id][center-id-1];
						else if (id > center)  distance = distances[center][id-center-1];
						if (distance < seconddistance) {
							seconddistance = distance;
							secondcentroid = z;
						}
					}
				}
			}
			/**Exclude centroids**/
			if (mindistance != -1)	insert_points(&(clusters[mincentroid].items),temp->key,mindistance,seconddistance,centroids[secondcentroid]);
			temp = temp->next;
		}	
	}
	for (i=0; i < k; i++) 	clusters[i].center = centroids[i];
	return clusters;
}

pcluster matrix_reverse_approach(pcluster clusters, int **distances, hash_table *htable, ghashp *g, centroid *centroids, pinfo info) {
	int i, j, radii, done, all, previous, pos;
	chainp **barriers;
	/**L tables of pointers**/
	barriers = malloc((info->L)*sizeof(chainp *));
	for (i=0; i < info->L; i++)	{
		barriers[i] = malloc(pow(2,info->num_of_hash)*sizeof(chainp));
		for (j=0; j < pow(2,info->num_of_hash); j++) barriers[i][j] = NULL;
	}
	radii = matrix_compute_start_radius(distances,centroids,info->k);
	done = all = 0;
	/**Range Search**/
	do {
		/**For each cluster**/
		previous = done;
		for (i=0; i < info->k; i++) {
			/**For each table**/
			for (j=0; j < info->L; j++) {
				pos = hash_func_MSearch(g[j],(int *)centroids[i].info,distances,info->num_of_hash,info->N);
				done += search_table_NNR(pos,&(htable[j]),(int *)centroids[i].info,radii,&(clusters[i].items),barriers[j],3,0,0,&all);			
			}
			clusters[i].center = centroids[i];
		}
		radii *= 2;
	}while ((done - previous > 1) && (all < info->L));
	/**If an item is in more than one clusters**/
	clusters = matrix_remove_clusters_duplicates(clusters,info->k);
	/**Find second best cluster so far for lsh items**/
	clusters = lsh_second_cluster(clusters,distances,info->k);
	/**Assign unassigned items**/
	clusters = matrix_assign_rest(clusters,distances,htable,g,barriers,info);
	for (i=0; i < info->L; i++)	free(barriers[i]);
	free(barriers);
	return clusters;
}

pcluster lsh_second_cluster(pcluster clusters, int **distances, int k) {
	int i, j, distance, mindistance, id, center, mincentroid;
	pointp temp;
	for (i=0; i < k; i++) {
		temp = clusters[i].items;
		while (temp != NULL) {
			id = make_item(temp->key) - 1;
			for (j=0; j < k; j++) {
				if (i != j) {
					center = (int)(intptr_t)clusters[j].center.center;
					if (id < center)  mindistance = distances[id][center-id-1];
					else if (id > center)  mindistance = distances[center][id-center-1];
					mincentroid = j;
					break;
				}
			}
			for (j=0; j < k; j++) {
				if (i == j) continue;
				center = (int)(intptr_t)clusters[j].center.center;
				if (id < center)  distance = distances[id][center-id-1];
				else if (id > center)  distance = distances[center][id-center-1];
				if (distance < mindistance) {
					mindistance = distance;
					mincentroid = j;
				}
			}
			temp->secdistance = mindistance;
			temp->second = clusters[mincentroid].center;
			temp = temp->next;
		}
	}
	return clusters;
}

pcluster matrix_remove_clusters_duplicates(pcluster clusters, int k) {
	int i, j, jump;
	pointp temp, temp1, temp2;
	/**For each cluster**/
	for (i=0; i < k; i++) {
		temp1 = clusters[i].items;
		/**For each item inside the cluster**/
		while (temp1 != NULL) {
			if (temp1->duplicate == 1) {
				temp1 = temp1->next;
				continue;
			}
			int jump = 0;
			/**For each item inside every other cluster**/
			for (j=i+1; j < k; j++) {
				temp2 = clusters[j].items;
				while (temp2 != NULL) {
					if (strcmp(temp1->key,temp2->key) == 0) {
						temp1->secdistance = -1;
						temp2->secdistance = -1;
						char *tkey = temp2->key;
						if (temp1->mindistance < temp2->mindistance) {
							if ((temp1->secdistance > temp2->mindistance) || (temp1->secdistance == -1)) {
								temp1->secdistance = temp2->mindistance;
								temp1->second = clusters[j].center;
							}
							temp2->duplicate = 1;
						}
						else  {	
							if ((temp2->secdistance > temp1->mindistance) || (temp2->secdistance == -1)) {
								temp2->secdistance = temp1->mindistance;
								temp2->second = clusters[i].center;
							}
							temp1->duplicate = 1;
							jump = 1;
						}
						break;	
					}
					temp2 = temp2->next;
				}
				if (jump == 1) break;
			}
			temp1 = temp1->next;
		}
	}
	for (i=0; i < k; i++) {
		temp = clusters[i].items;
		while (temp != NULL) {
			if (temp->duplicate == 1)	delete_from_chain(&(clusters[i].items),temp->key);
			temp = temp->next;
		}
	}
	return clusters;
}

pcluster matrix_assign_rest(pcluster clusters, int **distances, hash_table *htable, ghashp *g, chainp **barriers, pinfo info) {
	int i, j, z, id, pos, assigned, *qdata, distance, mincentroid, secondcentroid, seconddistance, center, mindistance;
	chainp temp, check;
	for (i=0; i < htable->size; i++) {
		temp = htable[0].table[i];
		/**For each item**/
		while (temp != NULL) {
			/**Stop when the barrier is found (All assigned points are placed after the barrier)**/
			if ((barriers[0][i] != NULL) && (strcmp(temp->key,barriers[0][i]->key) == 0))  	break;
			id = make_item(temp->key) - 1;
			/**Create info line with distances**/
			qdata = malloc((info->N)*sizeof(int));
			for (z=0; z < id; z++) qdata[z] = distances[z][id-z-1];
			qdata[id] = 0;
			for (z=(id+1); z < info->N; z++) qdata[z] = distances[id][z-id-1];
			assigned = 0;
			/**Hash in all L tables to check if item is assigned**/
			for (j=1; j < info->L; j++) {
				pos = hash_func_MSearch(g[j],qdata,distances,info->num_of_hash,info->N);
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
			secondcentroid = -1;
			for (j=1; j < info->k; j++) {
				center = (int)(intptr_t)clusters[j].center.center;
				if (id < center)  distance = distances[id][center-id-1];
				else if (id > center)  distance = distances[center][id-center-1];
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
				for (j=0; j < info->k; j++) {
					if (j != mincentroid) {
						center = (int)(intptr_t)clusters[j].center.center;
						if (id < center)  seconddistance = distances[id][center-id-1];
						else if (id > center)  seconddistance = distances[center][id-center-1];
						secondcentroid = j;
						break;
					}
				}
				/**Find minimum but exclude the minimum that has been found**/
				for (z=0; z < info->k; z++) {
					if (z != mincentroid) {
						center = (int)(intptr_t)clusters[z].center.center;
						if (id < center)  distance = distances[id][center-id-1];
						else if (id > center)  distance = distances[center][id-center-1];
						if (distance < seconddistance) {
							seconddistance = distance;
							secondcentroid = z;
						}
					}
				}
			}
			insert_points(&(clusters[mincentroid].items),temp->key,mindistance,seconddistance,clusters[secondcentroid].center);
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
	pointp temp;
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
