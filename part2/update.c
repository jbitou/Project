#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structure.h"

centroid *matrix_update_alaloyds(pcluster clusters, centroid * centroids, hash_table htable, int J, int **distances, pinfo info) {
	pointp medoid, temp, delete;
	int i, j, z, s, ci, m, id, id1, id2, tdistance, SDj = 0, *qdata, jump;
	centroid newcenter, mcenter;
	/**For each cluster**/
	for (i=0; i < info->k; i++) {
		printf("centroid is %d ",(int)(intptr_t)clusters[i].center.center);
		medoid = matrix_calculate_medoid(clusters[i].items,distances);
		if (medoid == NULL) continue;
		id1 = make_item(medoid->key) - 1;
		printf(" medoid is item%d\n",id1+1);
		/**Create distances of medoid with all the items**/
		qdata = malloc((info->N)*sizeof(int));
		for (z=0; z < id1; z++) qdata[z] = distances[z][id1-z-1];
		qdata[id1] = 0;
		for (z=(id1+1); z < info->N; z++) qdata[z] = distances[id1][z-id1-1];
		newcenter.center = (void *)(intptr_t)id1;
		newcenter.info = (void *) qdata;
		free(qdata);
		mcenter = clusters[i].center;	
		/**Store centroid where we are**/
		m = (int)(intptr_t)clusters[i].center.center;
		/**Swap m and t**/
		clusters[i].center = newcenter;
		/**Check all items using clusters**/
		for (j=0; j < info->k; j++) {
			ci = (int)(intptr_t)clusters[j].center.center;		
			temp = clusters[j].items;
			/**For each item in the cluster j**/
			while (temp != NULL) {
				id2 = make_item(temp->key) - 1;
				if (id1 == id2) {
					temp = temp->next;
					continue;
				}
				if (id1 < id2)  tdistance = distances[id1][id2-id1-1];
				else if (id1 > id2)  tdistance = distances[id2][id1-id2-1];
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
	}
	int totallen = 0;
	for (i=0; i < info->k; i++) {
		printf("\nCluster %d :",(int)(intptr_t)clusters[i].center.center);
		print_points(clusters[i].items);
		totallen += chain_length(clusters[i].items);
		printf("\n");
	}
	printf("total length = %d\n",totallen);
	return centroids;
}


pointp matrix_calculate_medoid(pointp items, int **distances) {
	int min, distance, sum, id1, id2;
	pointp temp, curr, first, medoid;
	temp = first = items;
	if (temp == NULL) return NULL;
	min = 0;
	/**For each item calculate total distance from first item**/
	id1 = make_item(first->key) - 1;
	medoid = temp;
	while (temp != NULL) {
		//printf("for item %s: \n",temp->key);
		id2 = make_item(temp->key) - 1;
		if (id1 < id2)  distance = distances[id1][id2-id1-1];
		else if (id2 < id1) distance = distances[id2][id1-id2-1];
		else 	distance = 0;
		min += distance;
		temp = temp->next;
	}
	/**For each item, beginning from second**/
	temp = first->next;
	while (temp != NULL) {
		//printf("for item %s: \n",temp->key);
		id1 = make_item(temp->key) - 1;
		sum = 0;
		curr = items;
		/**For each item calculate sum**/
		while (curr != NULL) {
			//printf("in loop %s, ",curr->key);
			id2 = make_item(curr->key) - 1;
			if (id1 < id2)  distance = distances[id1][id2-id1-1];
			else if (id2 < id1) distance = distances[id2][id1-id2-1];
			else 	distance = 0;
			sum += distance;
			curr = curr->next;
		}
		//printf("\n");
		if (sum < min) {
			min = sum;
			medoid = temp;
		}
		temp = temp->next;
	}
	printf("min =%d ",min);
	return medoid;
}
