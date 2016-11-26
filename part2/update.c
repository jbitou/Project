#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structure.h"

centroid *matrix_update_clarans(pcluster clusters, centroid *centroids, cpair pairs, int **distances, int J, pinfo info) {
	int i, j, z, SDj, J1, ci, m, id1, id2, pos, tdistance, *arr;
	pointp temp;
	/**For each pair**/
	for (i=0; i < info->fraction; i++) {
		SDj = 0;
		id1 = (int)(intptr_t)pairs[i].t.center;
		m = (int)(intptr_t)pairs[i].m.center;
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
		J1 = J + SDj;
		/**If J' < J, then swap centroid with medoid**/
		if (J1 < J) {
			for (j=0; j < info->k; j++) {
				if (m == (int)(intptr_t)centroids[j].center) {
					pos = j; 
					break;
				}
			}
			printf("centroid %d replaced by %d with pos %d\n",m,(int)(intptr_t)pairs[i].t.center,pos);
			centroids[pos].center = pairs[i].t.center;
			arr = (int *)centroids[pos].info;
			for (z=0; z < info->N; z++) arr[z] = ((int *)pairs[i].t.info)[z];
		}
	}
	return centroids;
}

centroid *matrix_update_alaloyds(pcluster clusters, centroid *centroids, int J, int **distances, pinfo info) {
	pointp medoid, temp, delete;
	int i, j, z, s, ci, m, id, id1, id2, tdistance, SDj, J1, *qdata, jump, *arr;
	centroid newcenter;
	/**For each cluster**/
	for (i=0; i < info->k; i++) {
		SDj = 0;
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
		/**Insert medoid's info to newcenter**/
		newcenter.center = (void *)(intptr_t)id1;
		newcenter.info = malloc(info->N*sizeof(int));
		arr = (int *)newcenter.info;
		for (z=0; z < info->N; z++) arr[z] = qdata[z];
		free(qdata);	
		/**Store centroid where we are**/
		m = (int)(intptr_t)clusters[i].center.center;
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
		J1 = J + SDj;
		/**If J' < J, then swap centroid with medoid**/
		if (J1 < J) {
			centroids[i].center = newcenter.center;
			arr = (int *)centroids[i].info;
			for (z=0; z < info->N; z++) arr[z] = ((int *)newcenter.info)[z];

		}
		free(newcenter.info);
	}
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
		id1 = make_item(temp->key) - 1;
		sum = 0;
		curr = items;
		/**For each item calculate sum**/
		while (curr != NULL) {
			id2 = make_item(curr->key) - 1;
			if (id1 < id2)  distance = distances[id1][id2-id1-1];
			else if (id2 < id1) distance = distances[id2][id1-id2-1];
			else 	distance = 0;
			sum += distance;
			curr = curr->next;
		}
		if (sum < min) {
			min = sum;
			medoid = temp;
		}
		temp = temp->next;
	}
	printf("min =%d ",min);
	return medoid;
}

int compare_centroids(centroid *centroids, centroid *previous, int k) {
	int i, diff = 0, id1, id2;
	for (i=0; i < k; i++) {
		id1 = (int)(intptr_t)centroids[i].center;	
		id2 = (int)(intptr_t)previous[i].center;	
		if (id1 != id2)	diff++;
	}
	return diff;
}

cpair select_pairs(centroid *centroids, int **distances, pinfo info) {
	int i, j, z, x, index, item, again, *arr, *qdata;
	cpair pairs;
	pairs = malloc(info->fraction*sizeof(pair));
	/**|Q| pairs**/
	for (i=0; i < info->fraction; i++) {
		/**Pick a uniformly distributed integer x**/
		x = (rand() / (RAND_MAX + 1.0)) * (info->k*info->N);
		index = mod(x,info->k);
		pairs[i].m.info = malloc(info->N*sizeof(int));
		arr = (int *)pairs[i].m.info;
		for (z=0; z < info->N; z++) arr[z] = ((int *)centroids[index].info)[z];
		pairs[i].m.center = centroids[index].center;
		item = x / info->k;
		again = 0;
		for (j=0; j < info->k; j++) {
			if (item == (int)(intptr_t)centroids[j].center) {
				i--;
				again = 1;
				free(pairs[i].m.info);
				break;
			}
		}
		/**If t isn't one of centroids**/
		if (again == 0) {
			if (item != 0)	pairs[i].t.center = (void *)(intptr_t)item - 1;
			else	pairs[i].t.center = (void *)(intptr_t)item;
			qdata = malloc((info->N)*sizeof(int));
			for (z=0; z < item; z++) qdata[z] = distances[z][item-z-1];
			qdata[item] = 0;
			for (z=(item+1); z < info->N; z++) qdata[z] = distances[item][z-item-1];
			pairs[i].t.info = malloc(info->N*sizeof(int));
			arr = (int *)pairs[i].t.info;
			for (z=0; z < info->N; z++) arr[z] = qdata[z];
			free(qdata);
		}
	}
	return pairs;
}
