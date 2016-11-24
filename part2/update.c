#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structure.h"

centroid *matrix_update_alaloyds(pcluster clusters, centroid * centroids, int J, int **distances, pinfo info) {
	pointp medoid;
	int i;
	/**For each cluster**/
	for (i=0; i < info->k; i++) {
		medoid = matrix_calculate_medoid(clusters[i].items,distances);
		printf("medoid is: %s\n",medoid->key);
	}
	return centroids;
}


pointp matrix_calculate_medoid(pointp items, int **distances) {
	pointp medoid;
	int min, distance, sum, id1, id2;
	pointp temp, curr, first;
	temp = first = items;
	if (temp == NULL) return NULL;
	medoid = malloc(sizeof(pointp));
	min = 0;
	/**For each item calculate total distance from first item**/
	id1 = make_item(first->key) - 1;
	medoid = temp;
	printf("id1: %d\n",id1);
	while (temp != NULL) {
		id2 = make_item(temp->key) - 1;
		if (id1 < id2)  distance = distances[id1][id2-id1-1];
		else if (id2 < id1) distance = distances[id2][id1-id2-1];
		else 	distance = 0;
		printf("distance from %d:%d\n",id2,distance);
		min += distance;
		temp = temp->next;
	}
	printf("start min: %d\n",min);
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
	printf("min: %d\n",min);
	return medoid;
}
