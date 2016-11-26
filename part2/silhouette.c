#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "silhouette.h"

void compute_silhouette(pcluster clusters, int **distances, pinfo info) {
	int i, j, id1, id2, idnext, s, c;
	double si, a, b, sum, suma, sumb, n, max, average;
	pointp object, temp, second;
	/**For every cluster**/
	for (i=0; i < info->k; i++) {
		object = clusters[i].items;
		c = (int)(intptr_t)clusters[i].center.center;	
		printf("for cluster %d\n",c);
		sum = average = 0;
		/**For each object of this cluster**/
		while (object != NULL) {
			temp = clusters[i].items;
			/**For each object i (itemid1) of this cluster**/
			id1 = make_item(object->key) - 1;
			/**Exclude centroid of cluster**/	
			if (id1 == c) {
				object = object->next;
				continue;
			}
			suma = 0;
			/**For each other object in same cluster (exclude distance from itself)**/
			while (temp != NULL) {
				id2 = make_item(temp->key) - 1;
				if (id1 == id2) {
					temp = temp->next;
					continue;
				}
				if (id1 < id2)  suma += distances[id1][id2-id1-1];
				else if (id1 > id2)  suma += distances[id2][id1-id2-1];
				temp = temp->next;
			}
			n = chain_length(clusters[i].items);
			/**Calculate average distance of i to other objects in same cluster**/
			a = suma / n;
			/**Find 2nd best neighbor of this object**/
			idnext = (int)(intptr_t)object->second.center;
			for (j=0; j < info->k; j++) {
				if (idnext == (int)(intptr_t)clusters[j].center.center) {
					s = j;
					break;
				}
			}
			sumb = 0;
			/**For each object in next best cluster**/
			second = clusters[s].items;
			while (second != NULL) {
				id2 = make_item(second->key) - 1;
				if (id1 < id2)  sumb += distances[id1][id2-id1-1];
				else if (id1 > id2)  sumb += distances[id2][id1-id2-1];
				second = second->next;
			}
			n = chain_length(clusters[s].items);
			/**Calculate average distance of i to other objects in same cluster**/
			b = sumb / n;
			if (b == a)	si = 0;
			else {
				max = a;
				if (b > a)	max = b;
				si = (b - a) / max;
			}
			//printf("si: %.5lf\n",si);
			sum += si;
			object = object->next;
		}
		/**Average of points in cluster i**/
		n = chain_length(clusters[i].items);
		average = sum / n;
		printf("Average for cluster %.5lf\n",average);
	}
}
