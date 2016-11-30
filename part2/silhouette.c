#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "silhouette.h"

double compute_silhouette(pcluster clusters, void *dis, pinfo info, int flag, FILE *fe) {
	int i, j, **distances, id1, id2, idnext, s, c;
	double **distancesv, si, a, b, sum, suma, sumb, n, max, average, total = 0.0;
	pointp object, temp, second;
	if ((flag == 3) || (flag == 0))		 distances = (int **)dis;
	else if ((flag == 1) || (flag == 2)) distancesv = (double **)dis;
	fprintf(fe,"Silhouette: [");
	/**For every cluster**/
	for (i=0; i < info->k; i++) {
		object = clusters[i].items;
		c = (int)(intptr_t)clusters[i].center.center;	
		sum = average = 0;
		/**For each object of this cluster**/
		while (object != NULL) {
			temp = clusters[i].items;
			/**For each object i (itemid1) of this cluster**/
			id1 = object->position;
			/**Exclude centroid of cluster**/	
			if (id1 == c) {
				object = object->next;
				continue;
			}
			suma = 0;
			/**For each other object in same cluster (exclude distance from itself)**/
			while (temp != NULL) {
				id2 = temp->position;
				if (id1 == id2) {
					temp = temp->next;
					continue;
				}
				if ((flag == 3) || (flag == 0)) {					
					if (id1 < id2)  suma += distances[id1][id2-id1-1];
					else if (id1 > id2)  suma += distances[id2][id1-id2-1];
				}
				else if ((flag == 1) || (flag == 2))  {
					if (id1 < id2)  suma += distancesv[id1][id2-id1-1];
					else if (id1 > id2)  suma += distancesv[id2][id1-id2-1];
				}
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
				id2 = second->position;
				if ((flag == 3) || (flag == 0)) {					
					if (id1 < id2)  sumb += distances[id1][id2-id1-1];
					else if (id1 > id2)  sumb += distances[id2][id1-id2-1];
				}
				else if ((flag == 1) || (flag == 2))  {
					if (id1 < id2)  sumb += distancesv[id1][id2-id1-1];
					else if (id1 > id2)  sumb += distancesv[id2][id1-id2-1];
				}
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
			sum += si;
			object = object->next;
		}
		/**Average of points in cluster i**/
		n = chain_length(clusters[i].items);
		if (n != 0) average = sum / n;
		else average = 0;
		total += average;
		fprintf(fe,"%.5lf, ",average);
	}
	return total;
}

