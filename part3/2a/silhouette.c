#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "silhouette.h"

double compute_silhouette(pcluster clusters, double **data, int N, int numConform, int k) {
	int i, j, id1, id2, idnext, s, c;
	double si, a, b, sum, suma, sumb, n, max, average, total = 0.0, totalsum = 0.0;
	pointp object, temp, second;
	/**For every cluster**/
	for (i=0; i < k; i++) {
		object = clusters[i].items;
		c = clusters[i].center;	
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
				suma += distanceCRMSD(data,N,id1,id2);
				temp = temp->next;
			}
			n = points_length(clusters[i].items);
			/**Calculate average distance of i to other objects in same cluster**/
			a = suma / n;
			/**Find 2nd best neighbor of this object**/
			idnext = object->second;
			for (j=0; j < k; j++) {
				if (idnext == clusters[j].center) {
					s = j;
					break;
				}
			}
			sumb = 0;
			/**For each object in next best cluster**/
			second = clusters[s].items;
			while (second != NULL) {
				id2 = second->position;
				sumb += distanceCRMSD(data,N,id1,id2);
				second = second->next;
			}
			n = points_length(clusters[s].items);
			/**Calculate average distance of i to other objects in same cluster**/
			b = sumb / n;
			if (b == a)	si = 0;
			else {
				max = a;
				if (b > a)	max = b;
				si = (b - a) / max;
			}
			totalsum += si;
			sum += si;
			object = object->next;
		}
		/**Average of points in cluster i**/
		n = points_length(clusters[i].items);
		if (n != 0) average = sum / n;
		else average = 0;
		total += average;
	}
	totalsum /= numConform;
	return totalsum;
}

