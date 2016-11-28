#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "algorithms.h"

pcluster IAU1(hash_table *htable, ghashp * g, pinfo info, void *dis, int ini, int assi, int flag) {
	int i, j, **p, times;
	double **dp, J = 0.0, prevJ;
	centroid *centroids, *previous;
	pcluster clusters;
	if (flag == 3)	p = (int **)dis;
	else if ((flag == 1) || (flag == 2))	dp = (double **)dis;
	/**k-medoids++ initialization**/
	if (ini == 1)  {
		if (flag == 3)	centroids = matrix_init_kmedoids(p,info,info->N);
		else if ((flag == 1) || (flag == 2))	centroids = vector_init_kmedoids(dp,info,info->N);
	}
	/**Park-Jun**/
	else {
		if (flag == 3)	centroids = matrix_init_concentrate(p,info,info->N);
		else if ((flag == 1) || (flag == 2))	centroids = vector_init_concentrate(dp,info,info->N);
	}
	printf("Initialization completed with success\n");
	times = 0;
	do {
		if (times > 0)	prevJ = J;
		if (times > 0) {
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
		}
		/**Allocate memory for clusters**/
		clusters = malloc((info->k)*sizeof(cluster));
		for (i=0; i < info->k; i++)	 clusters[i].items = NULL;
		/**PAM assignment**/
		if (assi == 1) {
			if (flag == 3)	clusters = matrix_simplest_assignment(clusters,p,htable[0],centroids,info->k);
			else if ((flag == 1) || (flag == 2))	clusters = vector_simplest_assignment(clusters,dp,htable[0],centroids,info->k);
		}
		/**LSH reverse approach**/
		else { 
			if (flag == 3)	clusters = matrix_reverse_approach(clusters,p,htable,g,centroids,info);
			if ((flag == 1) || (flag == 2))	clusters = vector_reverse_approach(clusters,dp,htable,g,centroids,info,flag);
		}
		printf("Assignment completed with success\n");
		if (flag == 3)	J = matrix_compute_objective_function(clusters,p,info->k);	
		else if ((flag == 1) || (flag == 2))	J = vector_compute_objective_function(clusters,dp,info->k);	
		if (times == 0)	prevJ = J;
		printf("J=%.5f\n",J);
		/**Update à la Lloyd’s**/
		if (flag == 3)	centroids = matrix_update_alaloyds(clusters,centroids,J,p,info);
		else if ((flag == 1) || (flag == 2))	centroids = vector_update_alaloyds(clusters,centroids,J,dp,info);
		times++;
	}while ((J < prevJ) || (times == 1));
	for (i=0; i < info->k; i++) free(centroids[i].info);
	free(centroids);
	return clusters;
}

pcluster IAU2(hash_table *htable, ghashp * g, pinfo info, void *dis, int ini, int assi, int flag) {
	int i, j, z, times, *arr, **p;
	double **dp, *arr1, J, prevJ, inprevJ;
	pcluster clusters, finalclusters;
	centroid *centroids, *previous;
	cpair pairs;
	finalclusters = malloc((info->k)*sizeof(cluster));
	for (j=0; j < info->k; j++)	 finalclusters[j].items = NULL;
	if (flag == 3)	p = (int **)dis;
	else if ((flag == 1) || (flag == 2))	dp = (double **)dis;
	J = 0.0;
	for (i=0; i < info->iterations; i++) {
		printf("for s=%d\n",i);
		/**k-medoids++ initialization**/
		if (ini == 1)  {
			if (flag == 3)	centroids = matrix_init_kmedoids(p,info,info->N);
			else if ((flag == 1) || (flag == 2))	centroids = vector_init_kmedoids(dp,info,info->N);
		}
		/**Park-Jun**/
		else {
			if (flag == 3)	centroids = matrix_init_concentrate(p,info,info->N);
			else if ((flag == 1) || (flag == 2))	centroids = vector_init_concentrate(dp,info,info->N);
		}
		printf("Initialization completed with success\n");
		prevJ = J;
		times = 0;
		do {
			if (times > 0)	inprevJ = J;
			if (times > 0) {
				for (j=0; j < info->k; j++) destroy_points(&(clusters[j].items));
				free(clusters);
			}
			clusters = malloc((info->k)*sizeof(cluster));
			for (j=0; j < info->k; j++)	 clusters[j].items = NULL;
			/**PAM assignment**/
			if (assi == 1) {
				if (flag == 3)		clusters = matrix_simplest_assignment(clusters,p,htable[0],centroids,info->k);
				else if ((flag == 1) || (flag == 2))	clusters = vector_simplest_assignment(clusters,dp,htable[0],centroids,info->k);
			}
			/**LSH reverse approach**/
			else { 
				if (flag == 3)	clusters = matrix_reverse_approach(clusters,p,htable,g,centroids,info);
				if ((flag == 1) || (flag == 2))	clusters = vector_reverse_approach(clusters,dp,htable,g,centroids,info,flag);
			}
			printf("Assignment completed with success\n");
			if (flag == 3)	J = matrix_compute_objective_function(clusters,p,info->k);	
			else if ((flag == 1) || (flag == 2))	J = vector_compute_objective_function(clusters,dp,info->k);	
			if (times == 0)	inprevJ = J;
			printf("J=%.5f\n",J);
			if (flag == 3) 		pairs = matrix_select_pairs(centroids,p,info);
			else if ((flag == 1) || (flag == 2)) pairs = vector_select_pairs(centroids,dp,info);
			if (flag == 3) 	centroids = matrix_update_clarans(clusters,centroids,pairs,p,J,info);
			else if ((flag == 1) || (flag == 2)) 	centroids = vector_update_clarans(clusters,centroids,pairs,dp,J,info);
			for (z=0; z < info->fraction; z++) {
				free(pairs[z].m.info);
				free(pairs[z].t.info);
			}
			free(pairs);
			times ++;
		}while ((J < inprevJ) || (times == 1));
		if ((J < prevJ) || (prevJ == 0)) {
			printf("this J=%.5f is better\n",J);
			for (z=0; z < info->k; z++) {
				/**Empty the final clusters to insert the new**/
				destroy_points(&(finalclusters[z].items));
				finalclusters[z].center = clusters[z].center;
				finalclusters[z].items = clone(clusters[z].items);
			}
		}
		for (j=0; j < info->k; j++) destroy_points(&(clusters[j].items));
		free(clusters);
		for (j=0; j < info->k; j++) free(centroids[j].info);
		free(centroids);
	}
	printf("dooooooone\n");
	return finalclusters;
}
