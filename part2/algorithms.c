#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "algorithms.h"

pcluster CLARA(pinfo info, int tableSize, ghashp *g, int n, void **fulldistances, FILE *fp, hash_table fulltable, int flag) {
	int i, j, z, y, *clarainsert, **distances, **dis1;
	double J, min, **dis2, **distances1;
	hash_table onetable;
	pcluster clusters, finalclusters;
	centroid *centroids;
	if ((flag == 0) || (flag == 3))	dis1 = (int **)fulldistances;
	if ((flag == 1) || (flag == 2))	dis2 = (double **)fulldistances;
	finalclusters = malloc((info->k)*sizeof(cluster));
	for (i=0; i < info->k; i++)	 finalclusters[i].items = NULL;
	min = J = -1.0;
	for (i=0; i < info->iterations; i++) {
		clusters = malloc((info->k)*sizeof(cluster));
		for (j=0; j < info->k; j++)	 clusters[j].items = NULL;
		clarainsert = choose_clara_sample(n,info->N);
		/**Keep n elements in a hash table**/
		init_table(info->num_of_hash,&onetable,tableSize);
		if (flag == 3)	{
			ghashp *smallg =  malloc(sizeof(ghashp));
			smallg[0] =  malloc(info->num_of_hash * sizeof(ghash));
			distances = create_small_distance_matrix(dis1,n,info->N,clarainsert);
			init_hash_matrix(smallg,distances,1,info->num_of_hash,n);
			onetable = matrix_one_hash(onetable,smallg,distances,clarainsert,info,n);
			centroids = PAM(onetable,g,info,distances,n,flag);
			free(smallg[0]);
			free(smallg);
		}
		else if (flag == 0)	{
			onetable = hamming_one_hash(onetable,g,fp,info,clarainsert,n);
			distances = create_ham_distance_table(onetable,n);
			centroids = PAM(onetable,g,info,distances,n,flag);
		}
		else {
			onetable = vector_one_hash(onetable,g,info,fp,n,clarainsert,flag);
			distances1 = create_vector_distance_table(onetable,info->N,info->d,flag);
			centroids = PAM(onetable,g,info,distances1,n,flag);
		}
		if ((flag == 0) || (flag == 3)) {
			clusters = matrix_simplest_assignment(clusters,dis1,fulltable,centroids,info->k);
			J = matrix_compute_objective_function(clusters,dis1,info->k);
		}
		else if ((flag == 1) || (flag == 2)) {
			clusters = vector_simplest_assignment(clusters,dis2,fulltable,centroids,info->k);
			J = vector_compute_objective_function(clusters,dis2,info->k);
		}
		if ((J < min) || (min == -1.0)) {
			min = J;
			for (j=0; j < info->k; j++) {
				/**Empty the final clusters to insert the new**/
				destroy_points(&(finalclusters[j].items));
				finalclusters[j].center = clusters[j].center;
				finalclusters[j].items = clone(clusters[j].items);
			}
		}
		for (j=0; j < info->k; j++) destroy_points(&(clusters[j].items));
		free(clusters);
		for (j=0; j < info->k; j++) free(centroids[j].info);
		free(centroids);
		free(clarainsert);
		if ((flag == 0) || (flag == 3)) {
			for(i=0; i < (n - 1); i++) free(distances[i]);	
			free(distances);
		}
		else if ((flag == 1) || (flag == 2)) {
			for(i=0; i < (n - 1); i++) free(distances1[i]);	
			free(distances1);
		}
	}
	destroy_table(&onetable,0);
	return finalclusters;
}

centroid *PAM(hash_table htable, ghashp * g, pinfo info, void *dis, int n, int flag) {
	int i, j, z, **p, times = 0, *arr;
	double **dp, J = 0.0, prevJ, *arr1;
	centroid *centroids, *previous;
	pcluster clusters;
	if ((flag == 3) || (flag == 0))	{
		p = (int **)dis;
		centroids = matrix_init_krandom(info,p,n);
	}
	else if ((flag == 1) || (flag == 2)) {
		dp = (double **)dis;
		centroids = vector_init_krandom(info,dp,n);
	}	
	do {
		if (times > 0) {
			prevJ = J;
			for (i=0; i < info->k; i++)	free(previous[i].info);
			free(previous);
		}
		/**Allocate memory for clusters**/
		clusters = malloc((info->k)*sizeof(cluster));
		for (i=0; i < info->k; i++)	 clusters[i].items = NULL;
		if ((flag == 3) || (flag == 0))	{
			clusters = matrix_simplest_assignment(clusters,p,htable,centroids,info->k);
			J = matrix_compute_objective_function(clusters,p,info->k);	
			if (times == 0)	prevJ = J;	
			previous = malloc(info->k*sizeof(centroid));
			for (i=0; i < info->k; i++) {
				previous[i].center = centroids[i].center;
				previous[i].info = malloc(n*sizeof(int));
				arr = (int *)previous[i].info;
				for (j=0; j < n; j++) arr[j] = ((int *)centroids[i].info)[j];
			}
			centroids = matrix_pam_update(clusters,centroids,J,p,info,n);		
		}
		else if ((flag == 1) || (flag == 2)) {
			clusters = vector_simplest_assignment(clusters,dp,htable,centroids,info->k);
			J = vector_compute_objective_function(clusters,dp,info->k);	
			if (times == 0)	prevJ = J;	
			previous = malloc(info->k*sizeof(centroid));
			for (i=0; i < info->k; i++) {
				previous[i].center = centroids[i].center;
				previous[i].info = malloc(n*sizeof(double));
				arr1 = (double *)previous[i].info;
				for (j=0; j < n; j++) arr1[j] = ((double *)centroids[i].info)[j];
			}
			centroids = vector_pam_update(clusters,centroids,J,dp,info,n);	
		}
		times++;
		for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
		free(clusters);
	}while ((J < prevJ) || (times == 1));
	for (i=0; i < info->k; i++) free(centroids[i].info);
	free(centroids);
	return previous;
}

pcluster IAU1(hash_table *htable, ghashp * g, pinfo info, void *dis, int ini, int assi, int flag) {
	int i, j, **p, times = 0;
	double **dp, J = 0.0, prevJ;
	centroid *centroids;
	pcluster clusters;
	if ((flag == 3) || (flag == 0))	p = (int **)dis;
	else if ((flag == 1) || (flag == 2))	dp = (double **)dis;
	/**k-medoids++ initialization**/
	if (ini == 1)  {
		if ((flag == 3) || (flag == 0))	centroids = matrix_init_kmedoids(p,info,info->N);
		else if ((flag == 1) || (flag == 2))	centroids = vector_init_kmedoids(dp,info,info->N);
	}
	/**Park-Jun**/
	else {
		if ((flag == 3) || (flag == 0))	centroids = matrix_init_concentrate(p,info,info->N);
		else if ((flag == 1) || (flag == 2))	centroids = vector_init_concentrate(dp,info,info->N);
	}
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
			if ((flag == 3) || (flag == 0))	clusters = matrix_simplest_assignment(clusters,p,htable[0],centroids,info->k);
			else if ((flag == 1) || (flag == 2))	clusters = vector_simplest_assignment(clusters,dp,htable[0],centroids,info->k);
		}
		/**LSH reverse approach**/
		else { 
			if ((flag == 3) || (flag == 0))	clusters = matrix_reverse_approach(clusters,p,htable,g,centroids,info,flag);
			if ((flag == 1) || (flag == 2))	clusters = vector_reverse_approach(clusters,dp,htable,g,centroids,info,flag);
		}
		if ((flag == 3) || (flag == 0))	{
			J = matrix_compute_objective_function(clusters,p,info->k);
			/**Update à la Lloyd’s**/	
			centroids = matrix_update_alaloyds(clusters,centroids,J,p,info);
		}
		else if ((flag == 1) || (flag == 2))	{
			J = vector_compute_objective_function(clusters,dp,info->k);	
			/**Update à la Lloyd’s**/
			centroids = vector_update_alaloyds(clusters,centroids,J,dp,info);
		}
		if (times == 0)	prevJ = J;	
		times++;
	}while ((J < prevJ) || (times == 1));
	for (i=0; i < info->k; i++) free(centroids[i].info);
	free(centroids);
	return clusters;
}

pcluster IAU2(hash_table *htable, ghashp * g, pinfo info, void *dis, int ini, int assi, int flag) {
	int i, j, z, times, *arr, **p;
	double **dp, *arr1, J = 0.0, prevJ, inprevJ;
	pcluster clusters, finalclusters;
	centroid *centroids;
	cpair pairs;
	finalclusters = malloc((info->k)*sizeof(cluster));
	for (j=0; j < info->k; j++)	 finalclusters[j].items = NULL;
	if ((flag == 3) || (flag == 0))			p = (int **)dis;
	else if ((flag == 1) || (flag == 2))	dp = (double **)dis;
	for (i=0; i < info->iterations; i++) {
		/**k-medoids++ initialization**/
		if (ini == 1)  {
			if ((flag == 3) || (flag == 0))	centroids = matrix_init_kmedoids(p,info,info->N);
			else if ((flag == 1) || (flag == 2))	centroids = vector_init_kmedoids(dp,info,info->N);
		}
		/**Park-Jun**/
		else {
			if ((flag == 3) || (flag == 0))	centroids = matrix_init_concentrate(p,info,info->N);
			else if ((flag == 1) || (flag == 2))	centroids = vector_init_concentrate(dp,info,info->N);
		}
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
				if ((flag == 3) || (flag == 0))		clusters = matrix_simplest_assignment(clusters,p,htable[0],centroids,info->k);
				else if ((flag == 1) || (flag == 2))	clusters = vector_simplest_assignment(clusters,dp,htable[0],centroids,info->k);
			}
			/**LSH reverse approach**/
			else { 
				if ((flag == 3) || (flag == 0))	clusters = matrix_reverse_approach(clusters,p,htable,g,centroids,info,flag);
				if ((flag == 1) || (flag == 2))	clusters = vector_reverse_approach(clusters,dp,htable,g,centroids,info,flag);
			}
			if ((flag == 3) || (flag == 0))	{
				J = matrix_compute_objective_function(clusters,p,info->k);	
				pairs = matrix_select_pairs(centroids,p,info);
				centroids = matrix_update_clarans(clusters,centroids,pairs,p,J,info);
			}
			else if ((flag == 1) || (flag == 2)) {
				J = vector_compute_objective_function(clusters,dp,info->k);	
				pairs = vector_select_pairs(centroids,dp,info);
				centroids = vector_update_clarans(clusters,centroids,pairs,dp,J,info);
			}
			if (times == 0)	inprevJ = J;
			for (z=0; z < info->fraction; z++) {
				free(pairs[z].m.info);
				free(pairs[z].t.info);
			}
			free(pairs);
			times ++;
		}while ((J < inprevJ) || (times == 1));
		if ((J < prevJ) || (prevJ == 0)) {
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
	return finalclusters;
}

int *choose_clara_sample(int n, int N) {
	int exists, i, j, *clarainsert;
	clarainsert = malloc(n*sizeof(int));
	for (i=0; i < n; i++) {
		exists = 1;
		while (exists == 1) {
			exists = 0;
			clarainsert[i] = (rand() / (RAND_MAX + 1.0)) * N;
			for (j=0; j < i; j++) {
				if (clarainsert[i] == clarainsert[j]) exists = 1;
			}
			if (exists == 0)  break;
		}				
	}
	return clarainsert;
}
