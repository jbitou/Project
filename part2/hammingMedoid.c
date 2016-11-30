#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "medoids.h"
#define ITEM_ID 15

void hamming_medoid(FILE *fp, FILE *fe, pinfo info, int ini, int assi, int upd, int clara) {
	int i, j, z, y, tableSize, **distances, numofitems = 0, n, *clarainsert, exists, **fulldistances;
	double total_t, totalS, J, min;
	char item[ITEM_ID], data[65], ms[14], space[10]; 
	hash_table *htable, fulltable;
	pcluster clusters;
	centroid *centroids;
	clock_t start_t, end_t;
	ghashp *g = malloc(info->L * sizeof(ghashp));		
	for (i=0; i < info->L; i++)  g[i] = malloc(info->num_of_hash * sizeof(ghash));	
	tableSize = 1 << (info->num_of_hash);
	/**Read  datasets' first line to compute 'N'(number of bits)**/
	fscanf(fp,"%s %s[^\n]",item,data);	
	/**Init hash functions**/
	init_hash_Ham(g,info->L,info->num_of_hash,data);
	fseek(fp,0,SEEK_SET);
	fscanf(fp,"%s%s[^\n]",ms,space);
	while (fscanf(fp,"%s %s[^\n]",item,data) != EOF) numofitems++;
	info->N = numofitems;
	if (clara == 1) {
		fprintf(fe,"CLARA\n");
		start_t = clock();
		n = 40 + 2*info->k;
		/**Keep all elements in one hash table**/
		init_table(info->num_of_hash,&fulltable,tableSize);
		fulltable = hamming_one_hash(fulltable,g,fp,info,NULL,-1);
		fulldistances = create_ham_distance_table(fulltable,info->N);
		clusters = CLARA(info,tableSize,g,n,(void *)fulldistances,fp,fulltable,0);
		end_t = clock();
		total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
		for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
		fprintf(fe,"Clustering time : %.6lf\n",total_t);
		totalS = compute_silhouette(clusters,fulldistances,info,0,fe);
		fprintf(fe,"%.6lf]\n\n",totalS);
		if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
		for (i=0; i < (info->N - 1); i++) free(fulldistances[i]);	
		free(fulldistances);
		destroy_table(&fulltable,0);
	} 
	else {
		/**Memory allocation for and tables**/
		htable = malloc(info->L * sizeof(hash_table));
		for (i=0; i < info->L; i++)	init_table(info->num_of_hash,&htable[i],tableSize);
		/**Insert to hash tables**/
		htable = hamming_insert_hash(htable,g,fp,info);
		distances = create_ham_distance_table(htable[0],info->N);
		if (upd == 1) {
			fprintf(fe,"Algorithm I%dA%dU%d\n",ini,assi,1);
			start_t = clock();
			clusters = IAU1(htable,g,info,distances,ini,assi,0);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,0,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
		}
		else if (upd == 2) {
			fprintf(fe,"Algorithm I%dA%dU%d\n",ini,assi,2);
			start_t = clock();
			clusters = IAU2(htable,g,info,distances,ini,assi,0);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,0,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
		}
		else {
			/**k-medoids++, simplest approach, a la Loyd's**/
			printf("------------------\nRunning k-medoids++, simplest approach, a la Loyd's...\n------------\n");
			fprintf(fe,"Algorithm I1A1U1\n");
			start_t = clock();
			clusters = IAU1(htable,g,info,distances,1,1,0);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,0,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**k-medoids++, simplest approach, clarans**/
			printf("------------------\nRunning k-medoids++, simplest approach, clarans...\n------------\n");
			fprintf(fe,"\nAlgorithm I1A1U2\n");
			start_t = clock();
			clusters = IAU2(htable,g,info,distances,1,1,0);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,0,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**k-medoids++, reverse assignment, a la Loyd's**/
			printf("------------------\nRunning k-medoids++, reverse assignment, a la Loyd's...\n------------\n");
			fprintf(fe,"\nAlgorithm I1A2U1\n");
			start_t = clock();
			clusters = IAU1(htable,g,info,distances,1,2,0);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,0,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**k-medoids++, reverse assignment, clarans**/
			printf("------------------\nRunning k-medoids++, reverse assignment, clarans...\n------------\n");
			fprintf(fe,"\nAlgorithm I1A2U2\n");
			start_t = clock();
			clusters = IAU2(htable,g,info,distances,1,2,0);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,0,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**concentrate, simplest approach, a la Loyd's**/
			printf("------------------\nRunning concentrate, simplest approach, a la Loyd's...\n------------\n");
			fprintf(fe,"\nAlgorithm I2A1U1\n");
			start_t = clock();
			clusters = IAU1(htable,g,info,distances,2,1,0);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,0,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**concentrate, simplest approach, clarans**/
			printf("------------------\nRunning concentrate, simplest approach, clarans...\n------------\n");
			fprintf(fe,"\nAlgorithm I2A1U2\n");
			start_t = clock();
			clusters = IAU2(htable,g,info,distances,2,1,0);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,0,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**concentrate, reverse assignment, a la Loyd's**/
			printf("------------------\nRunning concentrate, reverse assignment, a la Loyd's...\n------------\n");
			fprintf(fe,"\nAlgorithm I2A2U1\n");
			start_t = clock();
			clusters = IAU1(htable,g,info,distances,2,2,0);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,0,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**concentrate, reverse assignment, clarans**/
			printf("------------------\nRunning concentrate, reverse assignment, clarans...\n------------\n");
			fprintf(fe,"\nAlgorithm I2A2U2\n");
			start_t = clock();
			clusters = IAU2(htable,g,info,distances,2,2,0);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,0,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
		}
		for (i=0; i < info->L; i++) destroy_table(&htable[i],0);	
		free(htable);
		for(i=0; i < (info->N - 1); i++) free(distances[i]);	
		free(distances);
	}
	for (i=0; i < info->L; i++) free(g[i]);
	free(g);
	for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
	free(clusters);
}
