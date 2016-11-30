#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "medoids.h"
#define ITEM_ID 15

void vector_medoid(FILE *fp, FILE *fe, pinfo info, int ini, int assi, int upd, int clara, int flag) {
	int i, j, z, l, col = 1, lines = -2, c, t, itemid, euclID, ch, tableSize, *clarainsert, n;
	char metric[20], m[10], token[100], *eucldata, *ptr, delim[] = " \t", itemID[ITEM_ID], *pnt;
	double **distances, **fulldistances, totalS, total_t;
	hash_table *htable, fulltable;
	pcluster clusters;
	clock_t start_t, end_t;
	ghashp *g = malloc(info->L * sizeof(ghashp));		
	for(i = 0; i < info->L; i++)  g[i] = malloc(info->num_of_hash * sizeof(ghash));	
	/**Read first line to count the dimensions**/
	fscanf(fp,"%s",itemID);
	eucldata = inputString(fp,MAX_LINE);
	/**Count dimensions**/
	ptr = eucldata;
	strtok(ptr,delim);
	ptr = NULL;
	while ((pnt = strtok(ptr,delim)) != NULL)  {
		if (lines == -2)	col++;
		ptr = NULL; 
	}
	/**Dimensions found**/
	/**Return file to start**/
	fseek(fp,0,SEEK_SET);	
	/**Count lines**/
	while ((ch = fgetc(fp)) != EOF) {
		if(ch == '\n') lines++;
	}
	/**Lines found**/
	info->N = lines;
	info->d = col;
	/**Heuristic choice of tableSize**/
	if (flag == 1) 		tableSize = lines/8 + 1;
	else if (flag == 2)	tableSize = 1 << (info->num_of_hash);	
	/**Init hash functions**/	
	if (flag == 1) 		init_hash_Eucl(g,info->L,info->num_of_hash,info->d);	
	else if (flag == 2) init_hash_Cos(g,info->L,info->num_of_hash,info->d);		
	if (clara == 1) {
		fprintf(fe,"CLARA\n");
		start_t = clock();
		n = 40 + 2*info->k;
		/**Keep all elements in one hash table**/
		init_table(info->num_of_hash,&fulltable,tableSize);
		fulltable = vector_one_hash(fulltable,g,info,fp,-1,NULL,flag);
		fulldistances = create_vector_distance_table(fulltable,info->N,info->d,flag);
		clusters = CLARA(info,tableSize,g,n,(void *)fulldistances,fp,fulltable,flag);
		end_t = clock();
		total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
		for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
		fprintf(fe,"Clustering time : %.6lf\n",total_t);
		totalS = compute_silhouette(clusters,fulldistances,info,flag,fe);
		fprintf(fe,"%.6lf]\n\n",totalS);
		if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
		for (i=0; i < (info->N - 1); i++) free(fulldistances[i]);	
		free(fulldistances);
		destroy_table(&fulltable,0);
	} 
	else {
		/**Memory allocation for hash functions and tables**/
		htable = malloc(info->L * sizeof(hash_table));
		/**Init hash tables**/
		for (i=0; i < info->L; i++)	init_table(info->num_of_hash,&htable[i],tableSize);
		/**Insert to hash tables**/
		htable = vector_insert_hash(htable,g,info,fp,flag);
		distances = create_vector_distance_table(htable[0],info->N,info->d,flag);
		if (upd == 1) {
			fprintf(fe,"Algorithm I%dA%dU%d\n",ini,assi,1);
			start_t = clock();
			clusters = IAU1(htable,g,info,distances,ini,assi,flag);
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
			clusters = IAU2(htable,g,info,distances,ini,assi,flag);
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
			clusters = IAU1(htable,g,info,distances,1,1,flag);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,flag,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**k-medoids++, simplest approach, clarans**/
			printf("------------------\nRunning k-medoids++, simplest approach, clarans...\n------------\n");
			fprintf(fe,"\nAlgorithm I1A1U2\n");
			start_t = clock();
			clusters = IAU2(htable,g,info,distances,1,1,flag);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,flag,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**k-medoids++, reverse assignment, a la Loyd's**/
			printf("------------------\nRunning k-medoids++, reverse assignment, a la Loyd's...\n------------\n");
			fprintf(fe,"\nAlgorithm I1A2U1\n");
			start_t = clock();
			clusters = IAU1(htable,g,info,distances,1,2,flag);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,flag,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**k-medoids++, reverse assignment, clarans**/
			printf("------------------\nRunning k-medoids++, reverse assignment, clarans...\n------------\n");
			fprintf(fe,"\nAlgorithm I1A2U2\n");
			start_t = clock();
			clusters = IAU2(htable,g,info,distances,1,2,flag);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,flag,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**concentrate, simplest approach, a la Loyd's**/
			printf("------------------\nRunning concentrate, simplest approach, a la Loyd's...\n------------\n");
			fprintf(fe,"\nAlgorithm I2A1U1\n");
			start_t = clock();
			clusters = IAU1(htable,g,info,distances,2,1,flag);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,flag,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**concentrate, simplest approach, clarans**/
			printf("------------------\nRunning concentrate, simplest approach, clarans...\n------------\n");
			fprintf(fe,"\nAlgorithm I2A1U2\n");
			start_t = clock();
			clusters = IAU2(htable,g,info,distances,2,1,flag);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,flag,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**concentrate, reverse assignment, a la Loyd's**/
			printf("------------------\nRunning concentrate, reverse assignment, a la Loyd's...\n------------\n");
			fprintf(fe,"\nAlgorithm I2A2U1\n");
			start_t = clock();
			clusters = IAU1(htable,g,info,distances,2,2,flag);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,flag,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**concentrate, reverse assignment, clarans**/
			printf("------------------\nRunning concentrate, reverse assignment, clarans...\n------------\n");
			fprintf(fe,"\nAlgorithm I2A2U2\n");
			start_t = clock();
			clusters = IAU2(htable,g,info,distances,2,2,flag);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,distances,info,flag,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
		}
		for (i=0; i < info->L; i++) destroy_table(&htable[i],flag);	
		free(htable);
		for (i=0; i < (info->N - 1); i++) free(distances[i]);	
		free(distances);
	}
	free(eucldata);
	for (i=0; i < info->L; i++)
	{
		for (j=0; j < info->num_of_hash; j++)	free(g[i][j].v);
		free(g[i]);
	}
	free(g);
	for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
	free(clusters);
}
