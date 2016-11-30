#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "medoids.h"
#define ITEM_ID 15

void matrix_medoid(FILE *fp, FILE *fe, pinfo info, int ini, int assi, int upd, int clara) {
	char itemsline[7], *allitems, itemID[ITEM_ID];
	int numofitems, token, itemid, tableSize, i, j, pos, flag1, z, **fulldistances, n;
	double totalS, total_t;
	hash_table *htable, fulltable;
	pcluster clusters;
	clock_t start_t, end_t;
	ghashp *g = malloc(info->L * sizeof(ghashp));		
	for(i = 0; i < info->L; i++)  g[i] = malloc(info->num_of_hash * sizeof(ghash));		
	fscanf(fp,"%s",itemsline);
	allitems = inputString(fp,MAX_LINE);
	/**Read line with items to get the size of the matrix (numofitems x numofitems)**/
	i = 0;
	numofitems = 1;
	while (allitems[i] != '\0')	{
		if (allitems[i] == ',')   numofitems++;
		i++;	
	}
	info->N = numofitems;
	tableSize = 1 << (info->num_of_hash);
	/**Allocate memory for distance matrix**/
	int **p = malloc((numofitems-1)*sizeof(int *));	
	j = numofitems-1;
	for(i=0; i < numofitems-1; i++)  {
		p[i] = malloc(j*sizeof(int));
		j--;
	}
	/**Read matrix**/
	i = 0;
	while (fscanf(fp,"%d",&token) != EOF) {
		flag1 = z = 0;
		if (token == 0)  flag1 = 1;	
		/**From now on, read next distance and store it**/
		for(j=0; j < numofitems-1; j++) {
			fscanf(fp,"%d",&token);
			if(flag1) {
				p[i][z] = token;
				z++;	
			}	
			if (token == 0)  flag1=1;
		}
		i++;	
	}	
	if (clara == 1) {
		fprintf(fe,"CLARA\n");
		ghashp *fullg =  malloc(sizeof(ghashp));
		fullg[0] =  malloc(info->num_of_hash * sizeof(ghash));
		start_t = clock();
		n = 40 + 2*info->k;
		/**Keep all elements in one hash table**/
		init_table(info->num_of_hash,&fulltable,tableSize);
		init_hash_matrix(fullg,p,1,info->num_of_hash,numofitems);
		fulltable = matrix_one_hash(fulltable,g,p,NULL,info,-1);
		clusters = CLARA(info,tableSize,g,n,(void *)p,fp,fulltable,3);
		end_t = clock();
		total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
		for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
		fprintf(fe,"Clustering time : %.6lf\n",total_t);
		totalS = compute_silhouette(clusters,p,info,0,fe);
		fprintf(fe,"%.6lf]\n\n",totalS);
		if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
		destroy_table(&fulltable,0);
		free(fullg[0]);
		free(fullg);
	} 
	else {
		/**Allocate memory for tables**/
		htable = malloc(info->L * sizeof(hash_table));
		/**Init hash tables and lsh functions**/
		init_hash_matrix(g,p,info->L,info->num_of_hash,numofitems);
		for (i=0; i < info->L; i++)	 init_table(info->num_of_hash,&htable[i],tableSize);
		/**Insert data into hash tables**/
		htable = matrix_insert_hash(htable,g,p,info);
		if (upd == 1) {
			fprintf(fe,"Algorithm I%dA%dU%d\n",ini,assi,1);
			start_t = clock();
			clusters = IAU1(htable,g,info,p,ini,assi,3);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,p,info,3,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
		}
		else if (upd == 2) {
			fprintf(fe,"Algorithm I%dA%dU%d\n",ini,assi,2);
			start_t = clock();
			clusters = IAU2(htable,g,info,p,ini,assi,3);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,p,info,3,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);		
		} 
		else {
			/**k-medoids++, simplest approach, a la Loyd's**/
			printf("------------------\nRunning k-medoids++, simplest approach, a la Loyd's...\n------------\n");
			fprintf(fe,"Algorithm I1A1U1\n");
			start_t = clock();
			clusters = IAU1(htable,g,info,p,1,1,3);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,p,info,3,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**k-medoids++, simplest approach, clarans**/
			printf("------------------\nRunning k-medoids++, simplest approach, clarans...\n------------\n");
			fprintf(fe,"\nAlgorithm I1A1U2\n");
			start_t = clock();
			clusters = IAU2(htable,g,info,p,1,1,3);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,p,info,3,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**k-medoids++, reverse assignment, a la Loyd's**/
			printf("------------------\nRunning k-medoids++, reverse assignment, a la Loyd's...\n------------\n");
			fprintf(fe,"\nAlgorithm I1A2U1\n");
			start_t = clock();
			clusters = IAU1(htable,g,info,p,1,2,3);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,p,info,3,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**k-medoids++, reverse assignment, clarans**/
			printf("------------------\nRunning k-medoids++, reverse assignment, clarans...\n------------\n");
			fprintf(fe,"\nAlgorithm I1A2U2\n");
			start_t = clock();
			clusters = IAU2(htable,g,info,p,1,2,3);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,p,info,3,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**concentrate, simplest approach, a la Loyd's**/
			printf("------------------\nRunning concentrate, simplest approach, a la Loyd's...\n------------\n");
			fprintf(fe,"\nAlgorithm I2A1U1\n");
			start_t = clock();
			clusters = IAU1(htable,g,info,p,2,1,3);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,p,info,3,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**concentrate, simplest approach, clarans**/
			printf("------------------\nRunning concentrate, simplest approach, clarans...\n------------\n");
			fprintf(fe,"\nAlgorithm I2A1U2\n");
			start_t = clock();
			clusters = IAU2(htable,g,info,p,2,1,3);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,p,info,3,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**concentrate, reverse assignment, a la Loyd's**/
			printf("------------------\nRunning concentrate, reverse assignment, a la Loyd's...\n------------\n");
			fprintf(fe,"\nAlgorithm I2A2U1\n");
			start_t = clock();
			clusters = IAU1(htable,g,info,p,2,2,3);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,p,info,3,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
			for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
			free(clusters);
			/**concentrate, reverse assignment, clarans**/
			printf("------------------\nRunning concentrate, reverse assignment, clarans...\n------------\n");
			fprintf(fe,"\nAlgorithm I2A2U2\n");
			start_t = clock();
			clusters = IAU2(htable,g,info,p,2,2,3);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i < info->k; i++) fprintf(fe,"CLUSTER-%d {size: %d, centroid: %s}\n",i+1,chain_length(clusters[i].items),get_centroid_key(clusters[i].items)); 
			fprintf(fe,"Clustering time : %.6lf\n",total_t);
			totalS = compute_silhouette(clusters,p,info,3,fe);
			fprintf(fe,"%.6lf]\n\n",totalS);
			if (info->complete == 1) for (i=0; i < info->k; i++) print_cluster(i+1,clusters[i].items,fe);
		}
		for (i=0; i < info->L; i++) destroy_table(&htable[i],3);	
		free(htable);
	}
	/**Free memory**/
	for (i=0; i < info->L; i++) free(g[i]);
	free(g);	
	free(allitems);
	for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
	free(clusters);
	for(i=0; i < (info->N - 1); i++)	free(p[i]);	
	free(p);
}
