#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "medoids.h"
#define ITEM_ID 15

void hamming_medoid(FILE *fp, pinfo info, int ini, int assi, int upd, int clara) {
	int i, j, z, tableSize, **distances, numofitems = 0, n, *clarainsert;
	double total_t, totalS;
	char item[ITEM_ID], data[65], ms[14], space[10]; 
	hash_table *htable, onetable;
	pcluster clusters;
	clock_t start_t, end_t;
	/**Memory allocation for hash functions and tables**/
	htable = malloc(info->L * sizeof(hash_table));
	ghashp *g = malloc(info->L * sizeof(ghashp));		
	for(i = 0; i < info->L; i++)  g[i] = malloc(info->num_of_hash * sizeof(ghash));	
	tableSize = 1 << (info->num_of_hash);
	/**Read  datasets' first line to compute 'N'(number of bits)**/
	fscanf(fp,"%s %s[^\n]",item,data);	
	/**Init hash tables and lsh functions**/
	for (i=0; i < info->L; i++)	init_table(info->num_of_hash,&htable[i],tableSize);
	init_hash_Ham(g,info->L,info->num_of_hash,data);
	fseek(fp,0,SEEK_SET);
	fscanf(fp,"%s%s[^\n]",ms,space);
	while (fscanf(fp,"%s %s[^\n]",item,data) != EOF) { 
		numofitems++;
	}
	info->N = numofitems;
	printf("N = %d\n",info->N);
	if (clara == 1) {
		centroid *centroids;
		n = 40 + 2*info->k;
		for (i=0; i < info->iterations; i++) {
			clarainsert = malloc(n*sizeof(int));
			for (j=0; j < n; j++) clarainsert[j] = (rand() / (RAND_MAX + 1.0)) * info->N;
			init_table(info->num_of_hash,&onetable,tableSize);
			onetable = hamming_one_hash(onetable,g,fp,info,clarainsert,n);
			/*for (i=0; i < info->L; i++) {
				printf("Table %d:\n",i);
				for (j=0; j < tableSize; j++) {
					printf("Bucket %d:\n",j);
					print_chain(htable[i].table[j]);
					printf("\n");
				}	
			}*/
			distances = create_ham_distance_table(onetable,info->N);
			centroids = CLARA(onetable,g,info,distances,0);
			clusters = matrix_simplest_assignment(clusters,distances,onetable,centroids,info->k);
			for (j=0; j < info->k; j++) free(centroids[j].info);
			free(centroids);
		}
	} 
	else {
		/**Insert to hash tables**/
		htable = hamming_insert_hash(htable,g,fp,info);
		printf("Insertion completed with success\n");
		distances = create_ham_distance_table(htable[0],info->N);
		printf("Distance matrix created with success\n");
		start_t = clock();
		if (upd == 1) 			clusters = IAU1(htable,g,info,distances,ini,assi,0);
		else if (upd == 2) 		clusters = IAU2(htable,g,info,distances,ini,assi,0);
		/**else {
			
		}*/
		end_t = clock();
		total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
	}
	int totallen = 0;
	for (z=0; z < info->k; z++) {
		printf("\nCluster %d :",(int)(intptr_t)clusters[z].center.center);
		print_points(clusters[z].items);
		totallen += chain_length(clusters[z].items);
		printf("\n");
	}
	printf("total length = %d\n",totallen);
	printf("Time_clustering : %.6lf\n",total_t);
	totalS = compute_silhouette(clusters,distances,info,0);
	printf("Sum : %.6lf\n",totalS);
	/*z = info->N - 1;
	for (i=0; i < info->N-1; i++) {
		for (j=0; j < z; j++) printf("distances[%d][%d]=%d\n",i,j,distances[i][j]);
		z--;
	}*/
	for (i=0; i < info->L; i++) free(g[i]);
	free(g);
	for (i=0; i < info->L; i++) destroy_table(&htable[i],0);	
	free(htable);
	for(i=0; i < (info->N - 1); i++) free(distances[i]);	
	free(distances);
	for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
	free(clusters);
}
