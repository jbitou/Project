#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "medoids.h"
#define ITEM_ID 15

void vector_medoid(FILE *fp, pinfo info, int ini, int assi, int upd, int flag) {
	int i, j, z, l, col = 1, lines = -2, c, t, itemid, euclID, ch, tableSize;
	char metric[20], m[10], token[100], *eucldata, *ptr, delim[] = " \t", itemID[ITEM_ID], *pnt;
	double **distances, totalS;
	hash_table *htable;
	pcluster clusters;
	/**Memory allocation for hash functions and tables**/
	htable = malloc(info->L * sizeof(hash_table));
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
	printf("lines = %d\n",lines);
	printf("col = %d\n",col);
	/**Heuristic choice of tableSize**/
	if (flag == 1) 		tableSize = lines/8 + 1;
	else if (flag == 2)	tableSize = 1 << (info->num_of_hash);		
	/**Init hash tables and lsh functions**/
	for (i=0; i < info->L; i++)	init_table(info->num_of_hash,&htable[i],tableSize);
	if (flag == 1) 		init_hash_Eucl(g,info->L,info->num_of_hash,info->d);	
	else if (flag == 2) init_hash_Cos(g,info->L,info->num_of_hash,info->d);	
	/**Insert to hash tables**/
	htable = euclidean_insert_hash(htable,g,info,fp,tableSize);
	printf("Insertion completed with success\n");
	distances = create_distance_table(htable[0],info->N,info->d);
	printf("Distance matrix created with success\n");
	if (upd == 1)		clusters = IAU1(htable,g,info,distances,ini,assi,flag);
	else if (upd == 2)  clusters = IAU2(htable,g,info,distances,ini,assi,flag);
	/**else {
		
	}**/
	/*int totallen = 0;
	for (z=0; z < info->k; z++) {
		printf("\nCluster %d :",(int)(intptr_t)clusters[z].center.center);
		print_points(clusters[z].items);
		totallen += chain_length(clusters[z].items);
		printf("\n");
	}
	printf("total length = %d\n",totallen);*/
	totalS = compute_silhouette(clusters,distances,info,flag);
	printf("Sum : %.6lf\n",totalS);
	/*for (i=0; i < info->L; i++) {
		printf("Table %d:\n",i);
		for (j=0; j < tableSize; j++) {
			printf("Bucket %d:\n",j);
			print_chain(htable[i].table[j]);
			printf("\n");
		}	
	}*/
	/*z = info->N - 1;
	for (i=0; i < info->N-1; i++) {
		for (j=0; j < z; j++) printf("distances[%d][%d]=%.10f\n",i,j,distances[i][j]);
		z--;
	}*/
	free(eucldata);
	for (i=0; i < info->L; i++)
	{
		for (j=0; j < info->num_of_hash; j++)	free(g[i][j].v);
		free(g[i]);
	}
	free(g);
	for (i=0; i < info->L; i++) destroy_table(&htable[i],flag);	
	free(htable);
	for(i=0; i < (info->N - 1); i++) free(distances[i]);	
	free(distances);
	for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
	free(clusters);
}
