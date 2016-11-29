#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "medoids.h"
#define ITEM_ID 15

void matrix_medoid(FILE *fp, pinfo info, int ini, int assi, int upd, int clara) {
	char itemsline[7], *allitems, itemID[ITEM_ID];
	int numofitems, token, itemid, tableSize, i, j, pos, flag1, z;
	double totalS, total_t;
	hash_table *htable;
	pcluster clusters;
	clock_t start_t, end_t;
	/**Allocate memory for tables and g functions**/
	htable = malloc(info->L * sizeof(hash_table));
	ghashp *g = malloc(info->L * sizeof(ghashp));		
	for(i = 0; i < info->L; i++)  g[i] = malloc(info->num_of_hash * sizeof(ghash));		
	fscanf(fp,"%s",itemsline);
	allitems = inputString(fp,MAX_LINE);
	/**Read line with items to get the size of the matrix (numofitems x numofitems)**/
	i = 0;
	numofitems = 1;
	while( allitems[i] != '\0')	{
		if( allitems[i] == ',')   numofitems++;
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
	while(fscanf(fp,"%d",&token) != EOF) {
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
	/**Init hash tables and lsh functions**/
	for (i=0; i < info->L; i++)	 init_table(info->num_of_hash,&htable[i],tableSize);
	/**Initialize hash**/
	init_hash_matrix(g,p,info->L,info->num_of_hash,numofitems);
	/**Insert data into hash tables**/
	htable = matrix_insert_hash(htable,g,p,info);
	printf("Insertion completed with success\n");
	start_t = clock();
	if (upd == 1)		clusters = IAU1(htable,g,info,p,ini,assi,3);
	else if (upd == 2)  clusters = IAU2(htable,g,info,p,ini,assi,3);
	/**else {
		
	}**/
	end_t = clock();
	total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
	int totallen = 0;
	for (z=0; z < info->k; z++) {
		printf("\nCluster %d :",(int)(intptr_t)clusters[z].center.center);
		print_points(clusters[z].items);
		totallen += chain_length(clusters[z].items);
		printf("\n");
	}
	printf("total length = %d\n",totallen);
	totalS = compute_silhouette(clusters,p,info,3);
	printf("Time_clustering : %.6lf\n",total_t);
	printf("Sum : %.6lf\n",totalS);
	/**Free memory**/
	for (i=0; i < info->L; i++) free(g[i]);
	free(g);
	for (i=0; i < info->L; i++) destroy_table(&htable[i],3);	
	free(htable);		
	free(allitems);
	for(i=0; i < (info->N - 1); i++)	free(p[i]);	
	free(p);
	for (i=0; i < info->k; i++) destroy_points(&(clusters[i].items));
	free(clusters);
}
