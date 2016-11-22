#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrixMedoid.h"
#define ITEM_ID 15

void matrix_medoid(FILE *fp, pinfo info, int ini, int assi, int upd) {
	char itemsline[7], *allitems, itemID[ITEM_ID];
	int numofitems, token, itemid, tableSize, i, j, pos, flag1, z, y, J;
	centroid *centroids;
	hash_table *htable;
	/**Allocate memory for tables and g functions**/
	htable = malloc(info->L * sizeof(hash_table));
	ghashp *g = malloc(info->L * sizeof(ghashp));		
	for(i = 0; i < info->L; i++) 
		g[i] = malloc(info->num_of_hash * sizeof(ghash));		
	fscanf(fp,"%s",itemsline);
	allitems = inputString(fp,MAX_LINE);
	/**Read line with items to get the size of the matrix (numofitems x numofitems)**/
	i = 0;
	numofitems = 1;
	while( allitems[i] != '\0')	{
		if( allitems[i] == ',')   numofitems++;
		i++;	
	}
	tableSize = 1 << (info->num_of_hash);
	//tableSize = numofitems / 8 + 1;
	for (i=0; i < info->L; i++)	
			init_table(info->num_of_hash,&htable[i],tableSize);
	/**Read matrix**/
	int **p = malloc((numofitems-1)*sizeof(int*));	
	j = numofitems-1;
	for(i=0; i < numofitems-1; i++)  {
		p[i] = malloc(j*sizeof(int));
		j--;
	}
	i = 0;
	while(fscanf(fp,"%d",&token) != EOF) {
		flag1 = 0;
		z = 0;
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
	init_hash_matrix(g,p,info->L,info->num_of_hash,numofitems);
	/**Get centroids**/
	/**k-medoids++**/
	centroids = matrix_init_kmedoids(p, info, numofitems);
	/**Park-Jun**/
	/**centroids = matrix_init_concentrate(p, info, numofitems);**/
	/*for(i=0; i < info->k; i++) {
		printf("\ncentroids[%d]=%d\n",i,(int)(intptr_t)centroids[i].center);
		for (j=0; j < numofitems; j++) 	printf("%d\t",((int *)centroids[i].info)[j]);
	}
	printf("\n");*/
	/**Insert data into hash tables**/
	htable = matrix_insert_hash(htable,g,p,info->L,info->num_of_hash,numofitems);
	/*for (i=0; i < L; i++) {
		printf("Table %d:\n",i);
		for (j=0; j < tableSize; j++) {
			printf("Bucket %d:\n",j);
			print_chain(htable[i].table[j]);
			printf("\n");
		}	
	}*/
	/**Allocate memory for clusters**/
	pcluster clusters = malloc((info->k)*sizeof(cluster));
	for (i=0; i < info->k; i++) 
		clusters[i].items = NULL;
	/**Assignment**/
	/**clusters = matrix_simplest_assignment(clusters,p,htable[0],centroids,info->k);**/
	clusters = matrix_reverse_approach(clusters,p,htable,g,centroids,info->k,info->num_of_hash,numofitems,info->L);
	for (i=0; i < info->k; i++) {
		printf("\nCluster %d :",(int)(intptr_t)clusters[i].center.center);
		print_chain(clusters[i].items);
		printf("\n");
	}
	/**J = matrix_compute_objective_function(clusters,p,info->k);	
	printf("J=%d\n",J);**/
	/**Free memory**/
	for (i = 0; i < info->L; i++) 
		free(g[i]);
	free(g);
	for (i=0; i < info->L; i++)
		destroy_table(&htable[i],3);	
	free(htable);		
	free(allitems);
	for(i=0; i < numofitems-1; i++)	
		free(p[i]);	
	free(p);
	for (i=0; i < info->k; i++) {
		destroy_chain(&(clusters[i].items),3);
		free(centroids[i].info);
	}
	free(clusters);
	free(centroids);
}
