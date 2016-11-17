#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrixMedoid.h"
#define ITEM_ID 15

void print_chain(chainp );

void matrix_medoid(FILE *fp, pinfo info) {
	char itemsline[7], *allitems, itemID[ITEM_ID];
	int numofitems, token, itemid, tableSize, i, j, pos, *centroids, flag1, z, y;
	hash_table *htable;
	/**Allocate memory for tables and g functions**/
	htable = malloc(info->L * sizeof(hash_table));
	ghashp *g = malloc(info->L * sizeof(ghashp));		
	for(i = 0; i < info->L; i++) 
		g[i] = malloc(info->num_of_hash * sizeof(ghash));
		
	tableSize = 1 << (info->num_of_hash);
	fscanf(fp,"%s",itemsline);
	allitems = inputString(fp,MAX_LINE);
	/**Read line with items to get the size of the matrix (numofitems x numofitems)**/
	i = 0;
	numofitems = 1;
	while( allitems[i] != '\0')	{
		if( allitems[i] == ',')   numofitems++;
		i++;	
	}
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
	for(i=0; i < info->k; i++)
			printf("centroids[%d]=%d\n",i+1,centroids[i]);
	/**Park-Jun**/
	/*centroids = matrix_init_concentrate(p, info, numofitems);
	for(i=0; i < info->k; i++)
			printf("centroids[%d]=%d\n",i+1,centroids[i]);*/	
	/**Input phase**/	 
	for(i = 0; i < info->L; i++) { 
		for(j = 0; j < numofitems; j++)	{
			sprintf(itemID,"item%d",j+1);
			pos = hash_func_Matrix(g[i],j,p,info->num_of_hash,numofitems);
			insert_chain(itemID,NULL,&(htable[i].table[pos]),3,0,0);
		}
	}
	/**End of input phase**/
	/*for (i=0; i < info->L; i++) {
		printf("Table %d:\n",i);
		for (j=0; j < tableSize; j++) {
			printf("Bucket %d:\n",j);
			print_chain(htable[i].table[j]);
			printf("\n");
		}	
	}*/
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
	free(centroids);
}

void print_chain(chainp l) 
{
	while (l != NULL) 
	{
		printf("key: %s,",l->key);
		l = l->next;
	}
}
