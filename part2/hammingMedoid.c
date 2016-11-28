#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "medoids.h"
#define ITEM_ID 15

void hamming_medoid(FILE *fp, pinfo info, int ini, int assi, int upd) {
	int i, j, z, tableSize, **distances, numofitems = 0;
	char item[ITEM_ID], data[65], ms[14], space[10]; 
	hash_table *htable;
	pcluster clusters;
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
	/**Insert to hash tables**/
	htable = hamming_insert_hash(htable,g,fp,info);
	printf("Insertion completed with success\n");
	fseek(fp,0,SEEK_SET);
	fscanf(fp,"%s%s[^\n]",ms,space);
	while (fscanf(fp,"%s %s[^\n]",item,data) != EOF) { 
		numofitems++;
	}
	info->N = numofitems;
	printf("N = %d\n",info->N);
	/*for (i=0; i < info->L; i++) {
		printf("Table %d:\n",i);
		for (j=0; j < tableSize; j++) {
			printf("Bucket %d:\n",j);
			print_chain(htable[i].table[j]);
			printf("\n");
		}	
	}*/
	distances = create_ham_distance_table(htable[0],info->N);
	printf("Distance matrix created with success\n");
	z = info->N - 1;
	/*for (i=0; i < info->N-1; i++) {
		for (j=0; j < z; j++) printf("distances[%d][%d]=%d\n",i,j,distances[i][j]);
		z--;
	}*/
}
