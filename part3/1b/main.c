#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "cross-validation.h"

int main (int argc, char **argv) {
	int i, j, input, output, validate, numofusers, numofitems, P, k, bestk, *items;
	double mae, bestS;
	char p[4];
	FILE *fp, *fe;
	user *users;
	pcluster bestclusters;
	srand(time(NULL));
	input = output = -1;
	validate = numofusers = numofitems = 0;
	if (command_processing(argv,argc,&input,&output,&validate) == -1)	return -1;
	/**Open input file**/
	fp = fopen(argv[input],"r");
	if (fp == NULL) {
		perror("Error opening input_file");
		return -1;
	}
	/**Open output file**/
	fe = fopen(argv[output],"w+");
	if (fe == NULL) {
		perror("Error opening output_file");
		return -1;
	}
	/**Count data and create array of items**/
	items = create_items(fp,&P,&numofusers,&numofitems);
	printf("P = %d users = %d items = %d\n",P,numofusers,numofitems);
	/**Initialize k - number of clusters**/
	if (numofusers <= 10)	k = 2;
	else if (numofusers <= 100) k = 10;
	else k = numofusers * 0.1;
	/**Sort the array of items**/
	quickSort(items,0,numofitems-1);
	/**Set R'(u,i)**/
	users = create_users(fp,items,numofusers,numofitems);
	/*for (i=0; i < numofusers; i++) {
		printf("user %d: ",users[i].userId);
		for (j=0; j < numofitems; j++)	{
			if (users[i].ratings[j].rate != 0) printf("item[%d] %d rating %lf\n",j,users[i].ratings[j].itemId,users[i].ratings[j].rate);
		}
	}*/
	/**Clustering**/
	k = 270;
	bestclusters = k_clustering(users,numofusers,numofitems,k,&bestk,&bestS);
	//nnlsh_recommendation(users,users,20,numofusers,numofitems,P,fe,flag);
	/*mae = fcross_validation(users,numofusers,numofitems,P,flag);
	if (flag == 1)	fprintf(fe,"Euclidean LSH MAE: %lf\n",mae);
	else fprintf(fe,"Cosine LSH MAE: %lf\n",mae);*/
	for (i=0; i < numofusers; i++) free(users[i].ratings);
	free(users);
	free(items);
	/*for (i=0; i < bestk; i++) destroy_points(&(bestclusters[i].items));
	free(bestclusters);*/
	bestclusters = NULL;
	fclose(fp);
	fclose(fe);
	return 0;
}
