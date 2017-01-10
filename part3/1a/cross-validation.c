#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "cross-validation.h"

double fcross_validation(user *users, int numofusers, int numofitems, int P, int flag) {
	int i, j, z, y, numofsubs, *random, r;
	double averageMAE = 0.0;
	user **subsets;
	numofsubs = numofusers / f;
	random = malloc(numofusers*sizeof(int));
	for (i=0; i < numofusers; i++)	random[i] = -1;
	subsets = malloc(numofsubs*sizeof(user *));
	z = 0;
	/**Create subsets of input data**/
	for (i=0; i < numofsubs; i++) {
		subsets[i] = malloc(f*sizeof(user));
		for (j=0; j < f; j++) {
			r = (rand() / (RAND_MAX + 1.0)) * numofusers;
			y = 0;
			/**Pick unique pairs of (user,items)**/
			while (random[y] != -1 && random[y] == r)	{
				r = (rand() / (RAND_MAX + 1.0)) * numofusers;
				y++;
			}
			random[z] = r;
			subsets[i][j] = users[r];
			z++;
		}
	}
	/**Set each subset as validation set**/
	for (i=0; i < 1; i++) {
		printf("for Subset %d:\n",i);
		averageMAE += nnlsh_recommendation(users,subsets[i],f,numofusers,numofitems,P,NULL,flag);
		printf("AVERAGE MAE: %lf\n",averageMAE);
	}
	/*printf("z=%d\n",z);
	for (i=0; i < numofsubs; i++) {
		printf("Subset %d:\n",i);
		for (j=0; j < f; j++)	printf("user %d\t",subsets[i][j].userId);
		printf("\n");
	}*/
	for (i=0; i < numofsubs; i++)	free(subsets[i]);
	free(subsets);
	free(random);
	return averageMAE / numofsubs;
}
