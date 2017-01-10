#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "cross-validation.h"

int main (int argc, char **argv) {
	int i, j, input, output, validate, numofusers, numofitems, P, *items, flag;
	char p[4];
	FILE *fp, *fe;
	user *users;
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
	/**Sort the array of items**/
	quickSort(items,0,numofitems-1);
	/**Set R'(u,i)**/
	users = create_users(fp,items,numofusers,numofitems);
	fclose(fp);
	/*for (i=0; i < numofusers; i++) {
		printf("user %d: ",users[i].userId);
		for (j=0; j < numofitems; j++)	{
			if (users[i].ratings[j].rate != 0) printf("item[%d] %d rating %lf\n",j,users[i].ratings[j].itemId,users[i].ratings[j].rate);
		}
	}*/
	/**for (i=1; i < 3; i++)**/
	flag = 2;
		nnlsh_recommendation(users,users,20,numofusers,numofitems,P,fe,flag);
		fcross_validation(users,numofusers,numofitems,P,flag);
	/**}**/
	for (i=0; i < numofusers; i++) free(users[i].ratings);
	free(items);
	free(users);
	fclose(fe);
	return 0;
}
