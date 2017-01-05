#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "NNsearch.h"

int main (int argc, char **argv) {
	int i, j, input, output, validate, numofusers, numofitems, P;
	char p[4];
	FILE *fp, *fe;
	user *users;
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
	P = count_data(fp,&numofusers,&numofitems);
	printf("P = %d users = %d items = %d\n",P,numofusers,numofitems);
	/**Set R'(u,i)**/
	users = create_users(fp,numofusers,numofitems);
	for (i=0; i < numofusers; i++) {
		printf("user %d: ",users[i].userId);
		for (j=0; j < numofitems; j++)	printf("item %d rating %d\n",users[i].ratings[j].itemId,users[i].ratings[j].rate);
	}
	for (i=0; i < numofusers; i++) free(users[i].ratings);
	free(users);
	fclose(fp);
	fclose(fe);
	return 0;
}
