#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "inputProcessing.h"

int main (int argc, char **argv) {
	int i, input = -1, output = -1, validate = 0;
	FILE *fp, *fe;
	if (command_processing(argc) == -1)	return -1;
	for (i=1; i < (argc - 1); i+=2) {
		if (strcmp(argv[i],"-d") == 0) 			input = i+1;
		else if (strcmp(argv[i],"-o") == 0)		output = i+1;
	}
	if (strcmp(argv[argc-1],"-validate") == 0)	validate = 1;
	if ((input == -1) || (output == -1)) {
		printf("Arguments -d,-c,-o are required. Try again.\n");
		return -1;
	}
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
	fclose(fp);
	fclose(fe);
	return 0;
}
