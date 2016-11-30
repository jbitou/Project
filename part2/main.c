#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "medoids.h"
#define ITEM_ID 15

int main (int argc, char **argv) {
	FILE *fp, *fc, *fe;
	int i, complete = 0, input, config, output, ini = 0, assi = 0, upd = 0, clara = 0, ch, lines = -2, flag;
	char ms[14], space[10], m[10], metric[20];
	pinfo info;
	srand(time(NULL));		
	if (command_processing(argc) == -1)	return -1;
	input = config = output = -1;
	for (i=1; i < (argc - 1); i+=2) {
		if (strcmp(argv[i],"-d") == 0) 			input = i+1;
		else if (strcmp(argv[i],"-c") == 0) 	config = i+1;
		else if (strcmp(argv[i],"-o") == 0)		output = i+1;
	}
	if (strcmp(argv[argc-1],"-complete") == 0)	complete = 1;
	if ((input == -1) || (config == -1) || (output == -1)) {
		printf("Arguments -d,-c,-o are required. Try again.\n");
		return -1;
	}
	/**Open input file**/
	fp = fopen(argv[input],"r");
	if (fp == NULL) {
		perror("Error opening input_file");
		return -1;
	}
	/**Open configuration file**/
	fc = fopen(argv[config],"r");
	if (fc == NULL) {
		perror("Error opening configuration_file");
		return -1;
	}
	/**Open output file**/
	fe = fopen(argv[output],"w+");
	if (fe == NULL) {
		perror("Error opening output_file");
		return -1;
	}
	/**Ask user for combination of methods**/
	user_choice(&ini,&assi,&upd,&clara);
	/**Count lines of input file**/
	while ((ch = fgetc(fp)) != EOF) {
		if (ch == '\n') lines++;
	}
	/**Lines found**/
	fseek(fp,0,SEEK_SET);
	/**Read first line of input_file**/	
	fscanf(fp,"%s%s[^\n]",ms,space);	
	if (strcmp(space,"hamming") == 0) 		flag = 0;
	else if (strcmp(space,"matrix") == 0) 	flag = 3;
	else {
		/**Read second line of input_file**/
		fscanf(fp,"%s %s[^\n]",m,metric);
		if (strcmp(metric,"euclidean") == 0)	flag = 1;
		else  flag = 2;	
	}
	/**Get info from configuration file**/
	if (flag == 0) 	lines++;
	info = get_config_info(fc, lines);
	info->complete = complete;
	if (flag == 3) 	matrix_medoid(fp, fe, info, ini, assi, upd, clara);
	else if (flag == 0) 	hamming_medoid(fp, fe, info, ini, assi, upd, clara);
	else if ((flag == 1) || (flag == 2))	vector_medoid(fp, fe, info, ini, assi, upd, clara, flag);
	free(info);
	fclose(fp);
	fclose(fc);
	fclose(fe);
	return 0;
}

