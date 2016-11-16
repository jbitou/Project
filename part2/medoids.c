#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrixMedoid.h"
#define ITEM_ID 15

int main (int argc, char **argv)
{
	FILE *fp, *fc, *fe;
	int i, input, config, output, ch, lines = -2, flag;
	char ms[14], space[10], m[10], metric[20];
	pinfo info;
	
	if (command_processing(argc) == -1)	return -1;
	input = config = output = -1;
	for (i=1; i < (argc-1); i+=2)
	{
		if (strcmp(argv[i],"-d") == 0) 			input = i+1;
		else if (strcmp(argv[i],"-c") == 0) 	config = i+1;
		else if (strcmp(argv[i],"-o") == 0)		output = i+1;
	}
	if ((input == -1) || (config == -1) || (output == -1))
	{
		printf("Arguments -d,-c,-o are required. Try again.\n");
		return -1;
	}
	/**Open input file**/
	fp = fopen(argv[input],"r");
	if (fp == NULL)
	{
		perror("Error opening input_file");
		return -1;
	}
	/**Open configuration file**/
	fc = fopen(argv[config],"r");
	if (fc == NULL)
	{
		perror("Error opening configuration_file");
		return -1;
	}
	/**Open output file**/
	fe = fopen(argv[output],"w+");
	if (fe == NULL)
	{
		perror("Error opening output_file");
		return -1;
	}
	/**Count lines of input file**/
	while ((ch = fgetc(fp)) != EOF) {
		if (ch == '\n') lines++;
	}
	/**Lines found**/
	fseek(fp,0,SEEK_SET);
	/**Read first line of input_file**/	
	fscanf(fp,"%s%s[^\n]",ms,space);	
	if (strcmp(space,"hamming") == 0) 	flag = 0;
	else if (strcmp(space,"matrix") == 0) flag = 3;
	else {
		fscanf(fp,"%s %s[^\n]",m,metric);
		if (strcmp(metric,"euclidean") == 0)	flag = 1;
		else  flag = 2;	
	}
	/**Get info from configuration file**/
	if (flag == 0) 	lines++;
	printf("lines=%d\n",lines);
	info = get_config_info(fc, lines);
	printf("k=%d\n",info->k);
	printf("hash function=%d\n",info->num_of_hash);
	printf("hash tables=%d\n",info->L);
	printf("fraction=%d\n",info->fraction);
	printf("iterations=%d\n",info->iterations);
	if (flag == 3) 	matrix_medoid(fp, info);
	free(info);
	fclose(fp);
	fclose(fc);
	fclose(fe);
	return 0;
}

