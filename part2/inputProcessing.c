#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inputProcessing.h"

int command_processing(int argc) {
	if (argc > 7)
	{
		printf("Too many arguments. Try again.\n");
		return -1;
	}
	else if (argc < 7)
	{
		printf("Too few arguments. Try again.\n");
		return -1;
	}
	/**Every parameter has to be given after a recogniser**/
	if ((argc % 2) == 0)
	{
		printf("Wrong arguments. Try again.\n");	
		return -1;
	}
	return 0;	
}

pinfo get_config_info(FILE *fc, int N) {
	char infoword[50];
	int i, ch, lines = 0, temp;
	pinfo info = malloc(sizeof(conf_info));
	/**Initialize struct**/
	info->num_of_hash  = 0;
	info->L = 0;
	info->fraction = 0;
	info->iterations = 0;
	/**Count lines**/
	while ((ch = fgetc(fc)) != EOF) {
		if (ch == '\n') lines++;
	}
	/**Lines found**/
	fseek(fc,0,SEEK_SET);
	/**Always read number of clusters**/
	fscanf(fc, "%s %d", infoword, &(info->k));
	if (lines == 1) {
		info->num_of_hash  = 4;
		info->L = 5;
		info->fraction =  0.12*info->k*(N-info->k);
		if (info->fraction < 250)	info->fraction = 250;
		info->iterations = 2;
	}
	else {
		for (i=0; i<(lines-1); i++) {
			fscanf(fc, "%s %d", infoword, &temp);
			if (strcmp(infoword,"number_of_hash_functions:") == 0) info->num_of_hash = temp;
			else if (strcmp(infoword,"number_of_hash_tables:") == 0) info->L = temp;
			else if (strcmp(infoword,"clarans_set_fraction:") == 0) info->fraction = temp;
			else if (strcmp(infoword,"clarans_iterations:") == 0) info->iterations = temp;
		}
		if (info->num_of_hash == 0)  info->num_of_hash = 4;
		if (info->L == 0)  info->L = 5;
		if (info->fraction == 0)  {
			info->fraction = 0.12*info->k*(N-info->k);
			if (info->fraction < 250)	info->fraction = 250;
		}
		if (info->iterations == 0)  info->iterations = 2;
	}
	return info;
}

char *inputString(FILE *fp, size_t size) {
	char *str;
	int ch;
	size_t len = 0;
	str = realloc(NULL,size*sizeof(char));
	if(!str)  	return str;
	while ((ch = fgetc(fp)) != EOF && ch != '\n') 
	{
		str[len++] = ch;
		if(len == size) 
		{
			str = realloc(str,(size += (MAX_LINE / 2))*sizeof(char));
			if(!str) 	return str;
		}
	}
	str[len++] = '\0';
	return realloc(str,len*sizeof(char)); 
}

