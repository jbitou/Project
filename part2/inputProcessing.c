#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inputProcessing.h"

int command_processing(int argc) {
	/**Parameter -complete may be given alone at the end**/
	if (argc > 8) {
		printf("Too many arguments. Try again.\n");
		return -1;
	}
	else if (argc < 7) {
		printf("Too few arguments. Try again.\n");
		return -1;
	}
	/**Every parameter has to be given after a recogniser**/
	if (((argc % 2) == 0) && (argc != 8)) {
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

void user_choice(int *ini, int *assi, int *upd) {
	char answer[3];
	int counter;
	printf("Choose specific algorithms to be run? Y or N?  ");
	scanf("%s",answer);
	/**If user answered positively**/
	if ((strcmp(answer,"Y") == 0) || (strcmp(answer,"y") == 0)) {
		printf("\nChoose Initialization algorithm:\nFor K-medoids++ click 1.\nFor Park-Jun click 2.\n");
		counter = 0;
		do {
			if (counter > 0)	printf("Try again:\nFor K-medoids++ click 1.\nFor Park-Jun click 2.\n");
			scanf("%d",ini);
			counter++;
		}while ((*ini != 1) && (*ini != 2));
		printf("\nChoose Assignment algorithm:\nFor PAM assignment click 1.\nFor assignment by LSH/DBH click 2.\n");
		counter = 0;
		do {
			if (counter > 0)	printf("Try again:\nFor PAM assignment click 1.\nFor assignment by LSH/DBH click 2.\n");
			scanf("%d",assi);
			counter++;
		}while ((*assi != 1) && (*assi != 2));
		printf("\nChoose Update algorithm:\nFor Update a la Lloyd’s click 1.\nFor CLARANS click 2.\n");
		counter = 0;
		do {
			if (counter > 0)	printf("Try again:\nFor Update a la Lloyd’s click 1.\nFor CLARANS click 2.\n");
			scanf("%d",upd);
			counter++;
		}while ((*upd != 1) && (*upd != 2));
	}
	/**If user wants to see all combinations**/
	else printf("\nRunning all combinations...\n");
}

char *inputString(FILE *fp, size_t size) {
	char *str;
	int ch;
	size_t len = 0;
	str = realloc(NULL,size*sizeof(char));
	if(!str)  	return str;
	while ((ch = fgetc(fp)) != EOF && ch != '\n') {
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

