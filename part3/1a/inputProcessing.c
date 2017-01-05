#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inputProcessing.h"
#define STARTSIZE 10

int command_processing(char **argv, int argc, int *input, int *output, int *validate) {
	int i;
	if (argc > 6) {
		printf("Too many arguments. Try again.\n");
		return -1;
	}
	else if (argc < 5) {
		printf("Too few arguments. Try again.\n");
		return -1;
	}
	/**Check the position of input and output files**/
	for (i=1; i < (argc-1); i+=2) {
		if (strcmp(argv[i],"-d") == 0) 			*input = i+1;
		else if (strcmp(argv[i],"-o") == 0)		*output = i+1;
	}
	/**Parameter -validate may be given alone at the end**/
	if (strcmp(argv[argc-1],"-validate") == 0)	*validate = 1;
	/**Check if all parameters are given**/
	if ((*input == -1) || (*output == -1)) {
		printf("Arguments -d,-o are required. Try again.\n");
		return -1;
	}
	return 0;	
}


int count_data(FILE *fp, int *numofusers, int *numofitems) {
	int i, j, flag, userId, itemId, rate, previousu, *items, size, P;
	char p[4];
	previousu = -1;
	/**Read the number of NN, if given**/
	fscanf(fp,"%s",p);
	if (strcmp(p,"P:") == 0)	fscanf(fp,"%d[^\n]",&P);
	else {
		P = 20;
		fseek(fp,0,SEEK_SET);
	}
	items = malloc(STARTSIZE*sizeof(int));
	for (i=0; i < STARTSIZE; i++)	items[i] = -1;
	/**Find number of users and items**/
	size = STARTSIZE;
	while (fscanf(fp,"%d%d%d[^\n]",&userId,&itemId,&rate) != EOF) {
		if (previousu != userId) (*numofusers)++;
		if (items[size-1] == -1) {
			j = flag = 0;
			while (items[j] != -1) {
				if (items[j] == itemId) {
					flag = 1;
					break;
				}
				j++;
			}
			if (!flag) items[j] = itemId;
		}	
		else {
			items = realloc(items,(size+STARTSIZE)*sizeof(int));
			items[size] = itemId;
			for (j=size+1; j < size+STARTSIZE; j++)	items[j] = -1;
			size += STARTSIZE;
		}	
		previousu = userId;
	}
	if (items[size-1] != -1)	*numofitems = size;
	else {
		i = 0;
		while (items[i] != -1) {
			i++;
			(*numofitems)++;
		}
	}
	for (i=0; i < size; i++) printf("items[%d]=%d\n",i,items[i]);
	free(items);
	return P;
}

user *create_users(FILE *fp, int numofusers, int numofitems) {
	int i, j, num, avrate, previousu, userId, itemId, rate;
	char p[4];
	user *users;
	previousu = -1;
	num = avrate = 0;
	users = malloc(numofusers*sizeof(user));
	for (i=0; i < numofusers; i++) {
		users[i].ratings = malloc(numofitems*sizeof(rating));
		for (j=0; j < numofitems; j++) {
			users[i].ratings[j].itemId = 0;
			users[i].ratings[j].rate = 0;
		}
	}
	fseek(fp,0,SEEK_SET);
	fscanf(fp,"%s",p);
	if (strcmp(p,"P:") != 0)	fseek(fp,0,SEEK_SET);
	i = 0;
	while (fscanf(fp,"%d%d%d[^\n]",&userId,&itemId,&rate) != EOF) {
		if (previousu != userId && previousu != -1) {
			if (previousu == -1) users[i].userId = userId;
			else users[i].userId = previousu;
			avrate /= num;
			printf("user %d has average rate %d\n",previousu,avrate);
			/**Normalize ratings**/
			for (j=0; j < num; j++)	users[i].ratings[j].rate -= avrate;
			num = avrate = 0;
			i++;
		}
		users[i].ratings[num].itemId = itemId;
		users[i].ratings[num].rate = rate;
		avrate += rate;
		num++;
		previousu = userId;
	}
	users[i].userId = previousu;
	avrate /= num;
	printf("user %d has average rate %d\n",previousu,avrate);
	for (j=0; j < num; j++)	users[i].ratings[j].rate -= avrate;
	return users;
}

