#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inputProcessing.h"

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


int *create_items(FILE *fp, int *P, int *numofusers, int *numofitems) {
	int i, j, flag, userId, itemId, rate, previousu, size, *items;
	char p[4];
	previousu = -1;
	/**Read the number of NN, if given**/
	fscanf(fp,"%s",p);
	if (strcmp(p,"P:") == 0)	fscanf(fp,"%d[^\n]",P);
	else {
		*P = 20;
		fseek(fp,0,SEEK_SET);
	}
	items = malloc(STARTSIZE*sizeof(int));
	for (i=0; i < STARTSIZE; i++)	items[i] = -1;
	/**Find number of users and items**/
	size = STARTSIZE;
	while (fscanf(fp,"%d%d%d[^\n]",&userId,&itemId,&rate) != EOF) {
		if (previousu != userId) (*numofusers)++;
		if (items[size-1] == -1) {
			/**Check if item already exists**/
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
			/**Check if item already exists**/
			flag = 0;
			for (j=0; j < size; j++) {
				if (items[j] == itemId) {
					flag = 1;
					break;
				}
			}
			if (!flag) {
				items = realloc(items,(size+STARTSIZE)*sizeof(int));
				items[size] = itemId;
				for (j=size+1; j < size+STARTSIZE; j++)	items[j] = -1;
				size += STARTSIZE;
			}
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
	return items;
}

user *create_users(FILE *fp, int *items, int numofusers, int numofitems) {
	int i, j,tempP, num, previousu, userId, itemId;
	double rate, avrate;
	char p[4];
	user *users;
	previousu = -1;
	num = avrate = 0;
	users = malloc(numofusers*sizeof(user));
	for (i=0; i < numofusers; i++) {
		users[i].ratings = malloc(numofitems*sizeof(rating));
		for (j=0; j < numofitems; j++) {
			users[i].ratings[j].itemId = items[j];
			users[i].ratings[j].rate = 0;
		}
	}
	fseek(fp,0,SEEK_SET);
	fscanf(fp,"%s",p);
	if (strcmp(p,"P:") != 0)	fseek(fp,0,SEEK_SET);
	else fscanf(fp,"%d[^\n]",&tempP);
	i = 0;
	while (fscanf(fp,"%d%d%lf[^\n]",&userId,&itemId,&rate) != EOF) {
		if (previousu != userId && previousu != -1) {
			if (previousu == -1) users[i].userId = userId;
			else users[i].userId = previousu;
			avrate /= num;
			if (metric == 1) {
				/**Normalize ratings**/
				for (j=0; j < numofitems; j++)	{
					if (users[i].ratings[j].rate != 0) users[i].ratings[j].rate -= avrate;
				}
			}
			users[i].average = avrate;
			num = avrate = 0;
			i++;
		}
		for (j=0; j < numofitems; j++)	{
			if (users[i].ratings[j].itemId == itemId) {
				users[i].ratings[j].rate = rate;
				break;
			}
		}
		avrate += rate;
		num++;
		previousu = userId;
	}
	users[i].userId = previousu;
	avrate /= num;
	if (metric == 1) {
		for (j=0; j < numofitems; j++)	{
			if (users[i].ratings[j].rate != 0) users[i].ratings[j].rate -= avrate;
		}
	}
	users[i].average = avrate;
	return users;
}

void quickSort(int *arr, int low, int high) {
	int pi;
    if (low < high) {
        pi = partition(arr,low,high);
		/**Separately sort elements before partition and after partition**/
        quickSort(arr,low,pi-1);
        quickSort(arr,pi+1,high);
    }
}

void swap(int *a, int *b) {
    int t = *a;
    *a = *b;
    *b = t;
}
 
int partition(int *arr, int low, int high) {
	/**This function takes last element as pivot**/
    int j, pivot = arr[high];
    /**Index of smallest element**/
    int i = low - 1;  
    for (j=low; j <= high-1; j++) {
        /**If current element is smaller than pivot**/
        if (arr[j] < pivot) {
			/**Increment index of smaller element**/
            i++;    
            swap(&arr[i],&arr[j]);
        }
    }
    swap(&arr[i+1],&arr[high]);
    return (i + 1);
}



