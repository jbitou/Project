#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nnrlist.h"

void insertnnrlist(char * key, nnrp *pointer)
{
	nnrp temp;
	temp = *pointer;
	if (temp == NULL)	//if list is empty, put the first node
	{
		temp = malloc(sizeof(nnr));
		//temp->key = malloc((strlen(key)+1)*sizeof(char));
		strcpy(temp->key ,key);
		temp->next = NULL;
		*pointer = temp;
	}
	else
	{
		while(temp->next!=NULL)
		{
			if (strcmp(temp->key,key) == 0)	return;
			temp = temp->next;
		}
		if (strcmp(temp->key,key) == 0)	return;
		temp->next = malloc(sizeof(nnr));
		//temp->next->key = malloc((strlen(key)+1)*sizeof(char));
		strcpy(temp->next->key ,key);
		temp->next->next = NULL;
	}
}

void destroy_nnrlist(nnrp *l, FILE *fe) 
{
	nnrp temp, curr;
	curr = *l;
	if (curr == NULL)		return;		//for safety
	while (curr != NULL) 
	{
		//printf("key is %s\n",curr->key);
		fprintf(fe,"%s\n",curr->key);
		temp = curr;
		curr = curr->next;
		//free(temp->key);
		free(temp);
	}
	*l = NULL;
}

void print_nnrlist(nnrp l)
{
	while (l != NULL)
	{
		printf("key: %s ",l->key);
		l=l->next;
	}
	printf("\n");
}
