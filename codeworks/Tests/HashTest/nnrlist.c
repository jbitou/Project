#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nnrlist.h"

void insert_nnrlist(char * key, nnrp *pointer)
{
	nnrp temp;
	temp = *pointer;
	/*If list is empty, put the first node*/
	if (temp == NULL)	
	{
		temp = malloc(sizeof(nnr));
		strcpy(temp->key,key);
		temp->next = NULL;
		*pointer = temp;
	}
	/*If list isn't empty, put new node at the end*/
	else
	{
		while(temp->next != NULL)
		{
			/*Avoid duplicates*/
			if (strcmp(temp->key,key) == 0)	return; 	
			temp = temp->next;
		}
		if (strcmp(temp->key,key) == 0)	return;
		temp->next = malloc(sizeof(nnr));
		strcpy(temp->next->key,key);
		temp->next->next = NULL;
	}
}

void destroy_nnrlist(nnrp *l, FILE *fe) 
{
	nnrp temp, curr;
	curr = *l;
	if (curr == NULL)		return;		
	while (curr != NULL) 
	{
		fprintf(fe,"%s\n",curr->key);	//First print in file
		temp = curr;
		curr = curr->next;
		free(temp);		//Then free the node
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
