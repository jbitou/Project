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

void print_nnrlist(nnrp *l, FILE *fe) 
{
	nnrp temp, curr;
	curr = *l;
	if (curr == NULL)		return;		//For safety
	while (curr != NULL) 
	{
		fprintf(fe,"%s\n",curr->key);
		temp = curr;
		curr = curr->next;
		free(temp);
	}
	*l = NULL;
}

void combine_nnrlist(nnrp *l1, nnrp *l2)
{
	nnrp temp;
	while (*l2 != NULL)
	{
		insert_nnrlist((*l2)->key,l1);
		temp = *l2;
		(*l2) = (*l2)->next;
		free(temp);
	}
}

void destroy_nnrlist(nnrp *l) 
{
	nnrp temp, curr;
	curr = *l;
	if (curr == NULL)		return;		//For safety
	while (curr != NULL) 
	{
		temp = curr;
		curr = curr->next;
		free(temp);
	}
	*l = NULL;
}

void display_nnrlist(nnrp l)
{
	while (l != NULL)
	{
		printf("key: %s ",l->key);
		l=l->next;
	}
	printf("\n");
}
