#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "chain.h"
#include "distances.h"
#define ITEM_ID 15
#define TRICK 6

void insert_chain(char * key, void *v, chainp *pointer, int flag, int d, int id)
{
	chainp temp;
	int i;
	temp = *pointer;
	/*if list is empty, put the first node*/
	if (temp == NULL)		
	{
		if (!flag)		//Hamming
		{
			char *value = (char *)v;
			int size = strlen(value);
			char *end;
			temp = malloc(sizeof(chain));
			temp-> p = NULL;
			temp->key = malloc((strlen(key)+1)*sizeof(char));
			strcpy(temp->key,key);
			if ((size>32) && (size<=64))
			{
				temp->value = malloc(sizeof(uint64_t));
				*temp->value = strtoull(value,&end,2);
			}
			else if ((size>16) && (size<=32))
			{
				temp->value = malloc(sizeof(uint32_t));
				*temp->value = strtoul(value,&end,2);
			}
			else if ((size>8) && (size<=16))
			{
				temp->value = malloc(sizeof(uint16_t));
				*temp->value = strtoul(value,&end,2);
			}
			else if (size<=8)
			{
				temp->value = malloc(sizeof(uint8_t));
				*temp->value = strtoul(value,&end,2);
			}
		}
		else if(flag == 3) //Matrix
		{
			temp = malloc(sizeof(chain));
			temp->p = NULL;
			temp->key = malloc((strlen(key)+1)*sizeof(char));
			strcpy(temp->key,key);
			temp->value = NULL;	
		}
		else 	//Vectors
		{
			double *value = (double *)v;
			temp = malloc(sizeof(chain));
			temp->key = malloc((strlen(key)+1)*sizeof(char));
			strcpy(temp->key ,key);
			temp->p = malloc(d*sizeof(double));
			for(i=0; i < d; i++)
				temp->p[i] = value[i];
			
			if (flag != 2) 	temp->id = id;
			temp->value = NULL;
		}
		temp->next = NULL;
		*pointer = temp;
	}
	/*If list isn't empty, put new node at the end*/
	else
	{
		while(temp->next!=NULL)
			temp = temp->next;
		if (!flag)		//Hamming
		{
			char *value = (char *)v;
			int size = strlen(value);
			char *end;
			temp->next = malloc(sizeof(chain));
			temp->next->p = NULL;
			temp->next->key = malloc((strlen(key)+1)*sizeof(char));
			strcpy(temp->next->key,key);
			if ((size>32) && (size<=64))
			{
				temp->next->value = malloc(sizeof(uint64_t));
				*temp->next->value = strtoull(value,&end,2);
			}
			else if ((size>16) && (size<=32))
			{
				temp->next->value = malloc(sizeof(uint32_t));
				*temp->next->value = strtoul(value,&end,2);
			}
			else if ((size>8) && (size<=16))
			{
				temp->next->value = malloc(sizeof(uint16_t));
				*temp->next->value = strtoul(value,&end,2);
			}
			else if (size<=8)
			{
				temp->next->value = malloc(sizeof(uint8_t));
				*temp->next->value = strtoul(value,&end,2);
			}
		}
		else if(flag == 3) //Matrix
		{
			temp->next = malloc(sizeof(chain));
			temp->next->key = malloc((strlen(key)+1)*sizeof(char));
			strcpy(temp->next->key ,key);
			temp->next->p = NULL;
			temp->next->value = NULL;
		}
		else 	//Vectors
		{
			double *value = (double *)v;
			temp->next = malloc(sizeof(chain));
			temp->next->key = malloc((strlen(key)+1)*sizeof(char));
			strcpy(temp->next->key ,key);
			if (flag != 2) 	temp->next->id = id;
			temp->next->value = NULL;
			temp->next->p = malloc(d*sizeof(double));
			for(i=0; i < d; i++)
				temp->next->p[i] = value[i];
		}
		temp->next->next = NULL;
	}
}

int search_chain_NNR(chainp *b, void *qdata, double R, chainp *list, chainp *barrier, int flag, int euclID, int d, int *all)
{
	chainp temp, tmp;
	int done = 0, duplicate;
	temp = *b;
	if(!flag)
	{
		uint64_t num1, num2;
		char *end;
		char *q = (char *)qdata;
		int diff = 0;
		num2 = strtoull(q,&end,2);
		while (temp != NULL)
		{
			num1 = *(temp->value);   
			diff = distance_Hamming(num1,num2);
			if (diff <= R)
				insert_chain(temp->key,NULL,list,flag,0,0);
			temp = temp->next;
		}
	}
	else if(flag == 3)
	{
		int *q = (int *)qdata;
		int diff = -1, position;
		while (temp != NULL)
		{
			tmp = *list;
			/**If barrier points to first chain node**/
			if (((*barrier) != NULL) && (strcmp(temp->key,(*barrier)->key) == 0)) {
				if (diff == -1)	(*all)++;
				break;
			}
			position = make_item(temp->key);
			diff = q[position-1];
			if (diff <= R) 
			{
				/**If this item is already in this cluster don't insert it again**/
				duplicate = 0;
				while (tmp != NULL) {
					if (strcmp(temp->key, tmp->key) == 0) {
						duplicate = 1;
						break;
					}
					tmp = tmp->next;	
				}
				if (duplicate == 1) {
					temp = temp->next;
					continue;
				}/**Done checking duplicate**/
				if (done == 0)	done = 1;
				/**Insert in cluster**/
				insert_chain(temp->key,NULL,list,flag,0,0);
				/**Add barrier**/
				if ((*barrier) == NULL)	*barrier = temp;
				/**Move inserted item to the end of the chain**/
				move_chain_nodes(b,temp);
			}
			temp = temp->next;
		}
	}
	else
	{
		double diff;
		double *q = (double *)qdata;
		int exists = 0;
		/*Only for euclidean metric, count items for which ID(p) = ID(q)*/
		if(flag != 2) 
		{
			while (temp != NULL)
			{
				if(euclID == temp->id) 	//ID(p) = ID(q)
					exists++;
				temp = temp->next;
			}
		}
		temp = *b;
		/*For cosine metric or euclidean if exists=0*/
		if(exists <= 1) 
		{
			while (temp != NULL)
			{
				if(flag!=2)
					diff = distance_Euclidean(temp->p,q,d);
				else
					diff = distance_Cosine(temp->p,q,d);
				if (diff <= R)
					insert_chain(temp->key,NULL,list,flag,0,0);
				temp = temp->next;
			}
		}
		else
		{
			while (temp != NULL)
			{
				if(euclID == temp->id)  //ID(p) = ID(q)
				{
					diff = distance_Euclidean(temp->p,q,d);
					if (diff <= R)
						insert_chain(temp->key,NULL,list,flag,0,0);
				}
				temp = temp->next;
			}
		}
	}
	return done;
}

void move_chain_nodes(chainp *b, chainp temp) {
	chainp l = *b;
	/**If list is empty or only one item inside, then return it**/
	if ((l == NULL) || (l->next == NULL))	return;
	/**If item wanted is the first node**/
	if (l == temp) {
		/**New start is the second item**/
		*b = l->next;
		/**Just go to last node**/
		while (l->next != NULL) l = l->next;
	}
	/**If item wanted is after the first one**/
	else {
		while (l->next != NULL) {
			if (l->next == temp) {
				l->next = l->next->next;
				if (l->next == NULL) break;
			}
			l = l->next;
		}
	}
	l->next = temp;
	l->next->next = NULL;
}

void print_chain(chainp l) 
{
	while (l != NULL) 
	{
		printf("key: %s,",l->key);
		l = l->next;
	}
}

void destroy_chain(chainp *l, int flag) 
{
	chainp temp, curr;
	curr = *l;
	if (curr == NULL)		return;		
	while (curr != NULL) 
	{
		temp = curr;
		curr = curr->next;
		free(temp->key);
		if (flag == 0)	free(temp->value);		//Hamming 
		else if ((flag == 1) || (flag == 2))	free(temp->p);		//Euclidean or Cosine
		free(temp);
	}
	*l = NULL;	
}

/*void print_nnrlist(nnrp *l, FILE *fe) 
{
	nnrp temp, curr;
	curr = *l;
	if (curr == NULL)		return;		//For safety
	while (curr != NULL) 
	{
		fprintf(fe,"%s\n",curr->neighbor.key);
		temp = curr;
		curr = curr->next;
		free(temp);
	}
	*l = NULL;
}*/

int make_item(char *item)
{
	int key;
	/*If the first character of string is type of char*/
	if (!isdigit(item[0]))	
	{	
		char *id;
		int s = strlen(item) - 4;
		id = malloc(s+1);
		strncpy(id,item+4,s);
		key = atoi(id);		//Keep only K of item_idK
		free(id);
	}
	else 	key = atoi(item);
	return key;
}


