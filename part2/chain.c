#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "chain.h"
#define ITEM_ID 15
#define TRICK 6

void insert_chain(char * key, void *v, chainp *pointer, int flag, int d, int id, int position)
{
	chainp temp;
	int i;
	temp = *pointer;
	/**if list is empty, put the first node**/
	if (temp == NULL)		
	{
		if (!flag)		/**Hamming**/
		{
			char *value = (char *)v;
			int size = strlen(value);
			char *end;
			temp = malloc(sizeof(chain));
			temp-> p = NULL;
			temp->key = malloc((strlen(key)+1)*sizeof(char));
			strcpy(temp->key,key);
			temp->position = position;
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
		else if(flag == 3) /**Matrix**/
		{
			temp = malloc(sizeof(chain));
			temp->key = malloc((strlen(key)+1)*sizeof(char));
			strcpy(temp->key,key);
			temp->position = position;
			temp->value = NULL;	
			temp->p = NULL;
		}
		else 	/**Vectors**/
		{
			double *value = (double *)v;
			temp = malloc(sizeof(chain));
			temp->key = malloc((strlen(key)+1)*sizeof(char));
			strcpy(temp->key ,key);
			temp->position = position;
			temp->p = malloc(d*sizeof(double));
			for(i=0; i < d; i++)	temp->p[i] = value[i];		
			if (flag != 2) 	temp->id = id;
			temp->value = NULL;
		}
		temp->next = NULL;
		*pointer = temp;
	}
	/**If list isn't empty, put new node at the end**/
	else
	{
		while(temp->next!=NULL) {
			temp = temp->next;
		}
		if (!flag)		/**Hamming**/
		{
			char *value = (char *)v;
			int size = strlen(value);
			char *end;
			temp->next = malloc(sizeof(chain));
			temp->next->p = NULL;
			temp->next->position = position;
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
		else if(flag == 3) /**Matrix**/
		{
			temp->next = malloc(sizeof(chain));
			temp->next->key = malloc((strlen(key)+1)*sizeof(char));
			strcpy(temp->next->key ,key);
			temp->next->position = position;
			temp->next->p = NULL;
			temp->next->value = NULL;
		}
		else 	/**Vectors**/
		{
			double *value = (double *)v;
			temp->next = malloc(sizeof(chain));
			temp->next->key = malloc((strlen(key)+1)*sizeof(char));
			strcpy(temp->next->key ,key);
			temp->next->position = position;
			if (flag != 2) 	temp->next->id = id;
			temp->next->value = NULL;
			temp->next->p = malloc(d*sizeof(double));
			for(i=0; i < d; i++)	temp->next->p[i] = value[i];
		}
		temp->next->next = NULL;
	}
}

int search_chain_NNR(chainp *b, void * qdata, double R, pointp *list, chainp *barrier, int flag, int euclID, int d, int *all)
{
	chainp temp;
	int done = 0, duplicate;
	centroid center;
	center.center = (void *)-1;
	temp = *b;
	if((flag == 0) || (flag == 3)) {
		int *q = (int *)qdata;
		int diff = -1, position;
		while (temp != NULL)
		{
			/**If barrier points to first chain node**/
			if (((*barrier) != NULL) && (strcmp(temp->key,(*barrier)->key) == 0)) {
				if (diff == -1)	(*all)++;
				break;
			}
			position = temp->position;
			diff = q[position];
			if (diff <= R) {
				if (done == 0)	done = 1;
				/**Insert in cluster**/
				duplicate = insert_points(list,temp->key,diff,-1,center,position);
				/**If this item is already in this cluster don't insert it again**/
				if (duplicate == 1) {
					temp = temp->next;
					continue;	
				}/**Done checking duplicate**/
				/**Add barrier**/
				if ((*barrier) == NULL)	*barrier = temp;
				/**Move inserted item to the end of the chain**/
				move_chain_nodes(b,temp);
			}
			temp = temp->next;
		}
	}
	else {
		double diff = -1.0;
		double *q = (double *)qdata;
		int exists = 0, position;
		/**Only for euclidean metric, count items for which ID(p) = ID(q)**/
		if(flag != 2) 
		{
			while (temp != NULL)
			{
				/**ID(p) = ID(q)**/
				if(euclID == temp->id) exists++;
				temp = temp->next;
			}
		}
		temp = *b;
		/**For cosine metric or euclidean if exists=0**/
		if(exists <= 1) 
		{
			while (temp != NULL)
			{
				/**If barrier points to first chain node**/
				if (((*barrier) != NULL) && (strcmp(temp->key,(*barrier)->key) == 0)) {
					if (diff == -1.0)	(*all)++;
					break;
				}
				position = temp->position;
				diff = q[position];
				if (diff <= R)	{
					if (done == 0)	done = 1;
					/**Insert in cluster**/
					duplicate = insert_points(list,temp->key,diff,-1,center,temp->position);
					/**If this item is already in this cluster don't insert it again**/
					if (duplicate == 1) {
						temp = temp->next;
						continue;	
					}/**Done checking duplicate**/
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
			while (temp != NULL)
			{
				/**ID(p) = ID(q)**/
				if(euclID == temp->id)  
				{
					/**If barrier points to first chain node**/
					if (((*barrier) != NULL) && (strcmp(temp->key,(*barrier)->key) == 0)) {
						if (diff == -1)	(*all)++;
						break;
					}
					position = temp->position;
					diff = q[position];
					if (diff <= R)	{
						if (done == 0)	done = 1;
						/**Insert in cluster**/
						duplicate = insert_points(list,temp->key,diff,-1,center,temp->position);
						/**If this item is already in this cluster don't insert it again**/
						if (duplicate == 1) {
							temp = temp->next;
							continue;	
						}/**Done checking duplicate**/
						/**Add barrier**/
						if ((*barrier) == NULL)	*barrier = temp;
						/**Move inserted item to the end of the chain**/
						move_chain_nodes(b,temp);
					}
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
		printf("key: %s",l->key);
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
		/**Hamming **/
		if (flag == 0)	free(temp->value);		
		/**Euclidean or Cosine**/
		else if ((flag == 1) || (flag == 2))	free(temp->p);		
		free(temp);
	}
	*l = NULL;	
}

int insert_points(pointp *list, char *key, double distance1, double distance2, centroid second, int position) {
	pointp temp;
	temp = *list;
	/**If list is empty, put the first node**/
	if (temp == NULL) {
		temp = malloc(sizeof(point));
		temp->key = malloc((strlen(key)+1)*sizeof(char));
		strcpy(temp->key,key);
		temp->position = position;
		temp->mindistance = distance1;
		temp->secdistance = distance2;
		temp->duplicate = 0;
		temp->second = second;
		temp->next = NULL;
		*list = temp;
	}
	/**If list isn't empty, put new node at the end**/
	else {
		while(temp->next != NULL) {
			/**Avoid duplicates**/
			if (strcmp(temp->key,key) == 0)	return 1; 	
			temp = temp->next;
		}
		if (strcmp(temp->key,key) == 0)	return 1;
		temp->next = malloc(sizeof(point));
		temp->next->key = malloc((strlen(key)+1)*sizeof(char));
		strcpy(temp->next->key,key);
		temp->next->position = position;
		temp->next->mindistance = distance1;
		temp->next->secdistance = distance2;
		temp->next->duplicate = 0;
		temp->next->second = second;
		temp->next->next = NULL;
	}
	return 0;
}

int chain_length(pointp l) 
{
	int length = 0;
	while (l != NULL) {
		length++;
		l = l->next;
	}
	return length;
}

void print_points(pointp l) 
{
	while (l != NULL) 
	{
		printf("\nkey: %s, mindistance: %.0f ",l->key,l->mindistance);
		printf(" second: cluster %d, secdistance: %.0f ",(int)(intptr_t)l->second.center,l->secdistance);
		l = l->next;
	}
}

void delete_from_chain(pointp *list, char *item) {
	pointp temp = *list, previous;
	if (temp == NULL) 	return;
	/**If the item is the first element of the list**/
	if (strcmp(temp->key,item) == 0) {
		*list = temp->next;	
		free(temp->key);
		free(temp);
		return;
	}
	previous = temp;
	temp = temp->next;
	while (temp != NULL) {
		if (strcmp(temp->key,item) == 0) {
			previous->next = temp->next;
			free(temp->key);
			free(temp);
			return;
		}
		previous = temp;
		temp = temp->next;
	}
}

void destroy_points(pointp *l) 
{
	pointp temp, curr;
	curr = *l;
	if (curr == NULL)		return;		
	while (curr != NULL) 
	{
		temp = curr;
		curr = curr->next;
		free(temp->key);		
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

pointp clone(pointp list) {
    if (list == NULL) return NULL;
    pointp result = malloc(sizeof(point));
    result->key = malloc((strlen(list->key)+1)*sizeof(char));
	strcpy(result->key,list->key);
	result->mindistance = list->mindistance;
	result->secdistance = list->secdistance;
	result->duplicate = list->duplicate;
	result->second = list->second;
	result->position = list->position;
    result->next = clone(list->next);
    return result;
}

int make_item(char *item)
{
	int key;
	/**If the first character of string is type of char**/
	char *id;
	int s = strlen(item) - 4;
	id = malloc(s+1);
	strncpy(id,item+4,s);
	id[s] = '\0';
	/**Keep only K of itemK**/
	key = atoi(id);		
	free(id);
	return key;
}


