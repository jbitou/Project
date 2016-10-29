#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "chain.h"
#include "distances.h"
#define ITEM_ID 15

void insert_chain(char * key, void *v, chainp *pointer, int flag, int d, int id)
{
	chainp temp;
	int i;
	temp = *pointer;
	if (temp == NULL)		//if list is empty, put the first node
	{
		if (!flag)		//Hamming
		{
			char *value = (char *)v;
			int size = strlen(value);
			char *end;
			temp = malloc(sizeof(chain));
			temp-> p = NULL;
			temp->key = malloc((strlen(key)+1)*sizeof(char));
			strcpy(temp->key ,key);
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
		else if(flag == 3) //matrix
		{
			int *value = (int *)v;
			temp = malloc(sizeof(chain));
			temp->p = NULL;
			temp->id = d;
			temp->key = malloc((strlen(key)+1)*sizeof(char));
			strcpy(temp->key ,key);
			temp->value = NULL;
			if(d != 0) 
			{
				temp->distances = malloc(d*sizeof(int));
				for(i=0; i < d; i++)
					temp->distances[i] = value[i];
			}
			
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
			strcpy(temp->next->key ,key);
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
		else if(flag == 3) //matrix
		{
			int *value = (int *)v;
			temp->next = malloc(sizeof(chain));
			temp->next->key = malloc((strlen(key)+1)*sizeof(char));
			strcpy(temp->next->key ,key);
			temp->next->id = d;
			temp->next->p = NULL;
			temp->next->value = NULL;
			if(d != 0) 
			{
				temp->next->distances = malloc(d*sizeof(int));
				for(i=0; i < d; i++)
					temp->next->distances[i] = value[i];
			}
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

void search_chain_NNR(chainp b, void *qdata, double R, nnrp *nnrlist, int flag, int euclID, int d)
{
	chainp temp;
	temp = b;
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
				insertnnrlist(temp->key,nnrlist);
			temp = temp->next;
		}
	}
	else if(flag == 3)
	{
		int *q = (int *)qdata;
		int diff = 0, position;
		while (temp != NULL)
		{
			position = make_item(temp->key);
			diff = q[position-1];
			if (diff <= R)
			{
				insertnnrlist(temp->key,nnrlist);
			}
			temp = temp->next;
		}
	}
	else
	{
		double diff;
		double *q = (double *)qdata;
		int exists = 0;
		if(flag!=2) //only for euclidean metric
		{
			while (temp != NULL)
			{
				if(euclID == temp->id) { //ID(p) = ID(q)
					exists++;
				}
				temp = temp->next;
			}
		}
		temp = b;
		if(exists <= 1) 
		{
			while (temp != NULL)
			{
				if(flag!=2)
					diff = distance_Euclidean(temp->p,q,d);
				else
					diff = distance_Cosine(temp->p,q,d);
				//printf("diff between %s and query is: %0.20f\n",temp->key,diff);
				if (diff <= R)
					insertnnrlist(temp->key,nnrlist);
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
					//printf("ID is:%d ----diff between %s and query is: %0.20f\n",euclID,temp->key,diff);
					if (diff <= R)
						insertnnrlist(temp->key,nnrlist);
				}
				temp = temp->next;
			}
		}
	}
}

nn search_chain_NN(chainp b, void *q, int flag, int bruflag, int euclID, int d, int *counter, int L)
{
	chainp temp,tmp;
	nn lshnn;
	lshnn.key = malloc(ITEM_ID);
	if (!flag)
	{
		uint64_t num1, num2;
		char *end;
		int diff = 0, diff1;
		char key[ITEM_ID];
		char *qdata = (char *)q;
		num2 = strtoull(qdata,&end,2);
		/*bucket is empty (only in brute force)*/
		if (b == NULL)
		{
			lshnn.distance = -1;
			return lshnn;
		}
		temp = b;
		num1 = *(temp->value);   
		diff = distance_Hamming(num1,num2);
		strcpy(key,temp->key);
		temp = temp->next;
		while (temp != NULL)
		{
			if (!bruflag)
			{
				if(*counter > 3*L)
				{
					strcpy(lshnn.key,key);
					lshnn.distance = (double)diff;
					return lshnn;
				}
			}
			num1 = *(temp->value);   
			diff1 = distance_Hamming(num1,num2);
			if (!bruflag) *counter= *counter+1;
			if(((diff1 < diff) && (diff1 > 0)) || (diff <= 0)) 
			{
				diff = diff1;
				strcpy(key,temp->key);
			}
			temp = temp->next;
		}
		strcpy(lshnn.key,key);
		lshnn.distance = (double)diff;
	}
	else if(flag == 3)
	{
		int diff = 0, position, diff1;
		char key[ITEM_ID];
		int *qdata = (int *)q;
		/*bucket is empty (only in brute force)*/
		if (b == NULL)
		{
			lshnn.distance = -1;
			return lshnn;
		}
		temp = b;
		position = make_item(temp->key);
		diff = qdata[position-1];
		strcpy(key,temp->key);
		temp = temp->next;
		while (temp != NULL)
		{
			if (!bruflag)
			{
				if (*counter > 3*L) 
				{
					strcpy(lshnn.key,key);
					lshnn.distance = (double)diff;
					return lshnn;
				}
			}
			position = make_item(temp->key);
			diff1 = qdata[position-1];
			//printf("2.key=%s  diff=%d--distance=%f\n",temp->key,diff,*distance);
			if (!bruflag) *counter = *counter+1;
			if(((diff1 < diff) && (diff1 > 0)) || (diff <= 0)) 
			{
				diff = diff1;
				strcpy(key,temp->key);
			}
			temp = temp->next;
		}
		strcpy(lshnn.key,key);
		lshnn.distance = (double)diff;
	}
	else
	{
		double diff = 0.0, diff1;
		int exists = 0;
		char key[ITEM_ID];
		double *qdata = (double *)q;
		int i;
		/*bucket is empty (only in brute force)*/
		if (b == NULL)
		{
			lshnn.distance = -1;
			return lshnn;
		}
		temp = b;
		tmp = b;
		if(flag != 2)
		{
			while (tmp != NULL)
			{
				
				if(euclID == tmp->id)  //ID(p) = ID(q)
					exists++;
				tmp = tmp->next;
			}
		}
		if(exists <= 1)  
		{
			
			if(flag != 2)
				diff = distance_Euclidean(qdata,temp->p,d);
			else
				diff = distance_Cosine(qdata,temp->p,d);	
			strcpy(key,temp->key); 
			/*apo panw kalutera??????*/
			if (!bruflag)
			{
				*counter = *counter+1;
				if(*counter > 6*L)
				{
					strcpy(lshnn.key,key);
					lshnn.distance = (double)diff;
					return lshnn;
				}
			}
		}
		else 
		{
			if(euclID == temp->id)
			{
				diff = distance_Euclidean(qdata,temp->p,d);
				strcpy(key,temp->key);
				if (!bruflag)
				{
					*counter= *counter+1;
					if(*counter > 6*L)
					{
						strcpy(lshnn.key,key);
						lshnn.distance = (double)diff;
						return lshnn;
					}
				}
			}
		}	
		temp = temp->next;
		if(exists <= 1)  
		{
			while (temp != NULL)
			{ 
				if(flag != 2)
					diff1 = distance_Euclidean(qdata,temp->p,d);
				else
					diff1 = distance_Cosine(qdata,temp->p,d);
				if(((diff1 < diff) && (diff1 > 0)) || (diff <= 0)) 
				{
					diff = diff1;
					strcpy(key,temp->key);
					
				}
				if (!bruflag)
				{
					*counter= *counter+1;
					if(*counter > 6*L)
					{
						strcpy(lshnn.key,key);
						lshnn.distance = (double)diff;
						return lshnn;
					}
				}
				temp = temp->next;
			}
		}
		else 
		{
			while (temp != NULL)
			{ 
				if(euclID == temp->id)
				{
					diff1 = distance_Euclidean(qdata,temp->p,d);
					if(((diff1 < diff)  && (diff1 > 0)) || (diff <= 0))  
					{
						diff = diff1;
						strcpy(key,temp->key);
					}
					if (!bruflag)
					{
						*counter= *counter+1;
						if(*counter > 6*L)
						{
							strcpy(lshnn.key,key);
							lshnn.distance = (double)diff;
							return lshnn;
						}
					}
				}
				temp = temp->next;
			}
		}
		strcpy(lshnn.key,key);
		lshnn.distance = diff;
	}
	return lshnn;
}

int make_item(char *item)
{
	int key;
	if (!isdigit(item[0]))	//If the first character of string is type of char
	{	
		char *id;
		int s = strlen(item) - 4;
		id = malloc(s+1);
		strncpy(id,item+4,s);
		key = atoi(id);		//keep only K of item_idK
		free(id);
	}
	else 	key = atoi(item);		//else, just convert type char* to int *
	return key;
}

void destroy_chain(chainp *l, int flag) 
{
	chainp temp, curr;
	curr = *l;
	if (curr == NULL)		return;		//for safety
	while (curr != NULL) 
	{
		temp = curr;
		curr = curr->next;
		free(temp->key);
		if (flag == 0)	free(temp->value);	//Hamming 
		else if(flag == 3) free(temp->distances);	//Matrix
		else if ((flag == 1) || (flag == 2))	free(temp->p);	//Euclidean or Cosine
		free(temp);
	}
	*l = NULL;	//for safety
}

void print_chain(chainp l)
{
	//uint64_t d;
	int i;
	while (l != NULL)
	{
		printf("key: %s, ",l->key);
		/*for(i=0;i<l->id;i++)
			printf("distance[%d]=%d ",i,l->distances[i]);
		printf("\n");	*/
		//printf("value: %" PRIu64 "\n",*(l->value));
		//d = decimalToBinary(*(l->value));
		//printf("value bits: %" PRIu64 "\n",d);
		l=l->next;
	}
}
