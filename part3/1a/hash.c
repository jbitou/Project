#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hash.h"
#define ITEM_ID 15
#define WINDOW_SIZE 4
#define TRICK 6

int mod(int a, long long b)
{
	if (b < 0) 		 return mod(a,-b);
	int ret = a % b;
	if(ret < 0) 	ret += b;
	return ret;
}

void init_table(int k, hash_table *htable, int tableSize) 
{
	int i;
	htable->size = tableSize;
	htable->table = malloc((htable->size)*sizeof(chainp));
	for(i=0; i < htable->size; i++)	
		htable->table[i] = NULL;
	return;
}

void init_hash_Ham(ghashp *g, int L, int k, char *data)
{
	int i, j, M = 1, N = strlen(data);
	for (i=0; i < L; i++) 
	{
		for (j=0; j < k; j++)
		{ 
			/*Choose uniformly an h function (position of a bit)*/
			g[i][j].t = M + (rand() / (RAND_MAX + 1.0)) * (N-M+1);
			g[i][j].v = NULL;
		}
	}
}

void init_hash_Eucl(ghashp * g, int L, int k, int d)
{
	int i,j,z, M = -1, N = 1, window = WINDOW_SIZE;
	double x,w,u,v,mult;
	/*For each table*/
	for (i=0; i < L; i++) 
	{
		/*Create k times h()*/
		for (j=0; j < k; j++)	
		{ 
			g[i][j].v = malloc(d*sizeof(double));
			/*Create a vector~N(0,1)*/
			for (z=0; z < d; z++)
			{ 
				do
				{
					u = M + (rand() / (RAND_MAX + 1.0)) * (N-M);
					v = M + (rand() / (RAND_MAX + 1.0)) * (N-M);
					w = pow (u,2) + pow (v,2);
				}
				while (w >= 1 || w == 0);
				mult = sqrt ((-2 * log (w)) / w);
				x = u * mult;
				g[i][j].v[z] = x;
			}
			g[i][j].t = (rand() / (RAND_MAX + 1.0)) * window;
			/*Choose r in range 0-128*/
			g[i][j].r = (rand() / (RAND_MAX + 1.0)) * 128;
		}
	}
}

void init_hash_Cos(ghashp * g, int L, int k, int d)
{
	int i,j,z, M = -1, N = 1;
	double x,w,u,v,mult;
	/*For each table*/
	for (i=0; i < L; i++) 	
	{
		/*Create k times h()*/
		for (j=0; j < k; j++) 
		{ 
			/*At cosine: v is r because r is real*/
			g[i][j].v = malloc(d*sizeof(double));	
			for (z=0; z < d; z++)
			{ 
				do
				{
					u = M + (rand() / (RAND_MAX + 1.0)) * (N-M);
					v = M + (rand() / (RAND_MAX + 1.0)) * (N-M);
					w = pow (u,2) + pow (v,2);
				}
				while (w >= 1 || w == 0);
				mult = sqrt ((-2 * log (w)) / w);
				x = u * mult;
				g[i][j].v[z] = x;
			}
			g[i][j].t = 0;
			g[i][j].r = 0;
		}
	}
}

int hash_func_Ham(ghashp g, char *data, int k)
{
	int i, h;
	char *end;
	char *temp = malloc(k*sizeof(char));
	for (i=0; i < k; i++)
		temp[i] = data[g[i].t];
	h = strtol(temp,&end,2);
	free(temp);
	return h;
}

int hash_func_Eucl(ghashp g, user usr, int k, int d) {
	int i = 0, j, window = WINDOW_SIZE, h, rh, sum = 0;
	long long M;
	double inner;
	M = (1LL << 32) - 5;	
	/**For each h()**/
	for (i=0; i < k; i++) {
		inner = 0.0;
		for (j=0; j < d; j++)	inner += g[i].v[j]*usr.ratings[j].rate;
		inner += g[i].t;
		inner /= window;
		h = (int)inner;
		rh = h*g[i].r;
		sum += rh;
	}
	sum = mod(sum,M);
	return sum;		
}

int hash_func_Cos(ghashp g, double *x, int k, int d)
{
	double inner;
	int i,j,h;
	char *end;
	char *temp = malloc(k*sizeof(char));
	/*For each h()*/
	for (i=0; i < k; i++)	
	{
		inner = 0.0;
		for (j=0; j < d; j++)	//Inner product
			inner += g[i].v[j]*x[j];	//ri*x
		if (inner >= 0.0)
			temp[i] = '1';
		else
			temp[i] = '0';	
	}
	h = strtol(temp,&end,2);
	free(temp);
	return h;
}

hash_table *vector_hash_insert(ghashp *g, user *users, int L, int k, int numofusers, int numofitems) {
	int i, j, euclID, pos, tableSize;
	hash_table *htable;
	tableSize = numofusers/16 + 1;
	/**Allocate memory for hash tables and hash functions**/
	htable = malloc(L * sizeof(hash_table));	
	for(i=0; i < L; i++)	init_table(k,&htable[i],tableSize);
	/**Insert users to hash tables**/
	for(i = 0; i < numofusers; i++)	{
		for(j = 0; j < L; j++)	{
			euclID = hash_func_Eucl(g[j],users[i],k,numofitems);
			euclID = abs(euclID);
			pos = mod(euclID,tableSize);
			vector_insert_chain(&(htable[j].table[pos]),users[i],numofitems,euclID);
		}
	}
	return htable;
}

nnrlist search_table_NNR(hash_table *htable, ghashp *g, nnrlist list, user usr, double R1, double R2, int numofusers, int numofitems, int L, int k) {
	int i, euclID, pos, tableSize;
	tableSize = numofusers/16 + 1;
	for(i = 0; i < L; i++)	{
		euclID = hash_func_Eucl(g[i],usr,k,numofitems);
		euclID = abs(euclID);
		pos = mod(euclID,tableSize);
		list = search_chain_NNR(htable[i].table[pos],list,usr,R1,R2,euclID,numofitems);
	}
	return list;
}

void destroy_table(hash_table *htable) 
{
	int i;
	for(i=0; i < htable->size; i++)	destroy_chain(&(htable->table[i]));
	free(htable->table);
}
