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

void init_hash_matrix(ghashp * g, int **distances, int L, int k, int numofitems) 
{
	int i, j, y, M = 0, N = numofitems-1, sum, total;
	/*For each table*/
	for (i=0; i < L; i++) 
	{
		/*Create k times h()*/
		for (j=0; j < k; j++)
		{ 
			total = 0;
			g[i][j].t = M + (rand() / (RAND_MAX + 1.0)) * (N-M+1); //x1
			g[i][j].r = M + (rand() / (RAND_MAX + 1.0)) * (N-M+1); //x2
			while(g[i][j].r == g[i][j].t)
				g[i][j].r = M + (rand() / (RAND_MAX + 1.0)) * (N-M+1); //x2
			/*Calculate t1*/
			for(y=0; y < numofitems; y++)
			{
				sum = 0;
				if(y < g[i][j].t)  sum += pow(distances[y][g[i][j].t-y-1],2);
				else if(y > g[i][j].t) sum += pow(distances[g[i][j].t][y-g[i][j].t-1],2);
								
				if(y < g[i][j].r)  sum += pow(distances[y][g[i][j].r-y-1],2);
				else if(y > g[i][j].r) sum += pow(distances[g[i][j].r][y-g[i][j].r-1],2);
				
				if(g[i][j].t < g[i][j].r)
				{
					sum -= pow(distances[g[i][j].t][g[i][j].r-g[i][j].t-1],2);
					sum = sum / (2*distances[g[i][j].t][g[i][j].r-g[i][j].t-1]);
				}  
				else if(g[i][j].t > g[i][j].r) 
				{
					sum -= pow(distances[g[i][j].r][g[i][j].t-g[i][j].r-1],2);
					sum = sum / (2*distances[g[i][j].r][g[i][j].t-g[i][j].r-1]);
				}
				total += sum;
			}	
			g[i][j].t1 = total / numofitems; //t1
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

int hash_func_Eucl(ghashp g, double *p, int k, int d)
{
	int i = 0, j, window = WINDOW_SIZE, h, rh, sum = 0;
	long long M;
	double inner;
	M = (1LL << 32) - 5;	//M: known prime
	/*For each h()*/
	for (i=0; i < k; i++)	
	{
		inner = 0.0;
		for (j=0; j < d; j++)	//Inner product
			inner += g[i].v[j]*p[j];
		inner += g[i].t;
		inner /= window;
		h = (int)inner;
		rh = h*g[i].r;
		sum += rh;
	}
	sum = mod(sum,M);
	return sum;		//return ID(p)
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

/*Hash function for insert*/
int hash_func_Matrix(ghashp g, int x, int **distances, int k, int numofitems)
{
	int i,h = 0,sum;
	char *end;
	char *temp = malloc(k*sizeof(char));
	for (i=0; i < k; i++)		
	{
		sum = 0;
		if(x < g[i].t)  sum += pow(distances[x][g[i].t-x-1],2);
		else if(x > g[i].t) sum += pow(distances[g[i].t][x-g[i].t-1],2);
		
		if(x < g[i].r)   sum += pow(distances[x][g[i].r-x-1],2);
		else if(x > g[i].r)  sum += pow(distances[g[i].r][x-g[i].r-1],2);
		
		if(g[i].t < g[i].r)
		{
			sum -= pow(distances[g[i].t][g[i].r-g[i].t-1],2);
			sum = sum / (2*distances[g[i].t][g[i].r-g[i].t-1]);
		}  
		else if(g[i].t > g[i].r) 
		{
			sum -= pow(distances[g[i].r][g[i].t-g[i].r-1],2);
			sum = sum / (2*distances[g[i].r][g[i].t-g[i].r-1]);
		}
		if(sum >= g[i].t1)  temp[i]='1';
		else temp[i]='0';
	}	
	h = strtol(temp,&end,2);
	free(temp);
	return h;
}

/*Hash function for search*/
int hash_func_MSearch(ghashp g, int *qdata, int **distances, int k, int numofitems)
{
	int i,h = 0,sum;
	char *end;
	char *temp = malloc(k*sizeof(char));
	for (i=0; i < k; i++)		
	{
		sum = 0;
		sum += pow(qdata[g[i].t],2);
		
		sum += pow(qdata[g[i].r],2);
		
		if(g[i].t < g[i].r)
		{
			sum -= pow(distances[g[i].t][g[i].r-g[i].t-1],2);
			sum = sum / (2*distances[g[i].t][g[i].r-g[i].t-1]);
		}  
		else if(g[i].t > g[i].r) 
		{
			sum -= pow(distances[g[i].r][g[i].t-g[i].r-1],2);
			sum = sum / (2*distances[g[i].r][g[i].t-g[i].r-1]);
		}
		if(sum >= g[i].t1)  temp[i]='1';
		else temp[i]='0';
	}	
	h = strtol(temp,&end,2);
	free(temp);
	return h;
}

void search_table_NNR(int pos, hash_table htable, void *qdata, double rad, nnrp *nnrlist, int flag, int euclID, int d) 
{
	search_chain_NNR(htable.table[pos],qdata,rad,nnrlist,flag,euclID,d);
}

void search_table_NN(ghashp *g, hash_table *htable, void *qdata, int ** distances, int flag, int d, int k, int L, int tableSize, nnrp *list, double *diff)
{
	int i, counter = 0, bruflag = 0, euclID, pos;
	nnrp lshnn = NULL;
	double diff1;
	/*Find position in hash table according to metric*/
	if (!flag) pos = hash_func_Ham(g[0],qdata,k);
	else if (flag == 1)
	{
		euclID = hash_func_Eucl(g[0],qdata,k,d);
		euclID = abs(euclID);
		pos = mod(euclID, tableSize);
	}
	else if (flag == 2)  pos = hash_func_Cos(g[0],qdata,k,d);	
	else if(flag == 3)	 pos = hash_func_MSearch(g[0],qdata,distances,k,d);	
	/*Bring the first as minimum*/
	search_chain_NN(htable[0].table[pos],qdata,flag,bruflag,euclID,d,&counter,L,list,diff);	
	/*Heuristic choice for trick, Euclidean and Cosine*/
	if ((flag == 1) || (flag == 2))
	{
		if (counter > TRICK*L) 	return;
	}
	for (i=1; i < L; i++)
	{
		if (!flag) pos = hash_func_Ham(g[i],qdata,k);
		else if (flag == 1)
		{
			euclID = hash_func_Eucl(g[i],qdata,k,d);
			euclID = abs(euclID);
			pos = mod(euclID, tableSize);
		}
		else if (flag == 2)  pos = hash_func_Cos(g[i],qdata,k,d);
		else if(flag == 3)	 pos = hash_func_MSearch(g[i],qdata,distances,k,d);	
		search_chain_NN(htable[i].table[pos],qdata,flag,bruflag,euclID,d,&counter,L,&lshnn,&diff1);		
		/*Compare to get the new min*/
		if (((diff1 < *diff) && (diff1 > 0)) || (*diff == 0))
		{
			*diff = diff1;
			destroy_nnrlist(list);
			*list = lshnn;
			lshnn = NULL;
		}
		else if((diff1 == *diff) && (*diff > 0)) combine_nnrlist(list,&lshnn);
		else destroy_nnrlist(&lshnn);
		/*Heuristic choice for trick, Euclidean and Cosine*/
		if ((flag == 1) || (flag == 2))
		{
			if (counter > TRICK*L) 	return;
		}
	}
	return;
}

void brute_force_table(hash_table htable, void *q, int flag, int bruflag, int euclID, int d, int L, nnrp *list,double *diff)
{
	int i;
	nnrp tlist = NULL;
	double diff1;
	if (!flag)
	{
		char *qdata = (char *)q;
		/*Initialize minimum distance as something very big*/
		*diff= 65;		
		/*For each bucket*/
		for (i=0; i < htable.size; i++)
		{
			search_chain_NN(htable.table[i],qdata,flag,bruflag,euclID,d,NULL,L,&tlist,&diff1);
			if ((diff1 < *diff) && (diff1 > 0))
			{
				*diff = diff1;	
				destroy_nnrlist(list);
				*list = tlist;
				tlist = NULL;
			}
			else if ((diff1 == *diff) && (*diff > 0)) 	combine_nnrlist(list,&tlist);
			else	destroy_nnrlist(&tlist);
		}
	}
	else if (flag == 3)
	{
		int *qdata = (int *)q;
		/*Initialize minimum distance as something invalid*/
		*diff = -1;		
		/*For each bucket*/
		for (i=0; i < htable.size; i++)
		{
			search_chain_NN(htable.table[i],qdata,flag,bruflag,euclID,d,NULL,L,&tlist,&diff1);
			if (((diff1 < *diff) && (diff1 > 0)) ||  (*diff <= 0))
			{
				*diff = diff1;	
				destroy_nnrlist(list);
				*list = tlist;
				tlist = NULL;
			}
			else if ((diff1 == *diff) && (*diff > 0)) 	combine_nnrlist(list,&tlist);
			else	destroy_nnrlist(&tlist);
		}
	}
	else
	{
		double *qdata = (double *)q;
		/*Initialize minimum distance as something invalid*/
		*diff = -1;		
		/*For each bucket*/
		for (i=0; i < htable.size; i++)
		{
			search_chain_NN(htable.table[i],qdata,flag,bruflag,euclID,d,NULL,L,&tlist,&diff1);
			if (((diff1 < *diff) && (diff1 > 0)) ||  (*diff <= 0))
			{
				*diff = diff1;	
				destroy_nnrlist(list);
				*list = tlist;
				tlist = NULL;
			}
			else if ((diff1 == *diff) && (*diff > 0)) 	combine_nnrlist(list,&tlist);
			else	destroy_nnrlist(&tlist);
		}
	}
	return;
}

void destroy_table(hash_table *htable, int flag) 
{
	int i;
	/*For each bucket*/
	for(i=0; i < htable->size; i++) 	
	{
		destroy_chain(&(htable->table[i]),flag);
	}
	free(htable->table);
}
