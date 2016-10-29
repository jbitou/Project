#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hash.h"
#define ITEM_ID 15
#define WINDOW_SIZE 4

int mod(int a, long long b)
{
	if(b<0)
		return mod(a,-b);
	int ret=a%b;
	if(ret<0)
		ret+=b;
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
			g[i][j].t = M + (rand() / (RAND_MAX + 1.0)) * (N-M+1);
			g[i][j].v = NULL;
		}
	}
}

void init_hash_Eucl(ghashp * g, int L, int k, int d)
{
	int i,j,z, M = -1, N = 1, window = WINDOW_SIZE;
	double x,y,w,u,v,mult,w1;
	for (i=0; i < L; i++) //For each table
	{
		for (j=0; j < k; j++)	//Create k times h()
		{ 
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
			g[i][j].t = (rand() / (RAND_MAX + 1.0)) * window;
			g[i][j].r = (rand() / (RAND_MAX + 1.0)) * 128;
		}
	}
}

void init_hash_Cos(ghashp * g, int L, int k, int d)
{
	int i,j,z, M = -1, N = 1;
	double x,y,w,u,v,mult,w1;
	for (i=0; i < L; i++) 	//For each table
	{
		for (j=0; j < k; j++) 	//Create k times h()
		{ 
			g[i][j].v = malloc(d*sizeof(double));	//At cosine: v is r because r is real *
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

void init_hash_matrix(ghashp * g, int **distances, int L, int k, int numofitems) {
	int i, j, y, M = 0, N = numofitems-1, sum, total;
	for (i=0; i < L; i++) 
	{
		for (j=0; j < k; j++)
		{ 
			total = 0;
			g[i][j].t = M + (rand() / (RAND_MAX + 1.0)) * (N-M+1); //x1
			g[i][j].r = M + (rand() / (RAND_MAX + 1.0)) * (N-M+1); //x2
			while(g[i][j].r == g[i][j].t)
				g[i][j].r = M + (rand() / (RAND_MAX + 1.0)) * (N-M+1); //x2
		
			for(y=0; y<numofitems; y++)
			{
				sum = 0;
				//printf("y=%d--x1=%d--x2=%d\n",y,g[i][j].t,g[i][j].r);
				if(y < g[i][j].t)  sum+=pow(distances[y][g[i][j].t-y-1],2);
				else if(y > g[i][j].t) sum+=pow(distances[g[i][j].t][y-g[i][j].t-1],2);
								
				if(y < g[i][j].r)  sum+=pow(distances[y][g[i][j].r-y-1],2);
				else if(y > g[i][j].r) sum+=pow(distances[g[i][j].r][y-g[i][j].r-1],2);
				
				if(g[i][j].t < g[i][j].r)
				{
					sum-=pow(distances[g[i][j].t][g[i][j].r-g[i][j].t-1],2);
					sum = sum / (2*distances[g[i][j].t][g[i][j].r-g[i][j].t-1]);
				}  
				else if(g[i][j].t > g[i][j].r) 
				{
					sum-=pow(distances[g[i][j].r][g[i][j].t-g[i][j].r-1],2);
					sum = sum / (2*distances[g[i][j].r][g[i][j].t-g[i][j].r-1]);
				}
				total += sum;
			}	
			g[i][j].t1 = total/numofitems; //t1
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
	M = (1LL << 32) - 5;
	for (i=0; i < k; i++)		//For each h()
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
	for (i=0; i < k; i++)		//For each h()
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

int hash_func_Matrix(ghashp g, int x, int **distances, int k, int numofitems)
{
	int i,j,h = 0,sum;
	for (i=0; i < k; i++)		
	{
		sum = 0;
		//printf("x=%d--x1=%d--x2=%d\n",x,g[i].t,g[i].r);
		if(x < g[i].t)  sum+=pow(distances[x][g[i].t-x-1],2);
		else if(x > g[i].t) sum+=pow(distances[g[i].t][x-g[i].t-1],2);
		
		if(x < g[i].r)   sum+=pow(distances[x][g[i].r-x-1],2);
		else if(x > g[i].r)  sum+=pow(distances[g[i].r][x-g[i].r-1],2);
		
		if(g[i].t < g[i].r)
		{
			sum-=pow(distances[g[i].t][g[i].r-g[i].t-1],2);
			sum = sum / (2*distances[g[i].t][g[i].r-g[i].t-1]);
		}  
		else if(g[i].t > g[i].r) 
		{
			sum-=pow(distances[g[i].r][g[i].t-g[i].r-1],2);
			sum = sum / (2*distances[g[i].r][g[i].t-g[i].r-1]);
		}
		//printf("sum is %d\n",sum);
		if(sum >= g[i].t1)  h |= 1 << i;
	}	
	return h;
}

int hash_func_MSearch(ghashp g, int *qdata, int **distances, int k, int numofitems)
{
	int i,j,h = 0,sum;
	for (i=0; i < k; i++)		
	{
		sum = 0;
		//printf("--x1=%d--x2=%d\n",g[i].t,g[i].r);
		sum+=pow(qdata[g[i].t],2);
		
		sum+=pow(qdata[g[i].r],2);
		
		if(g[i].t < g[i].r)
		{
			sum-=pow(distances[g[i].t][g[i].r-g[i].t-1],2);
			sum = sum / (2*distances[g[i].t][g[i].r-g[i].t-1]);
		}  
		else if(g[i].t > g[i].r) 
		{
			sum-=pow(distances[g[i].r][g[i].t-g[i].r-1],2);
			sum = sum / (2*distances[g[i].r][g[i].t-g[i].r-1]);
		}
		//printf("sum is %d\n",sum);
		if(sum >= g[i].t1)  h |= 1 << i;
	}	
	return h;
}

void search_table_NNR(int pos, hash_table htable, void *qdata, double rad, nnrp *nnrlist, int flag, int euclID, int d) 
{
	search_chain_NNR(htable.table[pos],qdata,rad,nnrlist,flag,euclID,d);
}

nn search_table_NN(ghashp *g, hash_table *htable, void *qdata, int ** distances, int flag, int d, int k, int L, int tableSize)
{
	int i, counter = 0, bruflag = 0, euclID, pos;
	nn lshnn, lshnn1;
	if (!flag) pos = hash_func_Ham(g[0],qdata,k);
	else if (flag == 1)
	{
		euclID = hash_func_Eucl(g[0],qdata,k,d);
		euclID = abs(euclID);
		pos = mod(euclID, tableSize);
	}
	else if (flag == 2)  pos = hash_func_Cos(g[0],qdata,k,d);	
	else if(flag == 3)	 pos = hash_func_MSearch(g[0],qdata,distances,k,d);
	
	lshnn = search_chain_NN(htable[0].table[pos],qdata,flag,bruflag,euclID,d,&counter,L);
	if (counter > 6*L) 	return lshnn;
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
		lshnn1 = search_chain_NN(htable[i].table[pos],qdata,flag,bruflag,euclID,d,&counter,L);		//Search for NN
		if (((lshnn1.distance < lshnn.distance) && (lshnn1.distance != 0)) || (lshnn.distance == 0))
		{
			lshnn.distance = lshnn1.distance;	
			strcpy(lshnn.key,lshnn1.key);
		}
		if(counter > 6*L) 	return lshnn;
	}
	return lshnn;
}

nn brute_force_table(hash_table htable, void *q, int flag, int bruflag, int euclID, int d, int L)
{
	int i;
	nn tnn, tnn1;
	tnn.key = malloc(ITEM_ID);
	if (!flag)
	{
		char *qdata = (char *)q;
		tnn.distance = 65;		//Initialize minimum distance as something very big
		for (i=0; i<htable.size; i++)
		{
			tnn1 = search_chain_NN(htable.table[i],qdata,flag,bruflag,euclID,d,NULL,L);
			if (((tnn1.distance < tnn.distance) && (tnn1.distance > 0)) ||  (tnn.distance <= 0))
			{
				tnn.distance = tnn1.distance;
				strcpy(tnn.key,tnn1.key);
			}
			free(tnn1.key);
		}
	}
	else if (flag == 3)
	{
		int *qdata = (int *)q;
		tnn.distance = -1;		
		for (i=0; i < htable.size; i++)
		{
			tnn1 = search_chain_NN(htable.table[i],qdata,flag,bruflag,euclID,d,NULL,L);
			if (((tnn1.distance < tnn.distance) && (tnn1.distance > 0)) ||  (tnn.distance <= 0))
			{
				tnn.distance = tnn1.distance;
				strcpy(tnn.key,tnn1.key);
			}
			free(tnn1.key);
		}
	}
	else
	{
		double *qdata = (double *)q;
		int j = 0;
		tnn.distance = -1;
		//strcpy(tnn.key,"no");
		/*while(tnn.distance <= 0) 
		{
			if (flag == 2)	euclID = -1;
			tnn1 = search_chain_NN(htable.table[j],qdata,key,flag,bruflag,euclID,d,NULL,L);		//Keep the first item found as this with minimum distance
			if (((tnn1.distance < tnn.distance) && (tnn1.distance > 0)) ||  (tnn.distance <= 0))
			{
				tnn.distance = tnn1.distance;
				strcpy(tnn.key,tnn1.key);
			}
			j++;
		}*/
		for (i=0; i<htable.size; i++)
		{
			tnn1 = search_chain_NN(htable.table[i],qdata,flag,bruflag,euclID,d,NULL,L);
			if (((tnn1.distance < tnn.distance) && (tnn1.distance > 0)) ||  (tnn.distance <= 0))
			{
				tnn.distance = tnn1.distance;
				strcpy(tnn.key,tnn1.key);
			}
			free(tnn1.key);
		}
	}
	return tnn;
}

void print_table(hash_table htable)
{
	int i;
	for (i=0; i < htable.size; i++)
	{
		printf("Bucket: %d - ",i);
		print_chain(htable.table[i]);
		printf("\n");
	}	
}

void destroy_table(hash_table *htable, int flag) 
{
	int i;
	for(i=0; i < htable->size; i++) 	//For each bucket
	{
		destroy_chain(&(htable->table[i]),flag);
	}
	free(htable->table);
}
