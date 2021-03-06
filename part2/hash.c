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

int *sortMatrix(int *array, int N) {
	int i, j, a;
	for (i = 0; i < N; i++) {
		for (j = i + 1; j < N; j++) {
			if (array[i] > array[j]) {
				a =  array[i];
                array[i] = array[j];
                array[j] = a;
            }
        }
	}
	return array;
}

void init_table(int k, hash_table *htable, int tableSize) 
{
	int i;
	htable->size = tableSize;
	htable->table = malloc((htable->size)*sizeof(chainp));
	for(i=0; i < htable->size; i++)	htable->table[i] = NULL;
	return;
}

void init_hash_Ham(ghashp *g, int L, int k, char *data)
{
	int i, j, M = 1, N = strlen(data);
	for (i=0; i < L; i++) {
		for (j=0; j < k; j++) { 
			/**Choose uniformly an h function (position of a bit)**/
			g[i][j].t = M + (rand() / (RAND_MAX + 1.0)) * (N-M+1);
			g[i][j].v = NULL;
		}
	}
}

void init_hash_Eucl(ghashp * g, int L, int k, int d)
{
	int i,j,z, M = -1, N = 1, window = WINDOW_SIZE;
	double x,w,u,v,mult;
	/**For each table**/
	for (i=0; i < L; i++) {
		/**Create k times h()**/
		for (j=0; j < k; j++) { 
			g[i][j].v = malloc(d*sizeof(double));
			/**Create a vector~N(0,1)**/
			for (z=0; z < d; z++) { 
				do {
					u = M + (rand() / (RAND_MAX + 1.0)) * (N-M);
					v = M + (rand() / (RAND_MAX + 1.0)) * (N-M);
					w = pow (u,2) + pow (v,2);
				}while (w >= 1 || w == 0);
				mult = sqrt ((-2 * log (w)) / w);
				x = u * mult;
				g[i][j].v[z] = x;
			}
			g[i][j].t = (rand() / (RAND_MAX + 1.0)) * window;
			/**Choose r in range 0-128**/
			g[i][j].r = (rand() / (RAND_MAX + 1.0)) * 128;
		}
	}
}

void init_hash_Cos(ghashp * g, int L, int k, int d)
{
	int i,j,z, M = -1, N = 1;
	double x,w,u,v,mult;
	/**For each table**/
	for (i=0; i < L; i++) {
		/**Create k times h()**/
		for (j=0; j < k; j++) { 
			/**At cosine: v is r because r is real**/
			g[i][j].v = malloc(d*sizeof(double));	
			for (z=0; z < d; z++) { 
				do {
					u = M + (rand() / (RAND_MAX + 1.0)) * (N-M);
					v = M + (rand() / (RAND_MAX + 1.0)) * (N-M);
					w = pow (u,2) + pow (v,2);
				}while (w >= 1 || w == 0);
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
	/**For each table**/
	for (i=0; i < L; i++) {
		/**Create k times h()**/
		for (j=0; j < k; j++) { 
			total = 0;
			g[i][j].t = M + (rand() / (RAND_MAX + 1.0)) * (N-M+1); //x1
			g[i][j].r = M + (rand() / (RAND_MAX + 1.0)) * (N-M+1); //x2
			while(g[i][j].r == g[i][j].t) g[i][j].r = M + (rand() / (RAND_MAX + 1.0)) * (N-M+1); //x2
			/**Calculate t1**/
			for(y=0; y < numofitems; y++) {
				sum = 0;
				if(y < g[i][j].t)  sum += pow(distances[y][g[i][j].t-y-1],2);
				else if(y > g[i][j].t) sum += pow(distances[g[i][j].t][y-g[i][j].t-1],2);							
				if(y < g[i][j].r)  sum += pow(distances[y][g[i][j].r-y-1],2);
				else if(y > g[i][j].r) sum += pow(distances[g[i][j].r][y-g[i][j].r-1],2);				
				if(g[i][j].t < g[i][j].r) {
					sum -= pow(distances[g[i][j].t][g[i][j].r-g[i][j].t-1],2);
					sum = sum / (2*distances[g[i][j].t][g[i][j].r-g[i][j].t-1]);
				}  
				else if(g[i][j].t > g[i][j].r) {
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
	char *temp = malloc((k+1)*sizeof(char));
	for (i=0; i < k; i++)	temp[i] = data[g[i].t];
	temp[i] = '\0';	
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
	/**For each h()**/
	for (i=0; i < k; i++)	
	{
		inner = 0.0;
		/**Inner product**/
		for (j=0; j < d; j++)	inner += g[i].v[j]*p[j];
		inner += g[i].t;
		inner /= window;
		h = (int)inner;
		rh = h*g[i].r;
		sum += rh;
	}
	sum = mod(sum,M);
	/**return ID(p)**/
	return sum;	
}

int hash_func_Cos(ghashp g, double *x, int k, int d)
{
	double inner;
	int i,j,h;
	char *end;
	char *temp = malloc((k+1)*sizeof(char));
	/**For each h()**/
	for (i=0; i < k; i++) {
		inner = 0.0;
		/**Inner product**/
		for (j=0; j < d; j++)	
			inner += g[i].v[j]*x[j];	//ri*x
		if (inner >= 0.0)
			temp[i] = '1';
		else
			temp[i] = '0';	
	}
	temp[i] = '\0';
	h = strtol(temp,&end,2);
	free(temp);
	return h;
}

/**Hash function for insert**/
int hash_func_Matrix(ghashp g, int x, int **distances, int k, int numofitems)
{
	int i,h = 0,sum;
	char *end;
	char *temp = malloc((k+1)*sizeof(char));
	for (i=0; i < k; i++) {
		sum = 0;
		if(x < g[i].t)  sum += pow(distances[x][g[i].t-x-1],2);
		else if(x > g[i].t) sum += pow(distances[g[i].t][x-g[i].t-1],2);		
		if(x < g[i].r)   sum += pow(distances[x][g[i].r-x-1],2);
		else if(x > g[i].r)  sum += pow(distances[g[i].r][x-g[i].r-1],2);		
		if(g[i].t < g[i].r) {
			sum -= pow(distances[g[i].t][g[i].r-g[i].t-1],2);
			sum = sum / (2*distances[g[i].t][g[i].r-g[i].t-1]);
		}  
		else if(g[i].t > g[i].r) {
			sum -= pow(distances[g[i].r][g[i].t-g[i].r-1],2);
			sum = sum / (2*distances[g[i].r][g[i].t-g[i].r-1]);
		}
		if(sum >= g[i].t1)  temp[i]='1';
		else temp[i]='0';
	}	
	temp[i] = '\0';
	h = strtol(temp,&end,2);
	free(temp);
	return h;
}

/**Hash function for search**/
int hash_func_MSearch(ghashp g, int *qdata, int **distances, int k, int numofitems)
{
	int i,h = 0,sum;
	char *end;
	char *temp = malloc((k+1)*sizeof(char));
	for (i=0; i < k; i++) {
		sum = 0;
		sum += pow(qdata[g[i].t],2);
		
		sum += pow(qdata[g[i].r],2);
		
		if(g[i].t < g[i].r) {
			sum -= pow(distances[g[i].t][g[i].r-g[i].t-1],2);
			sum = sum / (2*distances[g[i].t][g[i].r-g[i].t-1]);
		}  
		else if(g[i].t > g[i].r) {
			sum -= pow(distances[g[i].r][g[i].t-g[i].r-1],2);
			sum = sum / (2*distances[g[i].r][g[i].t-g[i].r-1]);
		}
		if(sum >= g[i].t1)  temp[i]='1';
		else temp[i]='0';
	}	
	temp[i] ='\0';
	h = strtol(temp,&end,2);
	free(temp);
	return h;
}

int **create_small_distance_matrix(int **distances, int n, int N, int *insert) {
	int **small, i, j, z, y, pos;
	small = malloc((n-1)*sizeof(int *));	
	j = n - 1;
	for(i=0; i < (n - 1); i++)  {
		small[i] = malloc(j*sizeof(int));
		j--;
	}
	insert = sortMatrix(insert,n);
	for(i=0; i < (n-1); i++)  {
		pos = 0;
		for (j=0; j < (N - 1); j++) {
			/**Found item in big distance table(line)**/
			if (insert[i] == j) {
				/**For each item in sample find its distance with insert[i]**/
				for (z=0; z < n; z++) {
					for (y=0; y < (N - 1); y++) {
						if (insert[z] == y)	{
							if (j < y) {
								small[i][pos] = distances[j][y-j-1];
								pos++;
							}
						}
					}
				}
			}
		}
	}
	return small;
}

int **create_ham_distance_table(hash_table htable, int N) {
	chainp temp, other;
	int i, j, id1, id2, **distances;
	/**Allocate memory for vector matrix**/
	distances = malloc((N-1)*sizeof(int *));	
	j = N - 1;
	for(i=0; i < (N - 1); i++)  {
		distances[i] = malloc(j*sizeof(int));
		j--;
	}
	/**Compute distances for each item inside hash table 0**/
	for(i=0; i < htable.size; i++)  {
		temp = htable.table[i];
		while (temp != NULL) {
			id1 = temp->position;
			/**There is no line for last item**/
			if (id1 == (N-1)) {
				temp = temp->next;
				continue;
			}
			for(j=0; j < htable.size; j++)  {
				other = htable.table[j];
				while (other != NULL) {
					id2 = other->position;
					if (id1 < id2) distances[id1][id2-id1-1] = distance_Hamming(*temp->value,*other->value);
					other = other->next;
				}
			}
			temp = temp->next;
		}
	}
	return distances;
}

double **create_vector_distance_table(hash_table htable, int N, int d, int flag) {
	double **distances;
	chainp temp, other;
	int i, j, id1, id2;
	/**Allocate memory for vector matrix**/
	distances = malloc((N-1)*sizeof(double *));	
	j = N - 1;
	for (i=0; i < (N - 1); i++)  {
		distances[i] = malloc(j*sizeof(double));
		j--;
	}
	/**Compute distances for each item inside hash table 0**/
	for (i=0; i < htable.size; i++)  {
		temp = htable.table[i];
		while (temp != NULL) {
			id1 = temp->position;
			/**There is no line for last item**/
			if (id1 == (N-1)) {
				temp = temp->next;
				continue;
			}
			for (j=0; j < htable.size; j++)  {
				other = htable.table[j];
				while (other != NULL) {
					id2 = other->position;
					if (id1 < id2) {
						if (flag == 1)	distances[id1][id2-id1-1] = distance_Euclidean(temp->p,other->p,d);
						else 			distances[id1][id2-id1-1] = distance_Cosine(temp->p,other->p,d);
					}
					other = other->next;
				}
			}
			temp = temp->next;
		}
	}
	return distances;
}

int search_table_NNR(int pos, hash_table* htable, void *center, double rad, pointp *list, chainp *barriers, int flag, int euclID, int d, int *all) 
{
	return search_chain_NNR(&(htable->table[pos]),center,rad,list,&(barriers[pos]),flag,euclID,d,all);
}

double *find_vector_info(hash_table htable, int id) {
	int i;
	chainp temp;
	for (i=0; i < htable.size; i++) {
		temp = htable.table[i];
		while (temp != NULL) {
			if (id == temp->position) return temp->p;
			temp = temp->next;
		}
	}
}

char *find_ham_info(hash_table htable, int id) {
	int i, c;
	char *array;
	chainp temp;
	array = malloc(65);
	for (i=0; i < htable.size; i++) {
		temp = htable.table[i];
		while (temp != NULL) {
			if (id == temp->position) {
				for (i = 0; i < 64; i++) {
					c = (*temp->value) & (1LL << i) ? 1 : 0;
					if (!c)	array[63-i] = '0';
					else 	array[63-i] = '1';
				}
				array[64] = '\0';
				return array;
			}
			temp = temp->next;
		}
	}
}

void destroy_table(hash_table *htable, int flag) 
{
	int i;
	/**For each bucket**/
	for(i=0; i < htable->size; i++) 	destroy_chain(&(htable->table[i]),flag);
	free(htable->table);
}
