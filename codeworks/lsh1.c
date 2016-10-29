#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>
#include "hash.h"

#define ITEM_ID 10
#define MAX_LINE 3000

int make_item(char *item);

char* deblank(char* input)                                         
{
    int i,j;
    char *output=input;
    for (i = 0, j = 0; i<strlen(input); i++,j++)          
    {
        if (input[i]!=' ')                           
            output[j]=input[i];                     
        else
            j--;                                     
    }
    output[j]=0;
    return output;
}

int main(int argc, char **argv)
{
	FILE *fp, *fq, *fe;
	int i, foundK, foundL, f1, f2, f3, k, L, key, bs, pos, input, query, output;
	char item[ITEM_ID], qitem[ITEM_ID], ms[14], data[65], qdata[65], space[10], radius[8];
	double rad;
	hash_table *htable;
	nnrp nnrlist = NULL;
	clock_t start_t, end_t;
	double total_t, total_t1;
	
	if (argc > 11)
	{
		printf("Too many arguments. Try again.\n");
		return -1;
	}
	if ((argc % 2) == 0)
	{
		printf("Wrong arguments. Try again.\n");	//Every parameter has to be given after a recogniser
		return -1;
	}	
	srand(time(NULL));	//Intializes random number generator
	foundK = foundL = 0;
	f1 = f2 = f3 = 0;
	for (i=1; i<(argc-1); i+=2)
	{
		if (strcmp(argv[i],"-k") == 0)
		{
			foundK = 1;
			k = atoi(argv[i+1]);
			while (k > 10)
			{
				printf("Too big k. Try again: ");
				scanf("%d", &k);
			}  
			while (k <= 0)
			{
				printf("Negative k. Try again giving a positive one: ");
				scanf("%d", &k);
			}  
		}
		else if (strcmp(argv[i],"-L") == 0)
		{
			foundL = 1;
			L = atoi(argv[i+1]);
			while (L > 30)
			{
				printf("Too big L. Try again: ");
				scanf("%d", &L);
			}  
			while (L <= 0)
			{
				printf("Negative L. Try again giving a positive one: ");
				scanf("%d", &L);
			}
		}
		else if (strcmp(argv[i],"-d") == 0)
		{
			f1 = 1;
			input = i+1;
		}
		else if (strcmp(argv[i],"-q") == 0)
		{
			f2 = 1;
			query = i+1;
		}
		else if (strcmp(argv[i],"-o") == 0)
		{
			f3 = 1;
			output = i+1;
		}
	}
	if (!foundK)	k = 4;	//If -k isn't given, use default values
	if (!foundL)	L = 5;	//If -L isn't given, use default values
	
	fp = fopen(argv[input],"r");
	if (fp == NULL)
	{
		perror("Error:");
		return -1;
	}
	fq = fopen(argv[query],"r");
	if (fq == NULL)
	{
		perror("Error:");
		return -1;
	}
	fe = fopen(argv[output],"w+");
	if (fe == NULL)
	{
		perror("Error:");
		return -1;
	}
	
	htable = malloc(L * sizeof(hash_table));
	ghashp *g = malloc(L * sizeof(ghashp));
	for(i = 0; i < L; i++) 
		g[i] = malloc(k * sizeof(ghash));
	fscanf(fp,"%s%s[^\n]",ms,space);	//Read first line of input_file
	printf("metric space=%s\n",space);
	
	int tableSize;
	if (strcmp(space,"hamming") == 0)
	{
		tableSize = pow(2,k);
		for (i=0; i<L; i++)
			init_table(k,&htable[i],tableSize);
		fscanf(fp,"%s %s[^\n]",item,data);	//Read first datasets' line to compute 'N'
		int j;
		init_hash_Ham(g,L,k,data);
		for(i=0; i < L; i++) 
		{
			for(j=0; j < k; j++) 
				printf("%d ",g[i][j].t);		//Print g[i][j] to check
			printf("\n");
		}
		printf("/////\n");
		fseek(fp,0,SEEK_SET);
		fscanf(fp,"%s%s[^\n]",ms,space);
		/*Ιnput phase*/
		while (fscanf(fp,"%s %s[^\n]",item,data)!=EOF)		
		{   
			key = make_item(item);
			for(i = 0; i < L; i++)	
			{
				pos = hash_func_Ham(g[i],data,k);	//Find the right bucket according to function g for Hamming
				insert_chain(key,data,&(htable[i].table[pos]));
			}
		}
		/*End of Ιnput phase*/
		fscanf(fq,"%s%lf[^\n]",radius,&rad);	//Read first line of query_file
		printf("%s	'%lf'\n",radius,rad);
		/*Search phase*/
		while (fscanf(fq,"%s %s[^\n]",qitem,qdata)!=EOF)		
		{
			key = make_item(qitem);
			start_t = clock();
			nn tnn;
			tnn = brute_force_table(htable[0],key,qdata);
			end_t = clock();
			total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			for (i=0; i<L;i++)
			{
				pos = hash_func_Ham(g[i],qdata,k);	//Find the right bucket according to function g
				search_table_NNR(pos,htable[i],qdata,rad,&nnrlist);		//Search for NNRs
			}
			fprintf(fe,"Query: %s\nR-near neighbors:\n",qitem); 
			destroy_nnrlist(&nnrlist,fe);
			start_t = clock();
			nn lshnn, lshnn1;
			pos = hash_func_Ham(g[0],qdata,k);
			lshnn = search_table_NN(pos,htable[0],qdata,key);
			for (i=1; i<L;i++)
			{
				pos = hash_func_Ham(g[i],qdata,k);	//Find the right bucket according to function g
				lshnn1 = search_table_NN(pos,htable[i],qdata,key);		//Search for NN
				if (((lshnn1.distance < lshnn.distance) && (key != lshnn1.key)) || (lshnn.key == key))
					lshnn = lshnn1;	
			}
			end_t = clock();
			total_t1 = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			if (lshnn.key != key)
			{
				fprintf(fe,"Nearest neighbor: item%d\n",lshnn.key);
				fprintf(fe,"True neighbor: item%d\n",tnn.key);
				fprintf(fe,"distanceLSH: %d\n",lshnn.distance);
				fprintf(fe,"distanceTrue: %d\n",tnn.distance);
				fprintf(fe,"tLSH: %f\n",total_t1);
				fprintf(fe,"tTrue: %f\n\n",total_t);
			}
			else
			{
				fprintf(fe,"Nearest neighbor wasn't found\n");
				fprintf(fe,"True neighbor: item%d\n",tnn.key);
				fprintf(fe,"distanceTrue: %d\n",tnn.distance);
				fprintf(fe,"tTrue: %f\n\n",total_t);
			}
		}
		/*End of search phase*/
	}
	else if (strcmp(space,"matrix") == 0)
	{
		//do some work for matrix
	}
	else
	{
		int c;
		int j, l, col = 0, lines = -2, t;
		char metric[20], m[10];
		char eucldata[MAX_LINE], *ptr, delim[] = " \t", *itemID, *point;
		fscanf(fp,"\n");	//Read \n of first line, go to second line
		i = 0;
		while (i < 19)
		{
			c = fgetc(fp);
			if (c == '\n')	break;	//Read only first line
			metric[i++] = c;
			printf("|%c|",c);
		}
		metric[i] = '\0';
		printf("metric: |%s|\n",metric);
		/*Delete blank*/
		for (i = 0, j = 0; i<strlen(metric); i++,j++)          
		{
			if (metric[i] != ' ')                           
				m[j] = metric[i];                     
			else
				j--;                                     
		}
		m[j]=0;
		/*Blank deleted*/
		printf("real metric: |%s|\n",m);
		fgets(eucldata,MAX_LINE,fp);	//Read real first line of data
		ptr = eucldata;
		itemID = strtok(ptr,delim);		//Get first item_id to count dimensions
		ptr = NULL;
		/*Count dimensions*/		
		while ((point = strtok(ptr,delim)) != NULL)   
		{
			if (lines == -2)	col++;
			ptr = NULL; 
		}
		/*Dimensions found*/
		printf("item: |%s|\n",itemID);
		printf("dimensions: %d\n",col);
		fseek(fp,0,SEEK_SET);	//Return file to start
		while (fgets(eucldata,MAX_LINE,fp) != NULL) 	lines++;	//Count lines 
		printf("lines: %d\n",lines);
		int power = pow(2,k);
		if ((lines/2) > power)	tableSize = power;
		else 	tableSize = lines/2;
		for (i=0; i < L; i++)	//Allocate memory for tables
			init_table(k,&htable[i],tableSize);
		if (strcmp(m,"@metriceuclidean") == 0)
		{
			init_hash_Eucl(g,L,k,col);
		}
		double p[col];		//Allocate array of double to store points
		/*Input phase*/
        fseek(fp,0,SEEK_SET);
        fgets(eucldata,MAX_LINE,fp);	//Read first line of input_file
		fgets(eucldata,MAX_LINE,fp);		//Read second line 
		while (fgets(eucldata,MAX_LINE,fp) != NULL)		
		{
			ptr = eucldata;
			itemID = strtok(ptr,delim);
			l = 0;	
			ptr = NULL;		
			while ((point = strtok(ptr,delim)) != NULL)   
			{
				printf("|%s|\n",point);
				p[l] = atof(point);
				ptr = NULL; 
				l++;  
			}
			if (strcmp(m,"@metriceuclidean") == 0)
			{
				for(i = 0; i < L; i++)	
				{
					pos = hash_func_Eucl(g[i],p,k,col);
					printf("pos=%d\n",pos);
					pos = abs(pos);
					printf("pos=%d\n",pos);
					int logn;
					printf("lines:%d\n", lines);
					logn = (int)log(lines);
					printf("log: %d\n",logn);	
					pos = pos % (int)pow(2,k);
					printf("pos=%d\n",pos);
				}				
			}
		}
		/*End of Input phase*/
	}
	for(i = 0; i < L; i++)
		free(g[i]);
	free(g);
	for (i=0; i<L; i++)
		destroy_table(&htable[i]);	
	free(htable);
	fclose(fe);
	fclose(fp);	
	fclose(fq);	
		
	return 0;
}

int make_item(char *item)
{
	int key;
	if (!isdigit(item[0]))	//If the first character of string is type of char
	{	
		char *id;
		int s = strlen(item) - strlen("item");
		id = malloc(s+1);
		strncpy(id,item+strlen("item"),s);
		key = atoi(id);		//keep only K of item_idK
		free(id);
	}
	else 	key = atoi(item);		//else, just convert type char* to int
	return key;
}
