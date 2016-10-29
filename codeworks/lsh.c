#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include "hash.h"
#define ITEM_ID 15
#define MAX_LINE 1000

char *inputString(FILE *, size_t);

int main(int argc, char **argv)
{
	FILE *fp, *fq, *fe;
	int i, j, foundK, foundL, f1, f2, f3, k, L, key, bs, pos, input, query, output, flag, valid, tableSize, bruflag, files;
	char item[ITEM_ID], qitem[ITEM_ID], ms[14], data[65], qdata[65], space[10], radius[8], answer[2], fileName[50];
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
	fe = fopen(argv[output],"w+");
	if (fe == NULL)
	{
		perror("Error:");
		return -1;
	}
	
	htable = malloc(L * sizeof(hash_table));
	ghashp *g = malloc(L * sizeof(ghashp));		//Allocate memory for hash family G
	for(i = 0; i < L; i++) 
		g[i] = malloc(k * sizeof(ghash));
	fscanf(fp,"%s%s[^\n]",ms,space);	//Read first line of input_file
	if (strcmp(space,"hamming") == 0)
	{
		flag = 0;
		tableSize = 1 << k;
		for (i=0; i<L; i++)
			init_table(k,&htable[i],tableSize);
		fscanf(fp,"%s %s[^\n]",item,data);	//Read first datasets' line to compute 'N'
		init_hash_Ham(g,L,k,data);
		/*for(i=0; i < L; i++) 
		{
			for(j=0; j < k; j++) 
				printf("%d ",g[i][j].t);		//Print g[i][j] to check
			printf("\n");
		}*/
		fseek(fp,0,SEEK_SET);
		fscanf(fp,"%s%s[^\n]",ms,space);
		/*Ιnput phase*/
		while (fscanf(fp,"%s %s[^\n]",item,data)!=EOF)		
		{   
			//key = make_item(item);
			for(i = 0; i < L; i++)	
			{
				pos = hash_func_Ham(g[i],data,k);	//Find the right bucket according to function g for Hamming
				insert_chain(item,data,&(htable[i].table[pos]),flag,0,0);
			}
		}
		/*End of Ιnput phase*/
		/*for(i = 0; i < L; i++)	
		{
			printf("TABLE %d\n",i);
			print_table(htable[i]);
		}*/
		files = 1;
		do
		{
			if (files == 1) 	fq = fopen(argv[query],"r");
			else
			{
				printf("Give QueryFile name: ");
				scanf("%s",fileName);
				fq = fopen(fileName,"r");
				printf("Give OutputFile name: ");
				scanf("%s",fileName);
				fe = fopen(fileName,"w+");
			}
			if (fq == NULL)
			{
				perror("Error:");
				return -1;
			}
			fscanf(fq,"%s%lf[^\n]",radius,&rad);	//Read first line of query_file
			if (rad > 0)	valid = 1;		//If a positive Radius is given, then compute all neighbors
			else if (rad == 0) valid = 0;	//If Radius equals to zero, then find only THE nearest neighbor
			else 							//If Radius < 0, then report error
			{
				printf("Negative radius found. Try again\n");
				return -1;
			}
			/*Search phase*/
			while (fscanf(fq,"%s %s[^\n]",qitem,qdata)!=EOF)		
			{
				bruflag = 1;
				start_t = clock();
				nn tnn;
				tnn = brute_force_table(htable[0],qdata,flag,bruflag,0,0,L);
				end_t = clock();
				total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
				bruflag = 0;
				if (valid)
				{
					for (i=0; i < L;i++)
					{
						pos = hash_func_Ham(g[i],qdata,k);	//Find the right bucket according to function g
						search_table_NNR(pos,htable[i],qdata,rad,&nnrlist,flag,0,0);		//Search for NNRs
					}
					fprintf(fe,"Query: %s\nR-near neighbors:\n",qitem); 
					destroy_nnrlist(&nnrlist,fe);
				}
				start_t = clock();
				nn lshnn;			
				lshnn = search_table_NN(g,htable,qdata,NULL,flag,0,k,L,tableSize);
				end_t = clock();
				total_t1 = (double)(end_t - start_t) / CLOCKS_PER_SEC;
				fprintf(fe,"Nearest neighbor: %s\n",lshnn.key);
				fprintf(fe,"True neighbor: %s\n",tnn.key);
				fprintf(fe,"distanceLSH: %f\n",lshnn.distance);
				fprintf(fe,"distanceTrue: %f\n",tnn.distance);
				fprintf(fe,"tLSH: %f\n",total_t1);
				fprintf(fe,"tTrue: %f\n\n",total_t);
				free(lshnn.key);
				free(tnn.key);
			}	/*End of search phase*/
			printf("Search completed with success.\n");
			fclose(fq);	
			fclose(fe);
			printf("Would you like to continue the search? Y or N: ");
			scanf("%s",answer);
			files++;
		}while (strcmp(answer,"Y") == 0);
	}
	else if (strcmp(space,"matrix") == 0)
	{
		char itemsline[7],*allitems;
		int numofitems, token, itemid;
		tableSize = 1 << k;
		flag = 3;	
		fscanf(fp,"%s",itemsline);
		allitems = inputString(fp,MAX_LINE);
		/*read line with items to get the size of te matrix (numofitems x numofitems)*/
		i = 0;
		numofitems = 1;
		while( allitems[i] != '\0')		
		{
			if( allitems[i]== ',')
				numofitems++;
			i++;	
		}
		printf("numofitems=%d\n",numofitems);
		/*Allocate memory for tables*/
		for (i=0; i < L; i++)	
			init_table(k,&htable[i],tableSize);
		/*read table*/
		int **p = malloc((numofitems-1)*sizeof(int*));	
		j = numofitems-1;
		for(i=0; i < numofitems-1; i++) 
		{
			p[i] = malloc(j*sizeof(int));
			j--;
		}
		int flag1,z,y;
		i=0;
		while(fscanf(fp,"%d",&token) != EOF)		
		{
			flag1=0;
			z=0;
			if (token == 0)  flag1 = 1;	//From now and then, read next distance and store it
			for(j=0; j < numofitems-1; j++)
			{
				fscanf(fp,"%d",&token);
				if(flag1)	
				{
					p[i][z] = token;
					z++;	
				}	
				if (token == 0)  flag1=1;
			}
			i++;	
		}
		/*z = numofitems-1;
		for(i=0; i < numofitems-1; i++) {
			for(j=0; j < z; j++)
				printf("distance of %d from %d is %d\n",i,i+j+1,p[i][j]);
			z--;
		}*/
		/*create g functions*/	
		init_hash_matrix(g,p,L,k,numofitems);
		/*for (i=0; i < L; i++) 
		{
			for (j=0; j < k; j++)
			{ 
				printf("H%d,%d--x1=%d--x2=%d--t1=%d\n",i,j,g[i][j].t,g[i][j].r,g[i][j].t1);
			}
		}*/
		char itemID[ITEM_ID];
		/*start input*/
		for(i = 0; i < L; i++)	
		{ 
			//z = numofitems-1;
			for(j = 0; j < numofitems; j++)	
			{
				sprintf(itemID,"item%d",j+1);
				pos = hash_func_Matrix(g[i],j,p,k,numofitems);
				insert_chain(itemID,NULL,&(htable[i].table[pos]),flag,0,0);
			}
		}
		/*end of input phase*/
		/*for(i = 0; i < L; i++)	
		{
			printf("TABLE %d\n",i);
			print_table(htable[i]);	
		}*/
		files = 1;
		do
		{
			if (files == 1) 	fq = fopen(argv[query],"r");
			else
			{
				printf("Give QueryFile name: ");
				scanf("%s",fileName);
				fq = fopen(fileName,"r");
				printf("Give OutputFile name: ");
				scanf("%s",fileName);
				fe = fopen(fileName,"w+");
			}
			if (fq == NULL)
			{
				perror("Error:");
				return -1;
			}
			/*read file, get radius*/
			fscanf(fq,"%s%lf[^\n]",radius,&rad);
			if (rad > 0)	valid = 1;		
			else if (rad == 0) valid = 0;
			else 							
			{
				printf("Negative radius found. Try again\n");
				return -1;
			}
			/*search starts*/
			int *qdata = malloc(numofitems*sizeof(int));
			while(fscanf(fq,"%s",item) != EOF)		
			{
				for(i=0; i < numofitems; i++)
				{
					fscanf(fq,"%d",&token);
					qdata[i] = token;
				}
				if (valid)
				{
					for (i=0; i < L;i++)
					{
						pos = hash_func_MSearch(g[i],qdata,p,k,numofitems);
						search_table_NNR(pos,htable[i],qdata,rad,&nnrlist,flag,0,0);	
					}
					fprintf(fe,"Query: %s\nR-near neighbors:\n",item); 
					destroy_nnrlist(&nnrlist,fe);
				}
				/*brute starts*/
				bruflag = 1;
				start_t = clock();
				nn tnn;
				tnn = brute_force_table(htable[0],qdata,flag,bruflag,0,0,L);
				end_t = clock();
				total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
				bruflag = 0;
				/*brute ends*/
				/*NN starts*/
				start_t = clock();
				nn lshnn;		
				lshnn = search_table_NN(g,htable,qdata,p,flag,numofitems,k,L,tableSize);	
				end_t = clock();
				total_t1 = (double)(end_t - start_t) / CLOCKS_PER_SEC;
				fprintf(fe,"Nearest neighbor: %s\n",lshnn.key);
				fprintf(fe,"True neighbor: %s\n",tnn.key);
				fprintf(fe,"distanceLSH: %f\n",lshnn.distance);
				fprintf(fe,"distanceTrue: %f\n",tnn.distance);
				fprintf(fe,"tLSH: %f\n",total_t1);
				fprintf(fe,"tTrue: %f\n\n",total_t);
				free(lshnn.key);
				free(tnn.key);
			}
			free(qdata);
			printf("Search completed with success.\n");
			fclose(fq);	
			fclose(fe);
			printf("Would you like to continue the search? Y or N: ");
			scanf("%s",answer);
			files++;
		}while (strcmp(answer,"Y") == 0);
		free(allitems);
		for(i=0; i < numofitems-1; i++)	
			free(p[i]);
		free(p);
		
	}
	else
	{
		int l, col, lines = -2, c, t, itemid, euclID;
		char metric[20], m[10], token[100];
		char *eucldata, *ptr, delim[] = " \t", itemID[ITEM_ID], *point;
		fscanf(fp,"\n");	//Read \n of first line, go to second line
		i = 0;
		while (i < 19)
		{
			c = fgetc(fp);
			if (c == '\n')	break;	//Read only second line
			metric[i++] = c;
			printf("|%c|",c);
		}
		metric[i] = '\0';
		/*Delete blank*/
		for (i = 0, j = 0; i < strlen(metric); i++,j++)          
		{
			if (metric[i] != ' ')                           
				m[j] = metric[i];                     
			else
				j--;                                     
		}
		m[j]=0;
		/*Blank deleted*/
		fscanf(fp,"%s",itemID);
		eucldata = inputString(fp,MAX_LINE);
		/*Count dimensions*/
		ptr = eucldata;
		strtok(ptr,delim);		//Get first item_id to count dimensions
		ptr = NULL;
		col = 1;		
		while ((point = strtok(ptr,delim)) != NULL)   
		{
			if (lines == -2)	col++;	//Count dimensions of vectors
			ptr = NULL; 
		}
		/*Dimensions found*/
		printf("dimensions: %d\n",col);
		fseek(fp,0,SEEK_SET);	//Return file to start
		/*count lines*/
		int ch;
		while ((ch = fgetc(fp)) != EOF) 
		{
			if(ch == '\n') lines++;
		}
		printf("lines: %d\n",lines);
		/*lines found*/
		if (strcmp(m,"@metriceuclidean") == 0)
		{
			tableSize = lines/16 + 1;
			flag = 1;
		}
		else 
		{
			tableSize = 1 << k;
			flag = 2;
		}		
		for (i=0; i < L; i++)	//Allocate memory for tables
			init_table(k,&htable[i],tableSize);
		if (strcmp(m,"@metriceuclidean") == 0)	
		{
			init_hash_Eucl(g,L,k,col);	//Euclidean metric
		}
		else 	
		{
			init_hash_Cos(g,L,k,col);	//Cosine metric
		}
		/*Input phase*/
        fseek(fp,0,SEEK_SET);
        fgets(eucldata,MAX_LINE,fp);	//Read first line of input_file
		fgets(eucldata,MAX_LINE,fp);		//Read second line 
		double *p = malloc(col * sizeof(double));	//Allocate array to store coordinates
		while(fscanf(fp,"%s",item) != EOF)		//String returned is itemK
		{
			for(i=0; i < col; i++)
			{
				fscanf(fp,"%s",token);	//Read coordinates
				p[i] = atof(token);		//Convert coordinates to type double
			}
			if (strcmp(m,"@metriceuclidean") == 0)	//Euclidean metric
			{
				for(i = 0; i < L; i++)	
				{
					euclID = hash_func_Eucl(g[i],p,k,col);
					euclID = abs(euclID);
					pos = mod(euclID , tableSize);
					insert_chain(item,p,&(htable[i].table[pos]),flag,col,euclID);
				}	
			}
			else //Cosine metric
			{
				for(i = 0; i < L; i++)	
				{
					pos = hash_func_Cos(g[i],p,k,col);
					insert_chain(item,p,&(htable[i].table[pos]),flag,col,euclID);
				}	
			}
		}
		free(p);
		/*End of Input phase*/
		/*for(i = 0; i < L; i++)	
		{
			printf("TABLE %d\n",i);
			print_table(htable[3]);	//Print tables for check
		}*/
		files = 1;
		do
		{
			if (files == 1) 	fq = fopen(argv[query],"r");
			else
			{
				printf("Give QueryFile name: ");
				scanf("%s",fileName);
				fq = fopen(fileName,"r");
				printf("Give OutputFile name: ");
				scanf("%s",fileName);
				fe = fopen(fileName,"w+");
			}
			if (fq == NULL)
			{
				perror("Error:");
				return -1;
			}
			/*Search phase*/
			fscanf(fq,"%s%lf[^\n]",radius,&rad);	//Read first line of query_file
			printf("%s	'%lf'\n",radius,rad);
			if (rad > 0)	valid = 1;		//If a positive Radius is given, then compute all neighbors
			else if (rad == 0) valid = 0;	//If Radius equals to zero, then find only THE nearest neighbor
			else 							//If Radius < 0, then report error
			{
				printf("Negative radius found. Try again\n");
				return -1;
			}
			double *q = malloc(col * sizeof(double));
			while(fscanf(fq,"%s",qitem) != EOF)		//String returned is itemK
			{
				//itemid = make_item(qitem);
				for(i=0; i < col; i++)
				{
					fscanf(fq,"%s",token);	//Read coordinates
					q[i] = atof(token);		//Convert coordinates to double
				}
				if (strcmp(m,"@metriceuclidean") == 0)
				{
					/*Brute force starts*/
					bruflag = 1;
					start_t = clock();
					nn tnn;
					tnn = brute_force_table(htable[0],q,flag,bruflag,0,col,L);
					end_t = clock();
					total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
					bruflag = 0;
					/*Brute force ends*/
					/*NNR starts*/
					if (valid)
					{
						for(i = 0; i < L; i++)	
						{
							//printf("TABLE %d\n",i);
							euclID = hash_func_Eucl(g[i],q,k,col);
							euclID = abs(euclID);
							pos = mod(euclID , tableSize);
							//printf("pos=%d\n",pos);
							search_table_NNR(pos,htable[i],q,rad,&nnrlist,flag,euclID,col);		//Search for NNRs
						}
						fprintf(fe,"Query: %s\nR-near neighbors:\n",qitem); 
						destroy_nnrlist(&nnrlist,fe);	
					}
					/*NNR ends*/
					/*NN starts*/
					start_t = clock();
					nn lshnn;			
					lshnn = search_table_NN(g,htable,q,NULL,flag,col,k,L,tableSize);
					end_t = clock();
					total_t1 = (double)(end_t - start_t) / CLOCKS_PER_SEC;
					/*NN ends*/
					fprintf(fe,"Nearest neighbor: %s\n",lshnn.key);
					fprintf(fe,"True neighbor: %s\n",tnn.key);
					fprintf(fe,"distanceLSH: %f\n",lshnn.distance);
					fprintf(fe,"distanceTrue: %f\n",tnn.distance);
					fprintf(fe,"tLSH: %f\n",total_t1);
					fprintf(fe,"tTrue: %f\n\n",total_t);
					free(lshnn.key);
					free(tnn.key);
				}
				else 
				{
					/*Brute force starts*/
					bruflag = 1;
					start_t = clock();
					nn tnn;
					tnn = brute_force_table(htable[0],q,flag,bruflag,0,col,L);
					end_t = clock();
					total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
					bruflag = 0;
					/*Brute force ends*/
					/*NNR starts*/
					if (valid)
					{
						for(i = 0; i < L; i++)	
						{
							pos = hash_func_Cos(g[i],q,k,col);
							search_table_NNR(pos,htable[i],q,rad,&nnrlist,flag,0,col);		//Search for NNRs
						}
						fprintf(fe,"Query: %s\nR-near neighbors:\n",qitem); 
						destroy_nnrlist(&nnrlist,fe);	
					}
					/*NNR ends*/
					/*NN starts*/
					start_t = clock();
					nn lshnn;		
					lshnn = search_table_NN(g,htable,q,NULL,flag,col,k,L,tableSize);
					end_t = clock();
					total_t1 = (double)(end_t - start_t) / CLOCKS_PER_SEC;
					/*NN ends*/
					fprintf(fe,"Nearest neighbor: %s\n",lshnn.key);
					fprintf(fe,"True neighbor: %s\n",tnn.key);
					fprintf(fe,"distanceLSH: %f\n",lshnn.distance);
					fprintf(fe,"distanceTrue: %f\n",tnn.distance);
					fprintf(fe,"tLSH: %f\n",total_t1);
					fprintf(fe,"tTrue: %f\n\n",total_t);
					free(lshnn.key);
					free(tnn.key);
				}
			}
			free(q);
			/*End of search phase*/
			printf("Search completed with success.\n");
			fclose(fq);	
			fclose(fe);
			printf("Would you like to continue the search? Y or N: ");
			scanf("%s",answer);
			files++;
		}while (strcmp(answer,"Y") == 0);
	}
	/*Memory release*/
	for (i = 0; i < L; i++)
	{
		if ((flag == 1) || (flag == 2))
		{
			for (j=0; j < k; j++)	
				free(g[i][j].v);
		}
		free(g[i]);
	}
	free(g);
	for (i=0; i < L; i++)
		destroy_table(&htable[i],flag);	
	free(htable);
	/*Total memory was released*/
	fclose(fp);		//Close Input File
	return 0;
}

char *inputString(FILE *fp,size_t size) 
{
	char *str;
	int ch;
	size_t len = 0;
	str = realloc(NULL,size*sizeof(char));
	if(!str)  	return str;
	while ((ch = fgetc(fp)) != EOF && ch != '\n') 
	{
		str[len++] = ch;
		if(len == size) 
		{
			str = realloc(str,(size += (MAX_LINE / 2))*sizeof(char));
			if(!str) 	return str;
		}
	}
	str[len++] = '\0';
	return realloc(str,len*sizeof(char)); 
}
