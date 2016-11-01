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

int main (int argc, char **argv)
{
	FILE *fp, *fq, *fe;
	int i, j, foundK, foundL, f1, f2, f3, k, L, key, bs, pos, input, query, output, flag, valid, tableSize, bruflag, files;
	char item[ITEM_ID], qitem[ITEM_ID], ms[14], data[65], qdata[65], space[10], radius[8], answer[2], fileName[50];
	double rad, true_distance, min_distance;
	hash_table *htable;
	nnrp nnrlist = NULL, nnlist = NULL, tnnlist = NULL;
	clock_t start_t, end_t;
	double total_t, total_t1;
	
	if (argc > 11)
	{
		printf("Too many arguments. Try again.\n");
		return -1;
	}
	else if (argc < 6)
	{
		printf("Too few arguments. Try again.\n");
		return -1;
	}
	/*Every parameter has to be given after a recogniser*/
	if ((argc % 2) == 0)
	{
		printf("Wrong arguments. Try again.\n");	
		return -1;
	}	
	/*Initialize random number generator*/
	srand(time(NULL));	
	foundK = foundL = 0;
	f1 = f2 = f3 = 0;
	for (i=1; i < (argc-1); i+=2)
	{
		if (strcmp(argv[i],"-k") == 0)
		{
			foundK = 1;
			k = atoi(argv[i+1]); 
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
	ghashp *g = malloc(L * sizeof(ghashp));		
	for(i = 0; i < L; i++) 
		g[i] = malloc(k * sizeof(ghash));
	/*Read first line of input_file*/	
	fscanf(fp,"%s%s[^\n]",ms,space);	
	if (strcmp(space,"hamming") == 0)
	{
		flag = 0;
		tableSize = 1 << k;
		for (i=0; i<L; i++)
			init_table(k,&htable[i],tableSize);
		/*Read  datasets' first line to compute 'N'(number of bits)*/
		fscanf(fp,"%s %s[^\n]",item,data);	
		init_hash_Ham(g,L,k,data);
		/*Go back to the start*/
		fseek(fp,0,SEEK_SET);
		fscanf(fp,"%s%s[^\n]",ms,space);
		/*Ιnput phase*/
		while (fscanf(fp,"%s %s[^\n]",item,data) != EOF)		
		{   
			for(i = 0; i < L; i++)	
			{
				pos = hash_func_Ham(g[i],data,k);
				insert_chain(item,data,&(htable[i].table[pos]),flag,0,0);
			}
		}
		/*End of Ιnput phase*/
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
			/*Read first line of query_file*/
			fscanf(fq,"%s%lf[^\n]",radius,&rad);
			/*If a positive Radius is given, then compute all neighbors else only the nearest one*/	
			if (rad > 0)	valid = 1;		
			else if (rad == 0) valid = 0;	
			else 							
			{
				printf("Negative radius found. Try again\n");
				return -1;
			}
			/*Search phase*/
			while (fscanf(fq,"%s %s[^\n]",qitem,qdata) != EOF)		
			{
				/*Brute force starts*/
				bruflag = 1;
				start_t = clock();
				brute_force_table(htable[0],qdata,flag,bruflag,0,0,L,&tnnlist,&true_distance);
				end_t = clock();
				total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
				bruflag = 0;
				/*Brute force ends*/
				/*NNR starts*/
				if (valid)
				{
					for (i=0; i < L; i++)
					{
						pos = hash_func_Ham(g[i],qdata,k);	
						search_table_NNR(pos,htable[i],qdata,rad,&nnrlist,flag,0,0);
					}
					fprintf(fe,"Query: %s\nR-near neighbors:\n",qitem); 
					print_nnrlist(&nnrlist,fe);
				}
				/*NNR ends*/
				/*NN starts*/
				start_t = clock();
				search_table_NN(g,htable,qdata,NULL,flag,0,k,L,tableSize,&nnlist,&min_distance);
				end_t = clock();
				total_t1 = (double)(end_t - start_t) / CLOCKS_PER_SEC;
				/*NN ends*/
				fprintf(fe,"Nearest neighbors: \n");
				print_nnrlist(&nnlist,fe);
				fprintf(fe,"True neighbors: \n");
				print_nnrlist(&tnnlist,fe);
				fprintf(fe,"distanceLSH: %.0f\n",min_distance);
				fprintf(fe,"distanceTrue: %.0f\n",true_distance);
				fprintf(fe,"tLSH: %f\n",total_t1);
				fprintf(fe,"tTrue: %f\n\n",total_t);
			}	/*End of search phase*/
			printf("Search completed with success.\n");
			fclose(fq);	
			fclose(fe);
			destroy_nnrlist(&nnrlist);
			destroy_nnrlist(&nnlist);
			destroy_nnrlist(&tnnlist);
			printf("Would you like to continue the search? Y or N: ");
			scanf("%s",answer);
			files++;
		}while (strcmp(answer,"Y") == 0);
	}
	else if (strcmp(space,"matrix") == 0)
	{
		char itemsline[7],*allitems;
		int numofitems, token, itemid;
		while (k > 10)
		{
			printf("Too big k. Try again: ");
			scanf("%d", &k);
		} 
		while (L > 30)
		{
			printf("Too big L. Try again: ");
			scanf("%d", &L);
		} 
		tableSize = 1 << k;
		flag = 3;	
		fscanf(fp,"%s",itemsline);
		allitems = inputString(fp,MAX_LINE);
		/*Read line with items to get the size of te matrix (numofitems x numofitems)*/
		i = 0;
		numofitems = 1;
		while( allitems[i] != '\0')		
		{
			if( allitems[i] == ',')
				numofitems++;
			i++;	
		}
		for (i=0; i < L; i++)	
			init_table(k,&htable[i],tableSize);
		/*Read matrix*/
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
			if (token == 0)  flag1 = 1;	
			/*From now on, read next distance and store it*/
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
		init_hash_matrix(g,p,L,k,numofitems);
		/*Start input*/
		char itemID[ITEM_ID];
		for(i = 0; i < L; i++)	
		{ 
			for(j = 0; j < numofitems; j++)	
			{
				sprintf(itemID,"item%d",j+1);
				pos = hash_func_Matrix(g[i],j,p,k,numofitems);
				insert_chain(itemID,NULL,&(htable[i].table[pos]),flag,0,0);
			}
		}
		/*End of input phase*/
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
			/*Read file, get radius*/
			fscanf(fq,"%s%lf[^\n]",radius,&rad);
			if (rad > 0)	valid = 1;		
			else if (rad == 0) valid = 0;
			else 							
			{
				printf("Negative radius found. Try again\n");
				return -1;
			}
			/*Search starts*/
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
					/*NNR starts*/
					for (i=0; i < L;i++)
					{
						pos = hash_func_MSearch(g[i],qdata,p,k,numofitems);
						search_table_NNR(pos,htable[i],qdata,rad,&nnrlist,flag,0,0);	
					}
					fprintf(fe,"Query: %s\nR-near neighbors:\n",item); 
					print_nnrlist(&nnrlist,fe);
					/*NNR ends*/
				}
				/*Brute force starts*/
				bruflag = 1;
				start_t = clock();
				brute_force_table(htable[0],qdata,flag,bruflag,0,0,L,&tnnlist,&true_distance);
				end_t = clock();
				total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
				bruflag = 0;
				/*Brute force ends*/
				/*NN starts*/
				start_t = clock();
				search_table_NN(g,htable,qdata,p,flag,numofitems,k,L,tableSize,&nnlist,&min_distance);	
				end_t = clock();
				total_t1 = (double)(end_t - start_t) / CLOCKS_PER_SEC;
				/*NN ends*/
				fprintf(fe,"Nearest neighbors: \n");
				print_nnrlist(&nnlist,fe);
				fprintf(fe,"True neighbors: \n");
				print_nnrlist(&tnnlist,fe);
				fprintf(fe,"distanceLSH: %.0f\n",min_distance);
				fprintf(fe,"distanceTrue: %.0f\n",true_distance);
				fprintf(fe,"tLSH: %f\n",total_t1);
				fprintf(fe,"tTrue: %f\n\n",total_t);
			}
			free(qdata);
			printf("Search completed with success.\n");
			fclose(fq);	
			fclose(fe);
			destroy_nnrlist(&nnrlist);
			destroy_nnrlist(&nnlist);
			destroy_nnrlist(&tnnlist);
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
		while (k > 10)
		{
			printf("Too big k. Try again: ");
			scanf("%d", &k);
		} 
		while (L > 30)
		{
			printf("Too big L. Try again: ");
			scanf("%d", &L);
		} 
		/*Read \n of first line, go to second line*/
		fscanf(fp,"\n");
		/*Read only second line*/	
		i = 0;
		while (i < 19)
		{
			c = fgetc(fp);
			if (c == '\n')	break;	
			metric[i++] = c;
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
		/*Read first line to count the dimensions*/
		fscanf(fp,"%s",itemID);
		eucldata = inputString(fp,MAX_LINE);
		/*Count dimensions*/
		ptr = eucldata;
		strtok(ptr,delim);
		ptr = NULL;
		col = 1;		
		while ((point = strtok(ptr,delim)) != NULL)   
		{
			if (lines == -2)	col++;
			ptr = NULL; 
		}
		/*Dimensions found*/
		/*Return file to start*/
		fseek(fp,0,SEEK_SET);	
		/*Count lines*/
		int ch;
		while ((ch = fgetc(fp)) != EOF) 
		{
			if(ch == '\n') lines++;
		}
		/*Lines found*/
		if (strcmp(m,"@metriceuclidean") == 0)
		{
			/*Heuristic choice of tableSize*/
			tableSize = lines/16 + 1;
			flag = 1;
		}
		else 
		{
			tableSize = 1 << k;
			flag = 2;
		}		
		for (i=0; i < L; i++)	
			init_table(k,&htable[i],tableSize);
		if (strcmp(m,"@metriceuclidean") == 0)	
		{
			init_hash_Eucl(g,L,k,col);	
		}
		else 	
		{
			init_hash_Cos(g,L,k,col);	
		}
		/*Input phase*/
		/*Go back to start and read first and second line of input_file*/
        fseek(fp,0,SEEK_SET);
        fgets(eucldata,MAX_LINE,fp);	
		fgets(eucldata,MAX_LINE,fp);	
		/*Allocate array to store coordinates*/	
		double *p = malloc(col * sizeof(double));
		/*Fscanf will return itemK in item*/
		while(fscanf(fp,"%s",item) != EOF)		
		{
			/*Read coordinates*/
			for(i=0; i < col; i++)
			{
				fscanf(fp,"%s",token);	
				p[i] = atof(token);		
			}
			if (strcmp(m,"@metriceuclidean") == 0)	
			{
				for(i = 0; i < L; i++)	
				{
					euclID = hash_func_Eucl(g[i],p,k,col);
					euclID = abs(euclID);
					pos = mod(euclID , tableSize);
					insert_chain(item,p,&(htable[i].table[pos]),flag,col,euclID);
				}	
			}
			else 
			{
				for(i = 0; i < L; i++)	
				{
					pos = hash_func_Cos(g[i],p,k,col);
					insert_chain(item,p,&(htable[i].table[pos]),flag,col,euclID);
				}	
			}
		}
		free(p);
		files = 1;
		/*Search phase*/
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
			/*Read file, get radius*/
			fscanf(fq,"%s%lf[^\n]",radius,&rad);	
			if (rad > 0)	valid = 1;		
			else if (rad == 0) valid = 0;	
			else 							
			{
				printf("Negative radius found. Try again\n");
				return -1;
			}
			double *q = malloc(col * sizeof(double));
			/*Fscanf will return itemK in item*/
			while(fscanf(fq,"%s",qitem) != EOF)		
			{
				/*Read coordinates*/
				for(i=0; i < col; i++)
				{
					fscanf(fq,"%s",token);	
					q[i] = atof(token);		
				}
				if (strcmp(m,"@metriceuclidean") == 0)
				{
					/*Brute force starts*/
					bruflag = 1;
					start_t = clock();
					brute_force_table(htable[0],q,flag,bruflag,0,col,L,&tnnlist,&true_distance);
					end_t = clock();
					total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
					bruflag = 0;
					/*Brute force ends*/
					/*NNR starts*/
					if (valid)
					{
						for(i = 0; i < L; i++)	
						{
							euclID = hash_func_Eucl(g[i],q,k,col);
							euclID = abs(euclID);
							pos = mod(euclID , tableSize);
							search_table_NNR(pos,htable[i],q,rad,&nnrlist,flag,euclID,col);		
						}
						fprintf(fe,"Query: %s\nR-near neighbors:\n",qitem); 
						print_nnrlist(&nnrlist,fe);	
					}
					/*NNR ends*/
					/*NN starts*/
					start_t = clock();
					search_table_NN(g,htable,q,NULL,flag,col,k,L,tableSize,&nnlist,&min_distance);
					end_t = clock();
					total_t1 = (double)(end_t - start_t) / CLOCKS_PER_SEC;
					/*NN ends*/
					fprintf(fe,"Nearest neighbors: \n");
					print_nnrlist(&nnlist,fe);
					fprintf(fe,"True neighbors: \n");
					print_nnrlist(&tnnlist,fe);
					fprintf(fe,"distanceLSH: %f\n",min_distance);
					fprintf(fe,"distanceTrue: %f\n",true_distance);
					fprintf(fe,"tLSH: %f\n",total_t1);
					fprintf(fe,"tTrue: %f\n\n",total_t);
				}
				else 
				{
					/*Brute force starts*/
					bruflag = 1;
					start_t = clock();
					brute_force_table(htable[0],q,flag,bruflag,0,col,L,&tnnlist,&true_distance);
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
						print_nnrlist(&nnrlist,fe);	
					}
					/*NNR ends*/
					/*NN starts*/
					start_t = clock();
					search_table_NN(g,htable,q,NULL,flag,col,k,L,tableSize,&nnlist,&min_distance);
					end_t = clock();
					total_t1 = (double)(end_t - start_t) / CLOCKS_PER_SEC;
					/*NN ends*/
					fprintf(fe,"Nearest neighbors: \n");
					print_nnrlist(&nnlist,fe);
					fprintf(fe,"True neighbors: \n");
					print_nnrlist(&tnnlist,fe);
					fprintf(fe,"distanceLSH: %f\n",min_distance);
					fprintf(fe,"distanceTrue: %f\n",true_distance);
					fprintf(fe,"tLSH: %f\n",total_t1);
					fprintf(fe,"tTrue: %f\n\n",total_t);
				}
			}
			free(q);
			/*End of search phase*/
			printf("Search completed with success.\n");
			fclose(fq);	
			fclose(fe);
			destroy_nnrlist(&nnrlist);
			destroy_nnrlist(&nnlist);
			destroy_nnrlist(&tnnlist);
			printf("Would you like to continue the search? Y or N: ");
			scanf("%s",answer);
			files++;
		}while (strcmp(answer,"Y") == 0);
		free(eucldata);
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
	/*Close Input File*/
	fclose(fp);		
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
