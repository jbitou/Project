#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "NNsearch.h"
#define L 6
#define k 5
#define epsilon 0.0000000001

void nnlsh_recommendation(user *users, int numofusers, int numofitems, int P) {
	int i, j;
	double radius;
	hash_table *htable;
	ghashp *g;
	nnrlist neighbours;
	g = malloc(L * sizeof(ghashp));		
	for(i=0; i < L; i++)	g[i] = malloc(k * sizeof(ghash));
	init_hash_Eucl(g,L,k,numofitems);
	printf("start insert\n");
	htable = vector_hash_insert(g,users,L,k,numofusers,numofitems);
	printf("done insert\n");
	/*for(i = 0; i < L; i++)	{
		printf("table %d: \n",i);
		for(j = 0; j < numofusers/16+1; j++) {
			printf("bucket %d: \n",j);
			print_chain(htable[i].table[j]);
		}
	}*/
	/**Search P neighboors for each user**/
	radius = find_first_radius(users,numofusers,numofitems);
	printf("First radius %lf\n",radius);
	for (i=0; i < numofusers; i++) {
		neighbours.list = NULL;
		neighbours.length = 0;
		neighbours = binary_repeated_nnsearch(htable,g,neighbours,users,i,P,numofusers,numofitems,radius);
		printf("Neighbours of user %d\n",users[i].userId);
		//print_nnrlist(neighbours.list);
		
	}
	for (i=0; i < L; i++) {
		destroy_table(&htable[i]);
		for (j=0; j < k; j++)	free(g[i][j].v);
		free(g[i]);
	}
	free(g);
	free(htable);
}

nnrlist binary_repeated_nnsearch(hash_table *htable, ghashp *g, nnrlist neighbours, user *users, int i, int P, int numofusers, int numofitems, double radius) {
	int previous;
	double previousr, middler, tempr;
	nnrlist extra, nnlist;
	while (neighbours.length != P) {
		neighbours = search_table_NNR(htable,g,neighbours,users[i],0,radius,numofusers,numofitems,L,k);
		//printf("radius=%.30lf previousr=%.30lf start length = %d\n",radius,previousr,neighbours.length);
		if (neighbours.length < P) {
			if (previous > P) {
				tempr = previousr;
				do {
					destroy_nnrlist(&neighbours.list);
					neighbours.length = 0;
					middler = (radius+previousr) / 2;
					neighbours = search_table_NNR(htable,g,neighbours,users[i],0,middler,numofusers,numofitems,L,k);
					//printf("radius=%lf previousr=%lf middler=%lf 1a.length=%d\n",radius,previousr,middler,neighbours.length);
					previousr = middler;
					if (middler <= epsilon) {
						nnlist.list = NULL;
						nnlist.length = 0;
						combine_nnrlist(&nnlist,&neighbours,numofitems,P);
						return nnlist;
					}
				}while (neighbours.length > P);
				if (neighbours.length < P) {
					extra.list = NULL;
					extra.length = 0;
					extra = search_table_NNR(htable,g,neighbours,users[i],middler,tempr,numofusers,numofitems,L,k);
					combine_nnrlist(&neighbours,&extra,numofitems,P-neighbours.length);
					//printf("1.length=%d\n",neighbours.length);
					destroy_nnrlist(&extra.list);
				}
				break;
			} 
			else if (previous < P) {
				previousr = radius;
				radius *= 2;
			}
		}
		else if (neighbours.length > P) {
			if (previous < P) {
				tempr = radius;
				do {
					destroy_nnrlist(&neighbours.list);
					neighbours.length = 0;
					middler = (radius+previousr) / 2;
					neighbours = search_table_NNR(htable,g,neighbours,users[i],0,middler,numofusers,numofitems,L,k);
					//printf("radius=%lf previousr=%lf middler=%lf 2a.length=%d\n",radius,previousr,middler,neighbours.length);
					radius = middler;
					if (middler <= epsilon) {
						nnlist.list = NULL;
						nnlist.length = 0;
						combine_nnrlist(&nnlist,&neighbours,numofitems,P);
						return nnlist;
					}
				}while (neighbours.length > P);
				if (neighbours.length < P) {
					extra.list = NULL;
					extra.length = 0;
					extra = search_table_NNR(htable,g,neighbours,users[i],middler,tempr,numofusers,numofitems,L,k);
					combine_nnrlist(&neighbours,&extra,numofitems,P-neighbours.length);
					//printf("radius=%lf previousr=%lf middler=%lf 2b.length=%d extralen=%d\n",radius,previousr,middler,neighbours.length,extra.length);
					destroy_nnrlist(&extra.list);
				}
				break;
			}
			else if (previous > P) {
				previousr = radius;
				radius /= 2;
			}
		}
		previous = neighbours.length;
		if (radius <= epsilon) {
			nnlist.list = NULL;
			nnlist.length = 0;
			combine_nnrlist(&nnlist,&neighbours,numofitems,P);
			return nnlist;
		}
		destroy_nnrlist(&neighbours.list);
		neighbours.length = 0;
	}	
	return neighbours;
}


double find_first_radius(user *users, int numofusers, int numofitems) {
	int j, user1, user2;
	user1 =  (rand() / (RAND_MAX + 1.0)) * numofusers;
	user2 =  (rand() / (RAND_MAX + 1.0)) * numofusers;
	while (user2 == user1)	user2 =  (rand() / (RAND_MAX + 1.0)) * numofusers;
	return distance_Euclidean(users[user1].ratings,users[user2].ratings,numofitems) / 2;
}

 
