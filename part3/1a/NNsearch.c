#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "NNsearch.h"
#define FIRST_RADIUS 0.5
#define epsilon 0.00001

double nnlsh_recommendation(user *users, user *vset, int vsize, int numofusers, int numofitems, int P, FILE *fp, int flag) {
	int i, j, J = 0;
	double radius, z, total_t, sum = 0.0;
	hash_table *htable;
	ghashp *g;
	nnrlist neighbours;
	ratinglist recommendations, temp;
	clock_t start_t, end_t;
	if (fp != NULL)	{
		if (flag == 1)	fprintf(fp,"Euclidean LSH\n"); 
		else if (flag == 2)	fprintf(fp,"Cosine LSH\n"); 
	}
	/**Initialize g and hash functions**/
	g = malloc(L * sizeof(ghashp));		
	for(i=0; i < L; i++)	g[i] = malloc(k * sizeof(ghash));
	if (flag == 1) init_hash_Eucl(g,numofitems);
	else if (flag == 2)	init_hash_Cos(g,numofitems);
	/**Insert users into hash tables**/
	htable = vector_hash_insert(g,users,numofusers,numofitems,flag);
	if (fp != NULL) start_t = clock();
	/**Compute start radius**/
	radius = find_first_radius(vset,vsize,numofitems,flag);
	printf("First radius %lf\n",radius);
	/**Search P neighboors for each user**/
	for (i=0; i < vsize; i++) {
		neighbours = binary_repeated_nnsearch(htable,g,vset,i,P,numofusers,numofitems,radius,flag);
		if (vsize == 20) {
			recommendations = evaluate_rest_items(vset[i],neighbours,numofitems);
			/*printf("items for user %d\n",vset[i].userId);
			print_ratinglist(recommendations);*/
			printndestroy_ratinglist(&recommendations,vset[i].userId,fp);
		}
		else if (vsize == f) {
			recommendations = evaluate_items(vset[i],neighbours,numofitems);
			printf("items for user %d\n",vset[i].userId);
			print_ratinglist(recommendations);
			/**Compare predicted with real evaluation**/
			for (j=0; j < numofitems; j++) {
				if (vset[i].ratings[j].rate != 0) {
					J++;
					temp = recommendations;
					while (temp != NULL) {
						if (temp->itemId == vset[i].ratings[j].itemId) {
							printf("For user %d item %d real rating=%lf predicted=%lf\n",vset[i].userId,temp->itemId,vset[i].ratings[j].rate,temp->rate);
							sum += fabs(temp->rate  - vset[i].ratings[j].rate);
							break;
						}
						temp = temp->next;
					}
				}
			}
			destroy_ratinglist(&recommendations);
		}
		printf("Neighbours of user %d are %d with start %lf\n",vset[i].userId,neighbours.length,radius);
		print_nnrlist(neighbours.list);	
		destroy_nnrlist(&neighbours.list);
	}
	if (fp != NULL) {
		end_t = clock();
		total_t = ((double)(end_t - start_t) / CLOCKS_PER_SEC)*1000;
		fprintf(fp,"Execution Time: %.6lf\n",total_t);
	}
	for (i=0; i < L; i++) {
		destroy_table(&htable[i]);
		for (j=0; j < k; j++)	free(g[i][j].v);
		free(g[i]);
	}
	free(g);
	free(htable);
	printf("sum = %lf, J=%d\n",sum,J);
	return sum / J;
}

nnrlist binary_repeated_nnsearch(hash_table *htable, ghashp *g, user *users, int i, int P, int numofusers, int numofitems, double radius, int flag) {
	int previous = 0, same = 0;
	double previousr = 0.0, middler = 0.0, tempr;
	nnrlist neighbours, extra, nnlist;
	neighbours.list = NULL;
	neighbours.length = 0;
	while (neighbours.length != P) {
		neighbours = search_table_NNR(htable,g,users[i],0,radius,numofusers,numofitems,flag);
		printf("radius start %lf	length %d\n",radius,neighbours.length);
		/**If we cannot reach P(not enough users in buckets)**/
		if (neighbours.length == previous)	same++;
		else 	same = 0;
		if (same > 10) return neighbours;
		if (neighbours.length < P) {
			if (previous > P) {
				tempr = previousr;
				do {
					destroy_nnrlist(&neighbours.list);
					neighbours.length = 0;
					middler = (radius+previousr) / 2;
					neighbours = search_table_NNR(htable,g,users[i],0,middler,numofusers,numofitems,flag);
					printf("radius=%lf previousr=%lf middler=%lf 1a.length=%d\n",radius,previousr,middler,neighbours.length);
					/**If radius is approximately zero or previousr with middler have a negligible distance**/
					if (middler <= epsilon || fabs(previousr-middler) <= epsilon) {
						nnlist.list = NULL;
						nnlist.length = 0;
						nnlist = combine_nnrlist(nnlist,&neighbours,numofitems,P);
						return nnlist;
					}
					previousr = middler;
				}while (neighbours.length > P);
				if (neighbours.length < P) {
					extra.list = NULL;
					extra.length = 0;
					extra = search_table_NNR(htable,g,users[i],middler,tempr,numofusers,numofitems,flag);
					neighbours = combine_nnrlist(neighbours,&extra,numofitems,P-neighbours.length);
					printf("radius=%lf previousr=%lf middler=%lf 1b.length=%d extralen=%d\n",radius,previousr,middler,neighbours.length,extra.length);
				}
				break;
			} 
			else if (previous < P || previous == 0) {
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
					neighbours = search_table_NNR(htable,g,users[i],0,middler,numofusers,numofitems,flag);
					printf("radius=%lf previousr=%lf middler=%lf 2a.length=%d\n",radius,previousr,middler,neighbours.length);
					/**If radius is approximately zero or radius with middler have a negligible distance**/
					if (middler <= epsilon || fabs(radius-middler) <= epsilon) {
						nnlist.list = NULL;
						nnlist.length = 0;
						nnlist = combine_nnrlist(nnlist,&neighbours,numofitems,P);
						return nnlist;
					}
					radius = middler;
				}while (neighbours.length > P);
				if (neighbours.length < P) {
					extra.list = NULL;
					extra.length = 0;
					extra = search_table_NNR(htable,g,users[i],middler,tempr,numofusers,numofitems,flag);
					neighbours = combine_nnrlist(neighbours,&extra,numofitems,P-neighbours.length);
					printf("radius=%lf previousr=%lf middler=%lf 2b.length=%d extralen=%d\n",radius,previousr,middler,neighbours.length,extra.length);
				}
				break;
			}
			else if (previous > P) {
				previousr = radius;
				radius /= 2;
			}
		}
		else return neighbours;
		previous = neighbours.length;
		/**If radius is approximately zero**/
		if (radius <= epsilon) {
			nnlist.list = NULL;
			nnlist.length = 0;
			nnlist = combine_nnrlist(nnlist,&neighbours,numofitems,P);
			return nnlist;
		}
		destroy_nnrlist(&neighbours.list);
		neighbours.length = 0;
	}	
	return neighbours;
} 

ratinglist evaluate_rest_items(user usr, nnrlist neighbours, int numofitems) {
	int i;
	double z, sum, rate;
	nnrp temp;
	ratinglist recommendations = NULL;
	z = find_normalizing_factor(neighbours);
	for (i=0; i < numofitems; i++) {
		/**Exclude items that have been rated**/
		if (usr.ratings[i].rate != 0) continue;
		sum = 0.0;
		temp = neighbours.list;
		while (temp != NULL) {
			//printf("sum=%lf and rate=%lf item %d neighbour %d\n",sum,temp->usr.ratings[i].rate,temp->usr.ratings[i].itemId,temp->usr.userId);
			sum += (1 / (1 + temp->distance)) * temp->usr.ratings[i].rate;
			temp = temp->next;
		}
		rate = z * sum;
		if (rate != 0)	insert_ratinglist(&recommendations,usr.ratings[i].itemId,rate);
	}
	return recommendations;
}

ratinglist evaluate_items(user usr, nnrlist neighbours, int numofitems) {
	int i;
	double z, sum, rate;
	nnrp temp;
	ratinglist recommendations = NULL;
	z = find_normalizing_factor(neighbours);
	for (i=0; i < numofitems; i++) {
		/**Exclude items that haven't been rated**/
		if (usr.ratings[i].rate == 0) continue;
		sum = 0.0;
		temp = neighbours.list;
		while (temp != NULL) {
			//printf("sum=%lf and rate=%lf\n",sum,temp->usr.ratings[i].rate);
			sum += (1 / (1 + temp->distance)) * temp->usr.ratings[i].rate;
			temp = temp->next;
		}
		rate = z * sum;
		insert_ratinglist(&recommendations,usr.ratings[i].itemId,rate);
	}
	return recommendations;
}

double find_normalizing_factor(nnrlist neighbours) {
	double sum = 0.0, sim;
	nnrp temp = neighbours.list;
	while (temp != NULL) {
		sim = 1 / (1 + temp->distance);
		sum += fabs(sim);
		temp = temp->next;
	}
	return 1/sum;
}

double find_first_radius(user *users, int numofusers, int numofitems, int flag) {
	int j, user1, user2;
	double dis;
	do {
		user1 =  (rand() / (RAND_MAX + 1.0)) * numofusers;
		user2 =  (rand() / (RAND_MAX + 1.0)) * numofusers;
		while (user2 == user1)	user2 =  (rand() / (RAND_MAX + 1.0)) * numofusers;
		if (flag == 1)	dis = distance_Euclidean(users[user1].ratings,users[user2].ratings,numofitems) / 2;
		else if (flag == 2)	dis = distance_Cosine(users[user1].ratings,users[user2].ratings,numofitems) / 2;
		printf("dis=%lf\n",dis);
	}while (dis < FIRST_RADIUS);
	return dis;
}


 
