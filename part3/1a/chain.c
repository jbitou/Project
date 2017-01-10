#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "chain.h"
#define ITEM_ID 15
#define TRICK 6

void vector_insert_chain(chainp *pointer, user usr, int d, int id) {
	int i;
	chainp temp;
	temp = *pointer;
	/**if list is empty, put the first node**/
	if (temp == NULL) {
		temp = malloc(sizeof(chain));
		temp->usr.userId = usr.userId;
		temp->usr.ratings = malloc(d*sizeof(rating));
		for(i=0; i < d; i++) {
			temp->usr.ratings[i].itemId = usr.ratings[i].itemId;
			temp->usr.ratings[i].rate = usr.ratings[i].rate;
		}			
		temp->euclid = id;
		temp->next = NULL;
		*pointer = temp;
	}
	/**If list isn't empty, put new node at the end**/
	else {
		while (temp->next != NULL)	temp = temp->next;
		temp->next = malloc(sizeof(chain));
		temp->next->usr.userId = usr.userId;
		temp->next->usr.ratings = malloc(d*sizeof(rating));
		for(i=0; i < d; i++) {
			temp->next->usr.ratings[i].itemId = usr.ratings[i].itemId;
			temp->next->usr.ratings[i].rate = usr.ratings[i].rate;
		}			
		temp->next->euclid = id;
		temp->next->next = NULL;
	}
}

nnrlist search_chain_NNR(chainp b, nnrlist list, user usr, double R1, double R2, int euclID, int d, int flag) {
	int ok, exists = 0;
	double diff;
	chainp temp;
	if (flag == 1) {
		/**Count items for which ID(p) = ID(q) ---- Euclidean**/
		temp = b;
		while (temp != NULL) {
			if(euclID == temp->euclid) 	exists++;
			temp = temp->next;
		}
	}
	temp = b;
	/**If exists = 0 ---- Cosine & Euclidean**/
	if(exists <= 1) {
		while (temp != NULL) {
			if (usr.userId == temp->usr.userId) {
				temp = temp->next;
				continue;
			}
			if (flag == 1)	diff = distance_Euclidean(temp->usr.ratings,usr.ratings,d);
			else if (flag == 2)	diff = distance_Cosine(temp->usr.ratings,usr.ratings,d);
			if (diff >= R1 && diff <= R2) {
				ok = insert_nnrlist(&list.list,temp->usr,d,diff);
				if (ok == 1) (list.length)++;
			}
			temp = temp->next;
		}
	}
	else {
		while (temp != NULL) {
			if (usr.userId == temp->usr.userId) {
				temp = temp->next;
				continue;
			}
			if(euclID == temp->euclid) {
				diff = distance_Euclidean(temp->usr.ratings,usr.ratings,d);
				if (diff >= R1 && diff <= R2) {
					ok = insert_nnrlist(&list.list,temp->usr,d,diff);
					if (ok == 1) (list.length)++;
				}
			}
			temp = temp->next;
		}
	}
	return list;
}


void destroy_chain(chainp *l) {
	chainp temp, curr;
	curr = *l;
	if (curr == NULL)		return;		
	while (curr != NULL) {
		temp = curr;
		curr = curr->next;
		free(temp->usr.ratings);		
		free(temp);
	}
	*l = NULL;	
}

void print_chain(chainp l) {
	chainp temp = l;
	while (temp != NULL) {
		printf("user:%d\t",temp->usr.userId);
		temp = temp->next;
	}
	printf("\n");
}

int insert_nnrlist(nnrp *pointer, user usr, int d, double distance) {
	int i;
	nnrp temp, curr, prev;
	temp = curr = *pointer;
	/**If list is empty, put the first node**/
	if (temp == NULL) {
		temp = malloc(sizeof(nnr));
		temp->usr.userId = usr.userId;
		temp->usr.ratings = malloc(d*sizeof(rating));
		for(i=0; i < d; i++) {
			temp->usr.ratings[i].itemId = usr.ratings[i].itemId;
			temp->usr.ratings[i].rate = usr.ratings[i].rate;
		}	
		temp->distance = distance;
		temp->next = NULL;
		*pointer = temp;
	}
	/**If list isn't empty, put new node at the end**/
	else {
		/**Change start**/
		if (curr->distance > distance) {
			curr = malloc(sizeof(nnr));
			curr->usr.userId = usr.userId;
			curr->usr.ratings = malloc(d*sizeof(rating));
			for(i=0; i < d; i++) {
				curr->usr.ratings[i].itemId = usr.ratings[i].itemId;
				curr->usr.ratings[i].rate = usr.ratings[i].rate;
			}	
			curr->distance = distance;
			curr->next = *pointer;
			*pointer = curr;
			return 1;
		}
		if (curr->usr.userId == usr.userId)  return 0;
		prev = temp;
		temp = temp->next;
		while(temp != NULL && temp->distance <= distance) {
			if (temp->usr.userId == usr.userId)	return 0; 
			prev = temp;
			temp = temp->next;
		}
		curr = malloc(sizeof(nnr));
		curr->usr.userId = usr.userId;
		curr->usr.ratings = malloc(d*sizeof(rating));
		for(i=0; i < d; i++) {
			curr->usr.ratings[i].itemId = usr.ratings[i].itemId;
			curr->usr.ratings[i].rate = usr.ratings[i].rate;
		}	
		curr->distance = distance;
		prev->next = curr;
		curr->next = temp;
	}
	return 1;
}

int insert_ratinglist(ratinglist *pointer, int itemId, double rate) {
	int i;
	ratinglist temp, curr, prev;
	temp = curr = *pointer;
	/**If list is empty, put the first node**/
	if (temp == NULL) {
		temp = malloc(sizeof(recrating));	
		temp->itemId = itemId;
		temp->rate = rate;
		temp->next = NULL;
		*pointer = temp;
	}
	/**If list isn't empty, put new node at the end**/
	else {
		/**Change start**/
		if (curr->rate < rate) {
			curr = malloc(sizeof(recrating));
			curr->itemId = itemId;
			curr->rate = rate;
			curr->next = *pointer;
			*pointer = curr;
			return 1;
		}
		prev = temp;
		temp = temp->next;
		while(temp != NULL && temp->rate >= rate) {
			prev = temp;
			temp = temp->next;
		}
		curr = malloc(sizeof(recrating));
		curr->itemId = itemId;
		curr->rate = rate;
		prev->next = curr;
		curr->next = temp;
	}
	return 1;
}

nnrlist combine_nnrlist(nnrlist l1, nnrlist *l2, int d, int stop) {
	int i = 0;
	nnrp temp;
	while (l2->list != NULL) {
		if (i < stop) {
			insert_nnrlist(&l1.list,l2->list->usr,d,l2->list->distance);
			(l1.length)++;
		}
		temp = l2->list;
		l2->list = l2->list->next;
		free(temp->usr.ratings);
		free(temp);
		i++;
	}
	return l1;
}

int nnrlist_length(nnrp l) {
	int length = 0;
	nnrp temp = l;
	while (temp != NULL) {
		length++;
		temp = temp->next;
	}
	return length;
}

void destroy_nnrlist(nnrp *l)  {
	nnrp temp, curr;
	curr = *l;
	if (curr == NULL)		return;		
	while (curr != NULL) {
		temp = curr;
		curr = curr->next;
		free(temp->usr.ratings);
		free(temp);
	}
	*l = NULL;
}

void destroy_ratinglist(ratinglist *l) {
	ratinglist temp, curr;
	curr = *l;
	if (curr == NULL)		return;		
	while (curr != NULL) {
		temp = curr;
		curr = curr->next;
		free(temp);
	}
	*l = NULL;
}

void printndestroy_ratinglist(ratinglist *l, int id, FILE *fp) {
	int i = 0;
	ratinglist temp, curr;
	curr = *l;
	if (curr == NULL)		return;	
	fprintf(fp,"u%d\t",id);
	while (curr != NULL) {
		if (i < 5)	fprintf(fp,"%d\t",curr->itemId);
		temp = curr;
		curr = curr->next;
		free(temp);
		i++;
	}
	fprintf(fp,"\n");
	*l = NULL;
}

void print_nnrlist(nnrp l) {
	nnrp temp = l;
	while (temp != NULL) {
		printf("user: %d distance %.20lf\n",temp->usr.userId,temp->distance);
		temp = temp->next;
	}
	printf("\n");
}

void print_ratinglist(ratinglist l) {
	ratinglist temp = l;
	while (temp != NULL) {
		printf("item: %d rate %.20lf\n",temp->itemId,temp->rate);
		temp = temp->next;
	}
	printf("\n");
}
