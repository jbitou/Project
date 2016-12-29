#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "points_list.h"

int insert_points(pointp *list, double distance1, double distance2, centroid second, int position) {
	pointp temp;
	temp = *list;
	/**If list is empty, put the first node**/
	if (temp == NULL) {
		temp = malloc(sizeof(point));
		temp->position = position;
		temp->mindistance = distance1;
		temp->secdistance = distance2;
		temp->second = second;
		temp->next = NULL;
		*list = temp;
	}
	/**If list isn't empty, put new node at the end**/
	else {
		while (temp->next != NULL)	temp = temp->next;
		temp->next = malloc(sizeof(point));
		temp->next->position = position;
		temp->next->mindistance = distance1;
		temp->next->secdistance = distance2;
		temp->next->second = second;
		temp->next->next = NULL;
	}
	return 0;
}

pointp clone_points(pointp list) {
    if (list == NULL) return NULL;
    pointp result = malloc(sizeof(point));
	result->mindistance = list->mindistance;
	result->secdistance = list->secdistance;
	result->second = list->second;
	result->position = list->position;
    result->next = clone_points(list->next);
    return result;
}

int points_length(pointp l) {
	int length = 0;
	pointp temp = l;
	while (temp != NULL) {
		length++;
		temp = temp->next;
	}
	return length;
}

void destroy_points(pointp *l) {
	pointp temp, curr;
	curr = *l;
	if (curr == NULL)		return;		
	while (curr != NULL) {
		temp = curr;
		curr = curr->next;		
		free(temp);
	}
	*l = NULL;	
}

void printndestroy_points(pointp *l, FILE *fp) {
	pointp temp, curr;
	curr = *l;
	if (curr == NULL)		return;		
	while (curr != NULL) {
		fprintf(fp,"%d\t",curr->position+1);
		temp = curr;
		curr = curr->next;		
		free(temp);
	}
	fprintf(fp,"\n");
	*l = NULL;	
}

void print_points(pointp l) {
	pointp temp = l;
	while (temp != NULL) {
		printf("\nconform: %d, mindistance: %lf\n",temp->position,temp->mindistance);
		//printf(" second: cluster %d, secdistance: %lf\n",temp->second,temp->secdistance);
		temp = temp->next;
	}
}
