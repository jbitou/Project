#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "NNsearch.h"

void binary_repeated_nnsearch(user *users, int numofusers, int numofitems, int P) {
	int i, j;
	double radius;
	radius = find_first_radius(users, numofusers);
	/**Search P neighboors for each user**/
	/*for (i=0; i < numofusers; i++) {
	}*/
}


double find_first_radius(user *users, int numofusers) {
	int user1, user2;
	user1 =  1 + (rand() / (RAND_MAX + 1.0)) * numofusers;
	user2 =  1 + (rand() / (RAND_MAX + 1.0)) * numofusers;
	while (user2 == user1)	user2 =  1 + (rand() / (RAND_MAX + 1.0)) * numofusers;
	/**noooope**/
	return distance_Euclidean(users[user1].ratings,users[user2].ratings,numofitems);
}


double distance_Euclidean(user *v1, user *v2, int d) {
	int i;
	double distance = 0.0, diff = 0.0 ;
	for (i=0; i < d; i++) {
		diff = v1[i].rate - v2[i].rate;
		distance += diff*diff;
	}
	return sqrt(distance);
}
