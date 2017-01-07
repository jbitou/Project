#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "distances.h"

double distance_Euclidean(ratingp v1, ratingp v2, int d) {
	int i;
	double distance = 0.0, diff = 0.0 ;
	for (i=0; i < d; i++) {
		diff = v1[i].rate - v2[i].rate;
		distance += diff*diff;
	}
	return sqrt(distance);
}
