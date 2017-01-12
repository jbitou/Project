#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "distances.h"

double distance_Euclidean(ratingp v1, ratingp v2, int d) {
	int i;
	double distance = 0.0;
	/**For all dimensions of vectors v1 and v2**/
	for (i=0; i < d; i++) distance += (v1[i].rate - v2[i].rate) * (v1[i].rate - v2[i].rate);
	return sqrt(distance);
}

double distance_Cosine(ratingp v1, ratingp v2, int d) {
	int i;
	double inner = 0.0, normx = 0.0, normy = 0.0, cos;
	/**Inner product**/
	for (i=0; i < d; i++) {
		/**Numerator**/
		inner += v1[i].rate * v2[i].rate;	
		/**Euclidean norm x**/
		normx += v1[i].rate * v1[i].rate;
		/**Euclidean norm y**/
		normy += v2[i].rate * v2[i].rate;
	}
	normx = sqrt(normx);
	normy = sqrt(normy);
	cos = inner / (normx * normy);
	return (1 - cos);
}
