#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "distances.h"
  
/*Returns the different bits between two points*/
int distance_Hamming(uint64_t num1, uint64_t num2)		
{
	int diff = 0;
	uint64_t mask = 1;
	uint64_t num = num1^num2;
	while (num > 0)
	{
		if ((num & mask) == 1)
			diff++;
		num = num >> 1; 
	}
	return diff;
}

/*Returns the euclidean difference between two vectors*/
double distance_Euclidean(double *v1, double *v2, int d)		
{
	int i;
	double distance = 0.0, diff = 0.0 ;
	for (i=0; i < d; i++)
	{
		diff = v1[i] - v2[i];
		distance += pow(diff,2);
	}
	distance = pow(distance,0.5);	//Not necessary for the comparison
	return distance;
}

/*Returns the distance (1-cosine similarity) between two vectors*/
double distance_Cosine(double *v1, double *v2, int d)	
{
	int i;
	double inner = 0.0, normx = 0.0, normy = 0.0, cos = 0.0;
	for (i=0; i < d; i++)	//Inner product
	{
		inner += v1[i]*v2[i];	//Numerator
		normx += pow(v1[i],2); 	//Euclidean norm x
		normy += pow(v2[i],2);	//Euclidean norm y
	}
	normx = pow(normx,0.5);
	normy = pow(normy,0.5);
	cos = inner / (normx * normy);
	return (1 - cos);
}
