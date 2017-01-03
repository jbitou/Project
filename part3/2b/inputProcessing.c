#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inputProcessing.h"

int command_processing(int argc) {
	/**Input file has to be given after -d recogniser**/
	if (argc > 3) {
		printf("Too many arguments. Try again.\n");
		return -1;
	}
	else if (argc < 3) {
		printf("Too few arguments. Try again.\n");
		return -1;
	}
	return 0;	
}


int clean_memory(double ** vectors, dinfo *distances, pcluster clusters, int numConform, int k) {
	int i;
	if (vectors == NULL || distances == NULL || clusters  ==  NULL) return -1;
	/**Free allocated memory**/
	for (i=0; i < numConform; i++)	free(vectors[i]);
	free(vectors);
	vectors = NULL;
	free(distances);
	distances = NULL;
	for (i=0; i < k; i++) {
		free(clusters[i].center.vector);
		destroy_points(&(clusters[i].items));
	}
	free(clusters);
	clusters = NULL;
	return 0;
}
