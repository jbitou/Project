#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "inputProcessing.h"
#include "distanceDRMSD.h"
#define EXPERK 10

int main (int argc, char **argv) {
	FILE *fp, *fe;
	char read[12];
	int i, j, flag, input, T, r, numConform, N, size, k, bestk, counter = 0;
	double **data, v1, v2, v3, silhouette, previousS, bestS;
	dinfo *distances;
	srand(time(NULL));
	//pcluster clusters, bestclusters;
	if (command_processing(argc) == -1)	return -1;
	for (i=1; i < (argc - 1); i+=2) {
		if (strcmp(argv[i],"-d") == 0)	input = i+1;
		else if (strcmp(argv[i],"-T") == 0)	T = atoi(argv[i+1]);
	}
	while ((T < 1) || (T > 3)) {
		printf("For r smallest distances, choose 1.\nFor r largest distances, choose 2.\nFor randomly chosen distances, choose 3.\n");
		scanf("%d",&T);
	}
	/**Open input file**/
	fp = fopen(argv[input],"r");
	if (fp == NULL) {
		perror("Error opening input_file");
		return -1;
	}
	/**Open output file**/
	fe = fopen("experim.dat","w+");
	if (fe == NULL) {
		perror("Error opening output_file");
		return -1;
	}
	fscanf(fp,"%s%d\n",read,&numConform);
	fscanf(fp,"%s%d\n",read,&N);
	if (numConform > 100) k = 10;
	else k = 2;
	size = numConform*N;
	printf("numconform=%d and N=%d\n",numConform,N);
	data = malloc(size*sizeof(double *));
	for (i=0; i < size; i++)	data[i] = malloc(3*sizeof(double));
	i = 0;
	while (fscanf(fp,"%lf %lf %lf[^\n]",&v1,&v2,&v3) != EOF) {
		data[i][0] = v1;
		data[i][1] = v2;
		data[i][2] = v3;
		i++;
	}
	r = N;
	distances = create_distances(data,numConform,N,r,T);
	//for (i=0; i < size; i++)	printf("[%lf,%lf,%lf]\n",data[i][0],data[i][1],data[i][2]);
	for (i=0; i < size; i++)	free(data[i]);
	free(distances);
	free(data);
	fclose(fp);
	fclose(fe);
	return 0;
}
