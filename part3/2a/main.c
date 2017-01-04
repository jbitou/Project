#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "inputProcessing.h"
#include "silhouette.h"


int main (int argc, char **argv) {
	int i, j, numConform, N, size, k, bestk;
	double **data, v1, v2, v3, bestS;
	char read[12];
	FILE *fp, *fe;
	pcluster bestclusters;
	srand(time(NULL));
	if (command_processing(argc) == -1)	return -1;
	/**Open input file**/
	fp = fopen(argv[2],"r");
	if (fp == NULL) {
		perror("Error opening input_file");
		return -1;
	}
	/**Open output file**/
	fe = fopen("conform.dat","w+");
	if (fe == NULL) {
		perror("Error opening output_file");
		return -1;
	}
	fscanf(fp,"%d\n",&numConform);
	fscanf(fp,"%d\n",&N);
	if (numConform > 200) k = 100;
	else k = 2;
	size = numConform*N;
	data = malloc(size*sizeof(double *));
	for (i=0; i < size; i++)	data[i] = malloc(3*sizeof(double));
	i = 0;
	while (fscanf(fp,"%lf %lf %lf[^\n]",&v1,&v2,&v3) != EOF) {
		data[i][0] = v1;
		data[i][1] = v2;
		data[i][2] = v3;
		i++;
	}
	fclose(fp);
	/**Translate to common origin**/
	translation(data,numConform,N);
	/**Clustering**/
	bestclusters = k_clustering(data,numConform,N,k,&bestk,&bestS);
	/**Write in output file**/
	fprintf(fe,"k: %d\n",bestk);
	fprintf(fe,"s: %lf\n",bestS);
	for (i=0; i < bestk; i++) printndestroy_points(&(bestclusters[i].items),fe);
	free(bestclusters);
	bestclusters = NULL;	
	for (i=0; i < size; i++)	free(data[i]);
	free(data);
	data = NULL;
	fclose(fe);
	return 0;
}
