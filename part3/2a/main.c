#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inputProcessing.h"
#include "silhouette.h"
#define EXPERK 10

int main (int argc, char **argv) {
	FILE *fp, *fe;
	char read[12];
	int i, j, flag, numConform, N, size, k, bestk, counter = 0;
	double **data, v1, v2, v3, silhouette, previousS, bestS;
	pcluster clusters, bestclusters;
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
	fscanf(fp,"%s%d\n",read,&numConform);
	fscanf(fp,"%s%d\n",read,&N);
	if (numConform > 100) k = 10;
	else k = 4;
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
	for (i=0; i < size; i++)	printf("v1=%lf v2=%lf v3=%lf\n",data[i][0],data[i][1],data[i][2]);
	/**Translate to common origin**/
	translation(data,numConform,N);
	//distanceCRMSD(data,N,0,0);
	/**Clustering**/
	previousS = -1.1;
	while ((counter != -1) && (counter <= EXPERK)) {
		if (k == numConform / 2)	break;
		flag = 0;
		clusters = clustering(data,numConform,N,k);
		for (i=0; i < k; i++) {
			int n = points_length(clusters[i].items);
			if (n < 3) flag++;
		}
		if (flag == k) break;
		silhouette = compute_silhouette(clusters,data,N,numConform,k);
		printf("Silhouette = %lf for k = %d\n",silhouette,k);
		if (silhouette > previousS)	{
			if (counter != 0) {
				for (i=0; i < k; i++) destroy_points(&(bestclusters[i].items));
				free(bestclusters);
			}
			counter++;
			bestk = k;
			bestS = silhouette;
			bestclusters = malloc(k*sizeof(cluster));
			for (i=0; i < k; i++)	{
				bestclusters[i].items = NULL;	
				bestclusters[i].items = clone_points(clusters[i].items);
			}
		}
		else counter = -1;
		for (i=0; i < k; i++) destroy_points(&(clusters[i].items));
		free(clusters);
		previousS = silhouette;
		k++;
	}
	printf("bestk = %d\n",bestk);
	/**Write in output file**/
	fprintf(fe,"k: %d\n",bestk);
	fprintf(fe,"s: %lf\n",bestS);
	for (i=0; i < bestk; i++) printndestroy_points(&(bestclusters[i].items),fe);
	free(bestclusters);	
	for (i=0; i < size; i++)	free(data[i]);
	free(data);
	fclose(fp);
	fclose(fe);
	return 0;
}
