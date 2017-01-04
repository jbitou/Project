#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "inputProcessing.h"

int main (int argc, char **argv) {
	int i, j, T, r, numConform, N, size, k;
	double **data, v1, v2, v3;
	char read[12];
	FILE *fp, *fe;
	dinfo **alldistances;
	srand(time(NULL));
	if (command_processing(argc) == -1)	return -1;
	/**Open input file**/
	fp = fopen(argv[2],"r");
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
	/**Read number of conformations and N**/
	fscanf(fp,"%d\n",&numConform);
	fscanf(fp,"%d\n",&N);
	if (numConform > 100) k = 110;
	else k = 2;
	size = numConform*N;
	data = malloc(size*sizeof(double *));
	for (i=0; i < size; i++)	data[i] = malloc(3*sizeof(double));
	i = 0;
	/**Read points**/
	while (fscanf(fp,"%lf %lf %lf[^\n]",&v1,&v2,&v3) != EOF) {
		data[i][0] = v1;
		data[i][1] = v2;
		data[i][2] = v3;
		i++;
	}
	alldistances = get_all_distances(data,numConform,N);
	fprintf(fe,"r\tT\tk\tSilhouette\tTime\t\n");
	for (T=1; T <= 3; T++)  {
		r = N;
		experiment(data,alldistances,numConform,N,r,T,k,fe); 
		/*r = pow(N,1.5);
		if (r <= N*(N-1)/2) {
			experiment(data,alldistances,numConform,N,r,T,k,fe);
		}*/
	}
	//r = N*(N-1)/2;
	//experiment(data,alldistances,numConform,N,r,0,k,fe); 
	for (i=0; i < numConform; i++)	free(alldistances[i]);
	free(alldistances);
	alldistances = NULL;
	for (i=0; i < size; i++)	free(data[i]);
	free(data);
	data = NULL;
	fclose(fp);
	fclose(fe);
	return 0;
}
