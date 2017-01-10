#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lapacke.h"
#include "cblas.h"
#include "distanceCRMSD.h"
#define EPSILON 0.000000001

double distanceCRMSD(double **data, int N, int conf1, int conf2) {
	int i, j, flag = 0;
	lapack_int info;
	double *xtrs, *Y, *C, *singular, *stat, *V, *Vtrs, *tempQ, *Q, *U, *X, *XQ, *XQ_Y, det, crmsd, sum = 0.0;
	/**Get X and Y matrices from data**/
	X = get_pointset(data,N,conf1);
	Y = get_pointset(data,N,conf2);
	xtrs = find_transpose(X,N,3);
	C = LAPACKE_malloc(9*sizeof(double));
	singular = LAPACKE_malloc(3*sizeof(double));
	stat = LAPACKE_malloc(6*sizeof(double));
	V = LAPACKE_malloc(9*sizeof(double));
	for (i=0; i < 9; i++)	{
		V[i] = 0.0;
		C[i] = 0.0;
	}
	/**Multiply X^T with Y**/	
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,3,3,N,1.0,xtrs,N,Y,3,0.0,C,3);
	/**SVD**/
	info = LAPACKE_dgesvj(LAPACK_ROW_MAJOR,'G','U','V',3,3,C,3,singular,0,V,3,stat);
	Vtrs = find_transpose(V,3,3);
	U = C;
	/**Check singular**/
	for (i=0; i < 3; i++) {
		if (singular[i] <= EPSILON) flag = 1;
	}
	if (flag == 0) {
		/**Multiply U with V^T**/
		Q = LAPACKE_malloc(9*sizeof(double));
		tempQ = LAPACKE_malloc(9*sizeof(double));
		cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,3,3,3,1.0,U,3,Vtrs,3,0.0,Q,3);
		/**Find the determinant**/
		for (i=0; i < 9; i++)	tempQ[i] = Q[i];
		det = find_det(tempQ);
		LAPACKE_free(tempQ);
		/**If det < 0 change Q**/
		if (det < 0) {
			U[2] *= -1;
			U[5] *= -1;
			U[8] *= -1;
			for (i=0; i < 9; i++)	Q[i] = 0.0;
			cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,3,3,3,1.0,U,3,Vtrs,3,0.0,Q,3);
		}
		/**Calculate XQ-Y**/
		XQ = LAPACKE_malloc(N*3*sizeof(double));
		XQ_Y = LAPACKE_malloc(N*3*sizeof(double));
		for (i=0; i < N*3; i++)	XQ[i] = 0.0;
		cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N,3,3,1.0,X,3,Q,3,0.0,XQ,3);
		for (i=0; i < N*3; i++) XQ_Y[i] = XQ[i]-Y[i];
		/*for (i=0; i < N*3; i++) sum += XQ_Y[i]*XQ_Y[i];
		printf("1.frobenius = %lf\n",sqrt(sum));*/
		/**Find frobenius norm**/
		crmsd = find_frobenius_norm(XQ_Y,N);
		/**Calculate c-RMSD**/
		crmsd /= sqrt(N);
		/**Free allocated memory**/
		LAPACKE_free(Q);	
		LAPACKE_free(XQ);	
		LAPACKE_free(XQ_Y);	
	}
	else crmsd = 0.1;
	/**Free allocated memory**/
	LAPACKE_free(xtrs);
	LAPACKE_free(Y);
	LAPACKE_free(X);
	LAPACKE_free(C);
	LAPACKE_free(Vtrs);
	LAPACKE_free(V);
	LAPACKE_free(singular);
	LAPACKE_free(stat);
	return crmsd;
}

double find_frobenius_norm(double *XQ_Y, int N) {
	int i = 0;
	double trace = 0.0, *R, *XQ_Ytrs;
	/**	Calculate M^T x M**/
	XQ_Ytrs = find_transpose(XQ_Y,N,3);
	R = LAPACKE_malloc(N*N*sizeof(double));
	for (i=0; i < N*N; i++)	R[i] = 0.0;
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N,N,3,1.0,XQ_Y,3,XQ_Ytrs,N,0.0,R,N);
	i = 0;
	while (i < N*N) {
		trace += R[i];
		i += N+1;
	}
	LAPACKE_free(XQ_Ytrs);
	LAPACKE_free(R);
	return sqrt(trace);
}

double *get_pointset(double **data, int N, int conf) {
	int i, j = 0, start;
	double *X;
	X = LAPACKE_malloc(N*3*sizeof(double));
	start = conf*N;
	for (i=start; i < start+N; i++) {
		X[j] = data[i][0];
		j++;
		X[j] = data[i][1];
		j++;
		X[j] = data[i][2];
		j++;
	}
	return X;
}

double find_det(double *Q) {
	int i, pos = -1, pos1 = -1;
	double det, *real, *imaginary, *rvectors, *lvectors;
	real = LAPACKE_malloc(3*sizeof(double));
	imaginary = LAPACKE_malloc(3*sizeof(double));
	LAPACKE_dgeev(LAPACK_ROW_MAJOR,'N','N',3,Q,3,real,imaginary,lvectors,3,rvectors,3);
	if (imaginary[0] != 0.0 || imaginary[1] != 0.0 || imaginary[2] != 0.0)	{
		for (i=0; i < 3; i++) {
			if (abs(real[i]-1.0) <= EPSILON) {
				pos = i;
				break;
			}
			pos1 = i;
		}
		if (pos1 == -1) det = (real[1]*real[1]+imaginary[1]*imaginary[1])*real[pos];
		else det = (real[pos1]*real[pos1]+imaginary[pos1]*imaginary[pos1])*real[pos];
	}
	else 	det = real[0] * real[1] * real[2];
	LAPACKE_free(real);
	LAPACKE_free(imaginary);
	return det;
}

double *find_transpose(double *matrix, int N, int M) {
	int i, j, z = 0;
	double *trs = LAPACKE_malloc(M*N*sizeof(double));
	for (i=0; i < M; i++) {
		for (j=0; j < N; j++) {
			trs[z] = matrix[i+j*M];
			z++;
		}
	}
	return trs;
}

double *find_centroid(double **data, int N, int conf) {
	int i, start;
	double *result;
	result = malloc(3*sizeof(double));
	start = conf*N;
	result[0] = result[1] = result[2] = 0;
	for (i=start; i < start+N; i++) {
		result[0] += data[i][0];
		result[1] += data[i][1];
		result[2] += data[i][2];
	}
	result[0] /= N;
	result[1] /= N;
	result[2] /= N;
	return result;
}

void translation(double **data, int numConform, int N) {
	int i, j, start;
	double *xc;
	for (i=0; i < numConform; i++) {
		xc = find_centroid(data,N,i);
		start = i*N;
		for (j=start; j < start+N; j++) {
			data[j][0] -= xc[0];
			data[j][1] -= xc[1];
			data[j][2] -= xc[2];
		}
		free(xc);
	}
}

