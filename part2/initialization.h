#include "inputProcessing.h" 

typedef struct park_jun_info {
	int index;
	double v;
}pj_info;

int binarySearch(int, int, int *);
pj_info *sortArray(pj_info *, int);
int *matrix_init_kmedoids(int **, pinfo, int);
int *matrix_init_concentrate(int **, pinfo, int);
