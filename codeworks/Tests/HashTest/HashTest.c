#include <stdio.h>
#include <stdlib.h>
#include <CUnit/Basic.h>
#include <inttypes.h>
#include "hash.h"

/*
 * CUnit Test Suite
*/

int init_suite(void) 
{
    return 0;
}

int clean_suite(void) 
{
    return 0;
}

void test_hash_func_Ham() 
{
	int i;
	ghashp g = malloc(4*sizeof(ghash));
	for(i=0; i < 4; i++)
		g[i].t = 3 - i;	
 	CU_ASSERT(8 == hash_func_Ham(g,"00011110101110001100",4));
 	free(g);
}

void test_hash_func_Eucl() 
{
	int i,j;
	ghashp g = malloc(4*sizeof(ghash));
	for(i=0; i < 4; i++)
	{
		g[i].v = malloc(4*sizeof(double));	
		for(j=0; j < 4; j++)
			g[i].v[j] = 0.3*i + 0.8*j;
		g[i].t = 3 - i;	
		g[i].r = i;	
	}
	double *p = malloc(4*sizeof(double));	
	for(i=0; i < 4; i++)
		p[i] = i + 0.1234;
 	CU_ASSERT(23 == hash_func_Eucl(g,p,4,4));
 	for (i=0; i < 4; i++) free(g[i].v);
 	free(g);
 	free(p);
}

void test_hash_func_Cos() 
{
	int i,j;
	ghashp g = malloc(4*sizeof(ghash));
	for(i=0; i < 4; i++)
	{
		g[i].v = malloc(4*sizeof(double));	
		for(j=0; j < 4; j++)
			g[i].v[j] = 0.3*i + 0.8*j;	
	}
	double *p = malloc(4*sizeof(double));	
	for(i=0; i < 4; i++)
		p[i] = i + 0.1234;
	CU_ASSERT(15 == hash_func_Cos(g,p,4,4));
	for(i=0; i < 4; i++)
	{	
		for(j=0; j < 4; j++)
		{
			if (i == 0)	g[i].v[j] =  (-0.8)*j;
			else g[i].v[j] = 0.3*i + 0.8*j;	
		}
	}
	CU_ASSERT(7 == hash_func_Cos(g,p,4,4));
 	for (i=0; i < 4; i++) free(g[i].v);
 	free(g);
 	free(p);
}

void test_hash_func_Matrix() 
{
	int i,j,counter;
	ghashp g = malloc(4*sizeof(ghash));
	for(i=0; i < 4; i++)
	{
		g[i].t = 3 - i;
		g[i].r = i + 1;
		g[i].t1 = 2 * i;
	}
	int ** distances = malloc(4*sizeof(int*));
	counter = 4;
	for(i=0; i < 4; i++)
	{
		distances[i] = 	malloc(counter*sizeof(int));
		counter--;
	}
	counter = 4;
	for(i=0; i < 4; i++)
	{	
		for(j=0; j < counter; j++)
			distances[i][j] = i + j;
		counter--;	
	}
	CU_ASSERT(8 == hash_func_Matrix(g,2,distances,4,4));
	for (i=0; i < 4; i++) free(distances[i]);
 	free(distances);
 	free(g);
}

int main (void) 
{
    CU_pSuite pSuite = NULL;

    /* Initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    /* Add a suite to the registry */
    pSuite = CU_add_suite("HashTest", init_suite, clean_suite);
    if (NULL == pSuite) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    /* Add the tests to the suite */
    if ((NULL == CU_add_test(pSuite, "test_hash_func_Ham", test_hash_func_Ham))
     || (NULL == CU_add_test(pSuite, "test_hash_func_Eucl", test_hash_func_Eucl))
     || (NULL == CU_add_test(pSuite, "test_hash_func_Cos", test_hash_func_Cos))
     || (NULL == CU_add_test(pSuite, "test_hash_func_Matrix", test_hash_func_Matrix)))
    {
        CU_cleanup_registry();
        return CU_get_error();
    }

    /* Run all tests using the CUnit Basic interface */
    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();
    return  CU_get_error();
}
