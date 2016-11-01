#include <stdio.h>
#include <stdlib.h>
#include <CUnit/Basic.h>
#include <inttypes.h>
#include "distances.h"

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

void test_distance_Hamming()
{
	uint64_t num1, num2;
	num1 = 20;
	num2 = 30;
	CU_ASSERT(2 == distance_Hamming(num1,num2));
	num1 = 1896;
	num2 = 10523;
	CU_ASSERT(9 == distance_Hamming(num1,num2));
}

void test_distance_Euclidean()
{
	double *v1,*v2,epsilon = 0.000000001,distance;
	int d = 2,i;
	v1 = malloc(d*sizeof(double));
	v1[0] = 86.569874;
	v1[1] = 92.385705;
	v2 = malloc(d*sizeof(double));
	v2[0] = 69.736987;
	v2[1] = -5.333366;
	distance = distance_Euclidean(v1,v2,d);
	CU_ASSERT(fabs(distance - 99.15827208) < epsilon);
	free(v1);
	free(v2);
}

void test_distance_Cosine()
{
	double *v1,*v2,epsilon = 0.000000001,distance;
	int d = 2,i;
	v1 = malloc(d*sizeof(double));
	v1[0] = 86.569874;
	v1[1] = 92.385705;
	v2 = malloc(d*sizeof(double));
	v2[0] = 69.736987;
	v2[1] = -5.333366;
	distance = distance_Cosine(v1,v2,d);
	CU_ASSERT(fabs(distance - 0.3738690549) < epsilon);
	free(v1);
	free(v2);
}

int main (void) 
{
    CU_pSuite pSuite = NULL;

    /* Initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    /* Add a suite to the registry */
    pSuite = CU_add_suite("DistancesTest", init_suite, clean_suite);
    if (NULL == pSuite) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    /* Add the tests to the suite */
    if ((NULL == CU_add_test(pSuite, "test_distance_Hamming", test_distance_Hamming))
     || (NULL == CU_add_test(pSuite, "test_distance_Euclidean", test_distance_Euclidean))
     || (NULL == CU_add_test(pSuite, "test_distance_Cosine", test_distance_Cosine)))
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
