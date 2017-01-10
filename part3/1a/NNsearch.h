#include "hash.h"

double nnlsh_recommendation(user *, user *, int, int, int, int, FILE *, int);
nnrlist binary_repeated_nnsearch(hash_table *, ghashp *, user *, int, int, int, int, double, int);
ratinglist evaluate_rest_items(user, nnrlist, int);
ratinglist evaluate_items(user, nnrlist, int);
double find_normalizing_factor(nnrlist); 
double find_first_radius(user *, int, int, int);
