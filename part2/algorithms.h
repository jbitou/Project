#include "silhouette.h"

pcluster CLARA(pinfo, int, ghashp *, int, void **, FILE *, hash_table, int);
centroid *PAM(hash_table, ghashp *, pinfo, void *, int, int);
pcluster IAU1(hash_table *, ghashp *, pinfo, void *, int, int, int);
pcluster IAU2(hash_table *, ghashp *, pinfo, void *, int, int, int);
int *choose_clara_sample(int, int);
