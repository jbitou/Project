#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structure.h"

centroid *vector_pam_update(pcluster clusters, centroid *centroids, double J, double **distances, pinfo info, int n) {
	int i, j, z, id, m, ci, t, flag, exists;
	double SDj, prevJ, J1, tdistance, *arr, *qdata;
	pointp temp,temp1;
	centroid *newcenter = malloc(info->k*sizeof(centroid));
	/**For each centroid**/
	for (i=0; i < info->k; i++) {
		m = (int)(intptr_t) centroids[i].center;
		newcenter[i].info = malloc(n*sizeof(double));
		newcenter[i].center = (void *)(intptr_t)-1;
		/**Swap with all non-centroids**/
		for (j=0; j < info->k; j++) {	
			temp = clusters[j].items;
			/**We swap  m with t**/
			while (temp != NULL) {
				prevJ = J1;
				SDj = flag = 0; 
				t = temp->position;
				/**If t is a centroid break**/
				for (z=0; z < info->k; z++) {
					if (t == (int)(intptr_t)centroids[z].center) {
						flag = 1;
						break;
					}
				}
				if (flag == 1) {
					temp = temp->next;
					continue;
				}
				/**Run for all items**/
				for (z=0; z < info->k; z++) {	
					temp1 = clusters[z].items;
					ci = (int)(intptr_t)clusters[z].center.center;
					/**For each item in the cluster j-if we swap it with m**/
					while (temp1 != NULL) {
						id = temp1->position;
						if (id < t)  		tdistance = distances[id][t-id-1];
						else if (id > t) 	tdistance = distances[t][id-t-1];
						if (ci == m)  {
							/**If dist(i,t) > dist(i,c')**/
							if (tdistance > temp1->secdistance)	SDj += temp1->secdistance - temp1->mindistance;
							else  SDj += tdistance - temp1->mindistance;
						}
						else {
							/**if dist(i,t) < dist(i,c(i))**/
							if (tdistance < temp1->mindistance)	SDj += tdistance - temp1->mindistance;
						}
						temp1 = temp1->next;
					}
				}
				J1 = J + SDj;
				if (J1 < prevJ) {
					exists = 0;
					for (z=0; z < info->k; z++) {
						if (newcenter[z].center == (void *)(intptr_t)t) {
							exists = 1;
							break;
						}
					}
					if (exists == 0) {
						qdata = malloc((n)*sizeof(double));
						for (z=0; z < t; z++) qdata[z] = distances[z][t-z-1];
						qdata[t] = 0;
						for (z=(t+1); z < n; z++) qdata[z] = distances[t][z-t-1];
						/**Insert medoid's info to newcenter**/
						newcenter[i].center = (void *)(intptr_t)t;
						arr = (double *)newcenter[i].info;
						for (z=0; z < n; z++) arr[z] = qdata[z];
						free(qdata);
					}
				}
				temp = temp->next;
			}
		}
	}
	/**Change centroid**/
	for (i=0; i < info->k; i++) {
		if (newcenter[i].center == (void *)(intptr_t)-1) continue;
		centroids[i].center = newcenter[i].center;
		arr = (double *)centroids[i].info;
		for (z=0; z < n; z++) arr[z] = ((double *)newcenter[i].info)[z];
	}
	for (i=0; i < info->k; i++) free(newcenter[i].info);
	free(newcenter);
	return centroids;
}

centroid *matrix_pam_update(pcluster clusters, centroid *centroids, double J, int **distances, pinfo info, int n) {
	int i, j, z, id, m, ci, t, flag, tdistance, *arr, *qdata, exists;
	double SDj, prevJ, J1;
	pointp temp,temp1;
	centroid *newcenter = malloc(info->k*sizeof(centroid));
	/**For each centroid**/
	for (i=0; i < info->k; i++) {
		m = (int)(intptr_t) centroids[i].center;
		newcenter[i].info = malloc(n*sizeof(int));
		newcenter[i].center = (void *)(intptr_t)-1;
		/**Swap with all non-centroids**/
		for (j=0; j < info->k; j++) {	
			temp = clusters[j].items;
			/**We swap  m with t**/
			while (temp != NULL) {
				prevJ = J1;
				SDj = flag = 0; 
				t = temp->position;
				/**If t is a centroid break**/
				for (z=0; z < info->k; z++) {
					if (t == (int)(intptr_t)centroids[z].center) {
						flag = 1;
						break;
					}
				}
				if (flag == 1) {
					temp = temp->next;
					continue;
				}
				/**Run for all items**/
				for (z=0; z < info->k; z++) {	
					temp1 = clusters[z].items;
					ci = (int)(intptr_t)clusters[z].center.center;
					/**For each item in the cluster j-if we swap it with m**/
					while (temp1 != NULL) {
						id = temp1->position;
						if (id < t)  		tdistance = distances[id][t-id-1];
						else if (id > t) 	tdistance = distances[t][id-t-1];
						if (ci == m)  {
							/**If dist(i,t) > dist(i,c')**/
							if (tdistance > temp1->secdistance)	SDj += temp1->secdistance - temp1->mindistance;
							else  SDj += tdistance - temp1->mindistance;
						}
						else {
							/**if dist(i,t) < dist(i,c(i))**/
							if (tdistance < temp1->mindistance)	SDj += tdistance - temp1->mindistance;
						}
						temp1 = temp1->next;
					}
				}
				J1 = J + SDj;
				if (J1 < prevJ) {
					exists = 0;
					for (z=0; z < info->k; z++) {
						if (newcenter[z].center == (void *)(intptr_t)t) {
							exists = 1;
							break;
						}
					}
					if (exists == 0) {
						qdata = malloc((n)*sizeof(int));
						for (z=0; z < t; z++) qdata[z] = distances[z][t-z-1];
						qdata[t] = 0;
						for (z=(t+1); z < n; z++) qdata[z] = distances[t][z-t-1];
						/**Insert medoid's info to newcenter**/
						newcenter[i].center = (void *)(intptr_t)t;
						arr = (int *)newcenter[i].info;
						for (z=0; z < n; z++) arr[z] = qdata[z];
						free(qdata);
					}
				}
				temp = temp->next;
			}
		}
	}
	/**Change centroid**/
	for (i=0; i < info->k; i++) {
		if (newcenter[i].center == (void *)(intptr_t)-1) continue;
		centroids[i].center = newcenter[i].center;
		arr = (int *)centroids[i].info;
		for (z=0; z < n; z++) arr[z] = ((int *)newcenter[i].info)[z];
	}
	for (i=0; i < info->k; i++) free(newcenter[i].info);
	free(newcenter);
	return centroids;
}


centroid *matrix_update_alaloyds(pcluster clusters, centroid *centroids, double J, int **distances, pinfo info) {
	pointp medoid, temp, delete;
	int i, j, z, s, ci, m, id, id1, id2, tdistance, *qdata, jump, *arr;
	double  SDj, J1;
	centroid newcenter;
	/**For each cluster**/
	for (i=0; i < info->k; i++) {
		SDj = 0;
		medoid = matrix_calculate_medoid(clusters[i].items,distances);
		if (medoid == NULL) continue;
		id1 = medoid->position;
		/**Create distances of medoid with all the items**/
		qdata = malloc((info->N)*sizeof(int));
		for (z=0; z < id1; z++) qdata[z] = distances[z][id1-z-1];
		qdata[id1] = 0;
		for (z=(id1+1); z < info->N; z++) qdata[z] = distances[id1][z-id1-1];
		/**Insert medoid's info to newcenter**/
		newcenter.center = (void *)(intptr_t)id1;
		newcenter.info = malloc(info->N*sizeof(int));
		arr = (int *)newcenter.info;
		for (z=0; z < info->N; z++) arr[z] = qdata[z];
		free(qdata);	
		/**Store centroid where we are**/
		m = (int)(intptr_t)clusters[i].center.center;
		/**Check all items using clusters**/
		for (j=0; j < info->k; j++) {
			ci = (int)(intptr_t)clusters[j].center.center;		
			temp = clusters[j].items;
			/**For each item in the cluster j**/
			while (temp != NULL) {
				id2 = temp->position;
				if (id1 == id2) {
					temp = temp->next;
					continue;
				}
				if (id1 < id2)  tdistance = distances[id1][id2-id1-1];
				else if (id1 > id2)  tdistance = distances[id2][id1-id2-1];
				if (ci == m)  {
					/**If dist(i,t) > dist(i,c')**/
					if (tdistance > temp->secdistance)	SDj += temp->secdistance - temp->mindistance;
					else  SDj += tdistance - temp->mindistance;
				}
				else {
					/**if dist(i,t) < dist(i,c(i))**/
					if (tdistance < temp->mindistance)	SDj += tdistance - temp->mindistance;
				}
				temp = temp->next;
			}
		}
		J1 = J + SDj;
		/**If J' < J, then swap centroid with medoid**/
		if (J1 < J) {
			centroids[i].center = newcenter.center;
			arr = (int *)centroids[i].info;
			for (z=0; z < info->N; z++) arr[z] = ((int *)newcenter.info)[z];

		}
		free(newcenter.info);
	}
	return centroids;
}

centroid *vector_update_alaloyds(pcluster clusters, centroid *centroids, double J, double **distances, pinfo info) {
	pointp medoid, temp, delete;
	int i, j, z, s, ci, m, id, id1, id2, jump;
	double tdistance, SDj, J1, *qdata, *arr;
	centroid newcenter;
	/**For each cluster**/
	for (i=0; i < info->k; i++) {
		SDj = 0;
		medoid = vector_calculate_medoid(clusters[i].items,distances);
		if (medoid == NULL) continue;
		id1 = medoid->position;
		/**Create distances of medoid with all the items**/
		qdata = malloc((info->N)*sizeof(double));
		for (z=0; z < id1; z++) qdata[z] = distances[z][id1-z-1];
		qdata[id1] = 0;
		for (z=(id1+1); z < info->N; z++) qdata[z] = distances[id1][z-id1-1];
		/**Insert medoid's info to newcenter**/
		newcenter.center = (void *)(intptr_t)id1;
		newcenter.info = malloc(info->N*sizeof(double));
		arr = (double *)newcenter.info;
		for (z=0; z < info->N; z++) arr[z] = qdata[z];
		free(qdata);	
		/**Store centroid where we are**/
		m = (int)(intptr_t)clusters[i].center.center;
		/**Check all items using clusters**/
		for (j=0; j < info->k; j++) {
			ci = (int)(intptr_t)clusters[j].center.center;		
			temp = clusters[j].items;
			/**For each item in the cluster j**/
			while (temp != NULL) {
				id2 = temp->position;
				if (id1 == id2) {
					temp = temp->next;
					continue;
				}
				if (id1 < id2)  tdistance = distances[id1][id2-id1-1];
				else if (id1 > id2)  tdistance = distances[id2][id1-id2-1];
				if (ci == m)  {
					/**If dist(i,t) > dist(i,c')**/
					if (tdistance > temp->secdistance)	SDj += temp->secdistance - temp->mindistance;
					else  SDj += tdistance - temp->mindistance;
				}
				else {
					/**if dist(i,t) < dist(i,c(i))**/
					if (tdistance < temp->mindistance)	SDj += tdistance - temp->mindistance;
				}
				temp = temp->next;
			}
		}
		J1 = J + SDj;
		/**If J' < J, then swap centroid with medoid**/
		if (J1 < J) {
			centroids[i].center = newcenter.center;
			arr = (double *)centroids[i].info;
			for (z=0; z < info->N; z++) arr[z] = ((double *)newcenter.info)[z];

		}
		free(newcenter.info);
	}
	return centroids;
}

cpair sortPairs(cpair array, int N) {
	int i, j, m1, m2;
	pair a;
	for (i = 0; i < N; i++) {
		for (j = i + 1; j < N; j++) {
			m1 = (int)(intptr_t) array[i].m.center;
			m2 = (int)(intptr_t) array[j].m.center;
			if (m1 > m2) {
				a =  array[i];
                array[i] = array[j];
                array[j] = a;
            }
        }
	}
	return array;
}

centroid *matrix_update_clarans(pcluster clusters, centroid *centroids, cpair pairs, int **distances, int J, pinfo info) {
	int i, j, z, y, SDj, J1, ci, m, t, id, pos, tdistance, *arr, newm, in,change;
	pointp temp;
	centroid newcenter;
	pairs = sortPairs(pairs, info->fraction);	
	i = 0;
	while (i < info->fraction) {
		SDj = in = change = 0;
		m = (int)(intptr_t)pairs[i].m.center;
		newcenter.info = malloc((info->N)*sizeof(int));
		for (z=i; z < info->fraction; z++) {
			newm = (int)(intptr_t)pairs[z].m.center;
			if (newm != m) {
				in = 1;
				break;
			}
			t = (int)(intptr_t)pairs[z].t.center;
			/**Check all items using clusters**/
			for (j=0; j < info->k; j++) {
				ci = (int)(intptr_t)clusters[j].center.center;		
				temp = clusters[j].items;
				/**For each item in the cluster j**/
				while (temp != NULL) {
					id = temp->position;
					if (t == id) {
						temp = temp->next;
						continue;
					}
					if (t < id)  tdistance = distances[t][id-t-1];
					else if (t > id)  tdistance = distances[id][t-id-1];
					if (ci == m)  {
						/**If dist(i,t) > dist(i,c')**/
						if (tdistance > temp->secdistance)	SDj += temp->secdistance - temp->mindistance;
						else  SDj += tdistance - temp->mindistance;
					}
					else {
						/**if dist(i,t) < dist(i,c(i))**/
						if (tdistance < temp->mindistance)	SDj += tdistance - temp->mindistance;
					}
					temp = temp->next;
				}
			}	
			J1 = J + SDj;
			if (J1 < J) {
				change = 1;
				newcenter.center =  pairs[z].t.center;
				arr = (int *)newcenter.info;
				for (y=0; y < info->N; y++) arr[y] = ((int *)pairs[z].t.info)[y];
			}
		}
		if (change == 1) {
			for (y=0; y < info->k; y++) {
				if (m == (int)(intptr_t)centroids[y].center) {
					pos = y;
					break;
				}
			}
			centroids[pos].center = newcenter.center;
			arr = (int *)centroids[pos].info;
			for (y=0; y < info->N; y++) arr[y] = ((int *)newcenter.info)[y];
		}
		free(newcenter.info);
		if (!in) i++;
		else 	 i = z;
	}
	return centroids;
}

centroid *vector_update_clarans(pcluster clusters, centroid *centroids, cpair pairs, double **distances, double J, pinfo info) {
	int i, j, z, y, ci, m, t, id, pos, newm, in,change;
	double tdistance, *arr, SDj, J1;
	pointp temp;
	centroid newcenter;
	pairs = sortPairs(pairs, info->fraction);
	i = 0;
	while (i < info->fraction) {
		SDj = in = change = 0;
		m = (int)(intptr_t)pairs[i].m.center;
		newcenter.info = malloc((info->N)*sizeof(double));
		for (z=i; z < info->fraction; z++) {
			newm = (int)(intptr_t)pairs[z].m.center;
			if (newm != m) {
				in = 1;
				break;
			}
			t = (int)(intptr_t)pairs[z].t.center;
			/**Check all items using clusters**/
			for (j=0; j < info->k; j++) {
				ci = (int)(intptr_t)clusters[j].center.center;		
				temp = clusters[j].items;
				/**For each item in the cluster j**/
				while (temp != NULL) {
					id = temp->position;
					if (t == id) {
						temp = temp->next;
						continue;
					}
					if (t < id)  tdistance = distances[t][id-t-1];
					else if (t > id)  tdistance = distances[id][t-id-1];
					if (ci == m)  {
						/**If dist(i,t) > dist(i,c')**/
						if (tdistance > temp->secdistance)	SDj += temp->secdistance - temp->mindistance;
						else  SDj += tdistance - temp->mindistance;
					}
					else {
						/**if dist(i,t) < dist(i,c(i))**/
						if (tdistance < temp->mindistance)	SDj += tdistance - temp->mindistance;
					}
					temp = temp->next;
				}
			}	
			J1 = J + SDj;
			if (J1 < J) {
				change = 1;
				newcenter.center =  pairs[z].t.center;
				arr = (double *)newcenter.info;
				for (y=0; y < info->N; y++) arr[y] = ((double *)pairs[z].t.info)[y];
			}
		}
		if (change == 1) {
			for (y=0; y < info->k; y++) {
				if (m == (int)(intptr_t)centroids[y].center) {
					pos = y;
					break;
				}
			}
			centroids[pos].center = newcenter.center;
			arr = (double *)centroids[pos].info;
			for (y=0; y < info->N; y++) arr[y] = ((double *)newcenter.info)[y];
		}
		free(newcenter.info);
		if (!in) i++;
		else 	 i = z;
	}
	return centroids;
}

pointp matrix_calculate_medoid(pointp items, int **distances) {
	int min, distance, sum, id1, id2;
	pointp temp, curr, first, medoid;
	temp = first = items;
	if (temp == NULL) return NULL;
	min = 0;
	/**For each item calculate total distance from first item**/
	id1 = first->position;
	medoid = temp;
	while (temp != NULL) {
		id2 = temp->position;
		if (id1 < id2)  distance = distances[id1][id2-id1-1];
		else if (id2 < id1) distance = distances[id2][id1-id2-1];
		else 	distance = 0;
		min += distance;
		temp = temp->next;
	}
	/**For each item, beginning from second**/
	temp = first->next;
	while (temp != NULL) {
		id1 = temp->position;
		sum = 0;
		curr = items;
		/**For each item calculate sum**/
		while (curr != NULL) {
			id2 = curr->position;
			if (id1 < id2)  distance = distances[id1][id2-id1-1];
			else if (id2 < id1) distance = distances[id2][id1-id2-1];
			else 	distance = 0;
			sum += distance;
			curr = curr->next;
		}
		if (sum < min) {
			min = sum;
			medoid = temp;
		}
		temp = temp->next;
	}
	return medoid;
}

pointp vector_calculate_medoid(pointp items, double **distances) {
	int id1, id2;
	double min, distance, sum;
	pointp temp, curr, first, medoid;
	temp = first = items;
	if (temp == NULL) return NULL;
	min = 0;
	/**For each item calculate total distance from first item**/
	id1 = first->position;
	medoid = temp;
	while (temp != NULL) {
		id2 = temp->position;
		if (id1 < id2)  distance = distances[id1][id2-id1-1];
		else if (id2 < id1) distance = distances[id2][id1-id2-1];
		else 	distance = 0;
		min += distance;
		temp = temp->next;
	}
	/**For each item, beginning from second**/
	temp = first->next;
	while (temp != NULL) {
		id1 = temp->position;
		sum = 0;
		curr = items;
		/**For each item calculate sum**/
		while (curr != NULL) {
			id2 = curr->position;
			if (id1 < id2)  distance = distances[id1][id2-id1-1];
			else if (id2 < id1) distance = distances[id2][id1-id2-1];
			else 	distance = 0;
			sum += distance;
			curr = curr->next;
		}
		if (sum < min) {
			min = sum;
			medoid = temp;
		}
		temp = temp->next;
	}
	return medoid;
}


int compare_centroids(centroid *centroids, centroid *previous, int k) {
	int i, diff = 0, id1, id2;
	for (i=0; i < k; i++) {
		id1 = (int)(intptr_t)centroids[i].center;	
		id2 = (int)(intptr_t)previous[i].center;	
		if (id1 != id2)	diff++;
	}
	return diff;
}

cpair matrix_select_pairs(centroid *centroids, int **distances, pinfo info) {
	int i = 0, j, z, x, index, item, again, *arr, *qdata;
	cpair pairs;
	pairs = malloc(info->fraction*sizeof(pair));
	/**|Q| pairs**/
	while (i < info->fraction) {
		/**Pick a uniformly distributed integer x**/
		x = (rand() / (RAND_MAX + 1.0)) * (info->k*info->N);
		index = mod(x,info->k);
		item = x / info->k;
		if (item != 0)	item--;
		again = 0;
		for (j=0; j < info->k; j++) {
			if (item  == (int)(intptr_t)centroids[j].center) {
				again = 1;
				break;
			}
		}
		/**If t isn't one of centroids**/
		if (again == 0) {
			/**Make m**/
			pairs[i].m.info = malloc(info->N*sizeof(int));
			arr = (int *)pairs[i].m.info;
			for (z=0; z < info->N; z++) arr[z] = ((int *)centroids[index].info)[z];
			pairs[i].m.center = centroids[index].center;
			/**Make t**/
			pairs[i].t.center = (void *)(intptr_t)item;
			qdata = malloc((info->N)*sizeof(int));
			for (z=0; z < item; z++) qdata[z] = distances[z][item-z-1];
			qdata[item] = 0;
			for (z=(item+1); z < info->N; z++) qdata[z] = distances[item][z-item-1];
			pairs[i].t.info = malloc(info->N*sizeof(int));
			arr = (int *)pairs[i].t.info;
			for (z=0; z < info->N; z++) arr[z] = qdata[z];
			free(qdata);
			i++;
		}
	}
	return pairs;
}

cpair vector_select_pairs(centroid *centroids, double **distances, pinfo info) {
	int i = 0, j, z, x, index, item, again;
	double *arr, *qdata;
	cpair pairs;
	pairs = malloc(info->fraction*sizeof(pair));
	/**|Q| pairs**/
	while (i < info->fraction) {
		/**Pick a uniformly distributed integer x**/
		x = (rand() / (RAND_MAX + 1.0)) * (info->k*info->N);
		index = mod(x,info->k);
		item = x / info->k;
		if (item != 0)	item--;
		again = 0;
		for (j=0; j < info->k; j++) {
			if (item  == (int)(intptr_t)centroids[j].center) {
				again = 1;
				break;
			}
		}
		/**If t isn't one of centroids**/
		if (again == 0) {
			/**Make m**/
			pairs[i].m.info = malloc(info->N*sizeof(double));
			arr = (double *)pairs[i].m.info;
			for (z=0; z < info->N; z++) arr[z] = ((double *)centroids[index].info)[z];
			pairs[i].m.center = centroids[index].center;
			/**Make t**/
			pairs[i].t.center = (void *)(intptr_t)item;
			qdata = malloc((info->N)*sizeof(double));
			for (z=0; z < item; z++) qdata[z] = distances[z][item-z-1];
			qdata[item] = 0;
			for (z=(item+1); z < info->N; z++) qdata[z] = distances[item][z-item-1];
			pairs[i].t.info = malloc(info->N*sizeof(double));
			arr = (double *)pairs[i].t.info;
			for (z=0; z < info->N; z++) arr[z] = qdata[z];
			free(qdata);
			i++;
		}
	}
	return pairs;
}

