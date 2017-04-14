/*
 * get_list.c
 *
 *  Created on: Apr 14, 2017
 *      Author: fafa
 */

#include "libgtftk.h"

extern INDEX_ID *index_gtf(GTF_DATA *gtf_data, char *key);
extern COLUMN **column;

int n;
ROW_LIST *row_list;
TTEXT *vret;

static void action_row_list(const void *nodep, const VISIT which, const int depth) {
	ROW_LIST *datap;
	char tmp[100];

	switch (which) {
		case preorder:
			break;
		case postorder:
		case leaf:
			datap = *((ROW_LIST **)nodep);
			vret->data = (char ***)realloc(vret->data, (vret->size + 1) * sizeof(char **));
			vret->data[vret->size] = (char **)calloc(2, sizeof(char *));
			sprintf(tmp, "%d", datap->nb_row);
			vret->data[vret->size][0] = strdup(tmp);
			vret->data[vret->size][1] = strdup(datap->token);
			vret->size++;
			break;
		case endorder:
			break;
	}
}

__attribute__ ((visibility ("default")))
TTEXT *get_feature_list(GTF_DATA *gtf_data) {
	/*
	 * reserve memory for the TTEXT structure to return
	 */
	TTEXT *ret = (TTEXT *)calloc(1, sizeof(TTEXT));

	/*
	 * reserve memory for the local ROW_LIST
	 */
	row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));

	/*
	 * indexing the GTF_DATA with the feature column
	 */
	INDEX_ID *index_id = index_gtf(gtf_data, "feature");

	/*
	 * tree browsing of the feature index
	 */
	twalk(column[index_id->column]->index[index_id->index_rank].data, action_row_list);

	//ret->size = vret->size;
	//ret->data = (char ***)calloc(ret->size, sizeof(char **));
	int i;
	/*for (i = 0; i < ret->size; i++)	{
		ret->data[i] = (char **)calloc(2, sizeof(char *));
		ret->data[i][0] = vret->data[i][0];
		ret->data[i][1] = vret->data[i][1];
		//free(vret->data[i]);
	}*/
	//free(vret->data);
	for (i = 0; i < vret->size; i++)
		fprintf(stderr, "%s : %s\n", vret->data[i][0], vret->data[i][1]);
	return ret;
}

