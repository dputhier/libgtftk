/*
 * free_gtf_data.c
 *
 *  Created on: Mar 20, 2017
 *      Author: fafa
 */

#include "libgtftk.h"

extern int compare_row_list(const void *, const void *);

/*
 * global variables in libgtftk.c
 */
extern COLUMN **column;
extern int nb_column;

int update_index_table(COLUMN *col) {
	INDEX *idx;
	int i;

	if (col->index != NULL) {
		idx = col->index[0];
		col->index = (INDEX **)realloc(col->index, col->nb_index * sizeof(INDEX *));
		for (i = 0; i < col->nb_index; i++) {
			col->index[i] = idx;
			idx = idx->next;
		}
	}
	return 0;
}


__attribute__ ((visibility ("default")))
int free_gtf_data(GTF_DATA *gtf_data) {
	int i, j, c;
	GTF_ROW *row;
	INDEX *pindex, *pindex0;
	ROW_LIST *row_list;

	row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
	row_list->token = strdup("*");

	if (gtf_data != NULL) {
		for (i = 0; i < gtf_data->size; i++) {
			row = gtf_data->data[i];
			for (j = 0; j < 8; j++) free(row->field[j]);
			free(row->field);

			for (j = 0; j < row->attributes.nb; j++) {
				free(row->attributes.attr[j]->key);
				free(row->attributes.attr[j]->value);
				free(row->attributes.attr[j]);
			}
			free(row->attributes.attr);
			free(row);
		}
		free(gtf_data->data);
		gtf_data->data = NULL;
		//fprintf(stderr, "freed %d\n", gtf_data->size);
		/*for (c = 0; c < nb_column; c++) {
			//fprintf(stderr, "  col = %s %d\n", column[c]->name, column[c]->nb_index);
			if (column[c]->index != NULL)
				pindex = column[c]->index[0];
			else
				pindex = NULL;
			pindex0 = NULL;
			while (pindex != NULL) {
				if (pindex->gtf_data == gtf_data) {
					//fprintf(stderr, "    freeing index %s\n", pindex->key);
					tdelete(row_list, &(pindex->data), compare_row_list);
					//fprintf(stderr, "    OK\n");
					free(pindex->key);
					column[c]->nb_index--;
					if (pindex0 == NULL) {
						pindex0 = pindex->next;
						free(pindex);
						if (pindex == column[c]->index[0])
							column[c]->index[0] = pindex0;
						pindex = pindex0;
						pindex0 = NULL;
					}
					else {
						pindex0->next = pindex->next;
						free(pindex);
						if (pindex == column[c]->index[0])
							column[c]->index[0] = pindex0->next;
						pindex = pindex0->next;
					}
				}
				else {
					pindex0 = pindex;
					pindex = pindex->next;
				}
			}
			update_index_table(column[c]);
		}*/
		free(gtf_data);
		gtf_data = NULL;
	}
	free(row_list->token);
	free(row_list);
	return 0;
}
