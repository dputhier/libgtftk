/*
 * free_gtf_data.c
 *
 *  Created on: Mar 20, 2017
 *      Author: fafa
 */

#include "libgtftk.h"

/*
 * global variables in libgtftk.c
 */
extern COLUMN **column;
extern int nb_column;

__attribute__ ((visibility ("default")))
int free_gtf_data(GTF_DATA *gtf_data) {
	int i, j, c;
	GTF_ROW *row;

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
		free(gtf_data);
		for (c = 0; c < (nb_column - 1); c++) {
			for (i = 0; i < column[c]->nb_index; i++) {
				if (column[c]->index[i].gtf_data == gtf_data) {
					column[c]->index[i].gtf_data = NULL;
				}
			}
		}
		for (c = 0; c < column[8]->nb_index; c++) {
			if (column[8]->index[c].gtf_data == gtf_data) {
				column[8]->index[c].gtf_data = NULL;
			}
		}
		gtf_data = NULL;
	}
	return 0;
}
