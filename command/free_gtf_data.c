/*
 * free_gtf_data.c
 *
 *  Created on: Mar 20, 2017
 *      Author: fafa
 */

#include "libgtftk.h"

__attribute__ ((visibility ("default")))
int free_gtf_data(GTF_DATA *gtf_data) {
	int i, j;
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
		gtf_data = NULL;
	}
	return 0;
}
