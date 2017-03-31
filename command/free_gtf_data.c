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

	//fprintf(stderr, "%d\n", gtf_data->size);
	for (i = 0; i < gtf_data->size; i++) {
		row = &gtf_data->data[i];
		for (j = 0; j < 8; j++)	free(row->field[j]);
		free(row->field);

		for (j = 0; j < row->nb_attributes; j++) {
			free(row->key[j]);
			free(row->value[j]);
		}
		free(row->key);
		free(row->value);
	}
	free(gtf_data->data);
	free(gtf_data);
	return 0;
}
