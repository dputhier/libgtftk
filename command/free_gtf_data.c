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
	GTF_ROW *gtf_row;
	ATTRIBUTES *attributes;
	ATTRIBUTE *attribute;

	fprintf(stderr, "%d\n", gtf_data->size);
	for (i = 0; i < gtf_data->size; i++) {
		gtf_row = gtf_data->data[i];
		//for (j = 0; j < 8; j++)	{
		free((char *)(gtf_row->data[0]));
		free((char *)(gtf_row->data[1]));
		free((char *)(gtf_row->data[2]));
		free((int *)(gtf_row->data[3]));
		free((int *)(gtf_row->data[4]));
		free((float *)(gtf_row->data[5]));
		free((char *)(gtf_row->data[6]));
		free((int *)(gtf_row->data[7]));
		//}
		attributes = gtf_row->data[8];
		/*for (j = 0; j < attributes->nb; j++) {
			attribute = attributes->attr[j];
			//free(attribute->key);
			//free(attribute->value);
			free(attribute);
		}*/
		free(attributes);
		free(gtf_row);
	}
	free(gtf_data->data);
	free(gtf_data);
	return 0;
}
