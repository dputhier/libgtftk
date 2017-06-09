/*
 * del_attributes.c
 *
 *  Created on: Jun 8, 2017
 *      Author: fafa
 */

#include "libgtftk.h"

/*
 * external functions declaration
 */
extern GTF_DATA *clone(GTF_DATA *gtf_data);

__attribute__ ((visibility ("default")))
GTF_DATA *del_attributes(GTF_DATA *gtf_data, char *features, char *keys) {
	int i, ok, j;

	/*
	 * reserve memory for the GTF_DATA structure to return
	 */
	GTF_DATA *ret = clone(gtf_data);

	GTF_ROW *row;

	for (i = 0; i < ret->size; i++) {
		row = ret->data[i];
		ok = (features == NULL);
		if (!ok) ok = (strstr(features, row->field[2]) != NULL);
		if (ok) {
			/*for (j = 0; j < row->nb_attributes; j++) {
				if (strstr(keys, row->key[j])) {
					free(row->key[j]);
					free(row->value[j]);

				}
			}*/
		}
	}

	return ret;
}
