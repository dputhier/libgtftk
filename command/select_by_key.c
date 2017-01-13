/*
 * select_by_key.c
 *
 *  Created on: 13 janv. 2017
 *      Author: fafa
 */

#include <stdlib.h>
#include <libgtftk.h>

extern int nb_column;
extern COLUMN **column;
extern int split_ip(char ***tab, char *s, char *delim);
extern int compare_row_list(const void *p1, const void *p2);
extern int comprow(const void *m1, const void *m2);
extern int add_row_list(ROW_LIST *src, ROW_LIST *dst);

__attribute__ ((visibility ("default")))
GTF_DATA *select_by_key(GTF_DATA *gtf_data, char *key, char *value, int not) {
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));
	int i, j, k, p, found = 0, n = 0;
	ROW_LIST *test_row_list = calloc(1, sizeof(ROW_LIST)), *row_list = NULL, **find_row_list = NULL;
	char **values;
	int nb_value = split_ip(&values, value, ",");

	row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
	for (i = 0; i < nb_column; i++)
		if (!strcmp(column[i]->name, key)) {
			for (p = 0; p < nb_value; p++) {
				test_row_list->token = values[p];
				find_row_list = (ROW_LIST **)tfind(test_row_list, &(column[i]->index[0]->data), compare_row_list);
				if (find_row_list != NULL) add_row_list(*find_row_list, row_list);
			}
			qsort(row_list->row, row_list->nb_row, sizeof(int), comprow);
			if (!not) {
				ret->size = row_list->nb_row;
				ret->data = (GTF_ROW **)calloc(ret->size, sizeof(GTF_ROW *));
				for (j = 0; j < ret->size; j++) ret->data[j] = gtf_data->data[row_list->row[j]];
			}
			else {
				j = 0;
				ret->size = gtf_data->size - row_list->nb_row;
				ret->data = (GTF_ROW **)calloc(ret->size, sizeof(GTF_ROW *));
				for (k = 0; k < gtf_data->size; k++)
					if (k < row_list->row[j]) {
						ret->data[n] = gtf_data->data[k];
						n++;
					}
					else if (k == row_list->row[j])
						j++;
				if (n != ret->size) {
					for (k = row_list->row[row_list->nb_row - 1] + 1; k < gtf_data->size; k++) {
						ret->data[n] = gtf_data->data[k];
						n++;
					}
				}
			}
			found = 1;
			break;
		}
	if (!found)
		if ((value != NULL) && (column[8]->index != NULL)) {
			for (p = 0; p < nb_value; p++) {
				test_row_list->token = values[p];
				for (i = 0; i < column[8]->nb_index; i++)
					if (!strcmp(column[8]->index[i]->key, key)) {
						find_row_list = tfind(test_row_list, &(column[8]->index[i]->data), compare_row_list);
						break;
					}
				if (find_row_list != NULL) add_row_list(*find_row_list, row_list);
			}
			if (row_list != NULL) {
				qsort(row_list->row, row_list->nb_row, sizeof(int), comprow);
				if (!not) {
					ret->size = row_list->nb_row;
					ret->data = (GTF_ROW **)calloc(ret->size, sizeof(GTF_ROW *));
					for (j = 0; j < row_list->nb_row; j++) ret->data[j] = gtf_data->data[row_list->row[j]];
				}
				else {
					j = 0;
					ret->size = gtf_data->size - row_list->nb_row;
					ret->data = (GTF_ROW **)calloc(ret->size, sizeof(GTF_ROW *));
					for (k = 0; i < gtf_data->size; k++)
						if (k < row_list->row[j]) {
							ret->data[n] = gtf_data->data[k];
							n++;
						}
						else if (k == row_list->row[j])
							j++;
					if (n != ret->size) {
						for (k = row_list->row[row_list->nb_row - 1] + 1; k < ret->size; k++) {
							ret->data[n] = gtf_data->data[k];
							n++;
						}
					}
				}
			}
		}
	return ret;
}

