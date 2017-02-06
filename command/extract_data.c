/*
 * extract_data.c
 *
 *  Created on: Mar 31, 2016
 *      Author: fafa
 *
 * The extract_data function returns a RAW_DATA structure that contains a
 * selection of information from a GTF_DATA structure (some columns or some
 * attributes).
 */

#include  "libgtftk.h"

/*
 * external functions declaration
 */
extern int split_ip(char ***tab, char *s, char *delim);
extern STRING_LIST *get_attribute_list(GTF_DATA *gtf_data);

/*
 * global variables declaration
 */
extern COLUMN **column;
extern int nb_column;

int is_in_columns(COLUMN **column, int nb_column, char *key) {
	int ret = -1, i;
	for (i = 0; i < nb_column; i++)
		if (!strcmp(column[i]->name, key)) {
			ret = i;
			break;
		}
	return ret;
}

int is_in_attrs(ATTRIBUTES *attrs, char *key) {
	int ret = -1, i;
	for (i = 0; i < attrs->nb; i++)
		if (!strcmp(attrs->attr[i]->key, key)) {
			ret = i;
			break;
		}
	return ret;
}

__attribute__ ((visibility ("default")))
RAW_DATA *extract_data(GTF_DATA *gtf_data, char *key) {
	RAW_DATA *ret = (RAW_DATA *)calloc(1, sizeof(RAW_DATA));
	int i, k, n;
	STRING_LIST *attributes;
	ATTRIBUTES *attrs;

	if (!strcmp(key, "all")) {
		attributes = get_attribute_list(gtf_data);
		ret->column_name = (char **)malloc((8 + attributes->nb) * sizeof(char *));
		ret->nb_columns = 0;
		for (i = 0; i < 8; i++) ret->column_name[ret->nb_columns++] = column[i]->name;
		for (i = 0; i < attributes->nb; i++) ret->column_name[ret->nb_columns++] = strdup(attributes->list[i]);
		free(attributes->list);
	}
	else
		ret->nb_columns = split_ip(&(ret->column_name), key, ",");

	ret->data = (char ***)calloc(gtf_data->size, sizeof(char **));
	ret->nb_rows = gtf_data->size;
	for (k = 0; k < gtf_data->size; k++) {
		attrs = (ATTRIBUTES *)(gtf_data->data[k]->data[8]);
		ret->data[k] = (char **)calloc(ret->nb_columns, sizeof(char *));
		for (i = 0; i < ret->nb_columns - 1; i++) {
			if ((n = is_in_columns(column, nb_column, ret->column_name[i])) != -1)
				ret->data[k][i] = column[n]->convert_to_string(gtf_data->data[k]->data[n], column[n]->default_value);
			else if ((n = is_in_attrs(attrs, ret->column_name[i])) != -1)
				ret->data[k][i] = attrs->attr[n]->value;
			else
				ret->data[k][i] = strdup(".");
		}
		if ((n = is_in_columns(column, nb_column, ret->column_name[i])) != -1)
			ret->data[k][i] = column[n]->convert_to_string(gtf_data->data[k]->data[n], column[n]->default_value);
		else if ((n = is_in_attrs(attrs, ret->column_name[i])) != -1)
			ret->data[k][i] = attrs->attr[n]->value;
		else
			ret->data[k][i] = strdup(".");
	}
	return ret;
}
