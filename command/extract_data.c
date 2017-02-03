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
extern int get_attribute_list(FILE *output);

/*
 * global variables declaration
 */
extern COLUMN **column;

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

RAW_DATA *extract_data(GTF_DATA *gtf_data, char *key) {
	RAW_DATA *ret = (RAW_DATA *)calloc(1, sizeof(RAW_DATA));
/*	int k, i, n, nb;
	char **keys;
	char *line = (char *)calloc(10000, sizeof(char));
	int nb_keys;
	char *buffer;
	size_t size = 0;
	FILE *out = NULL, *in;

	if (!strcmp(key, "all")) {
		nb_keys = 8;
		keys = (char **)malloc(nb_keys * sizeof(char *));
		n = 0;
		for (i = 0; i < nb_column; i++)
			if (column[i]->type != 'A')
				keys[n++] = column[i]->name;
		size = 0;
		out = open_memstream(&buffer, &size);
		n = get_attribute_list(out);
		fflush(out);
		fclose(out);
		keys = (char **)realloc(keys, (nb_keys + n) * sizeof(char *));
		in = fmemopen(buffer, size, "r");
		while (fgets(line, 9999, in) != NULL) {
			*(strchr(line, (int)':')) = 0;
			keys[nb_keys++] = strdup(line);
		}
		free(buffer);
	}
	else
		nb_keys = split_ip(&keys, key, ",");

	fprintf(output, "# %s", keys[0]);
	for (i = 1; i < nb_keys; i++) fprintf(output, "\t%s", keys[i]);
	fprintf(output, "\n");
	nb = 1;
	for (k = 0; k < nb_row; k++) {
		ATTRIBUTE **attr = ((ATTRIBUTES *)(data[k]->data[8]))->attr;
		for (i = 0; i < nb_keys - 1; i++) {
			if ((n = is_in_columns(column, nb_column, keys[i])) != -1)
				column[n]->print(data[k]->data[n], output, column[n], '\t');
			else if ((n = is_in_attrs(data[k]->data[8], keys[i])) != -1)
				fprintf(output, "%s\t", attr[n]->value);
			else
				fprintf(output, ".\t");
		}
		if ((n = is_in_columns(column, nb_column, keys[i])) != -1)
			column[n]->print(data[k]->data[n], output, column[n], 0);
		else if ((n = is_in_attrs(data[k]->data[8], keys[i])) != -1)
			fprintf(output, "%s", attr[n]->value);
		else
			fprintf(output, ".");
		fprintf(output, "\n");
		nb++;
	}
	free(keys);*/
	return ret;
}
