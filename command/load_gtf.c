/*
 * load_gtf.c
 *
 *  Created on: Jan 12, 2017
 *      Author: fafa
 */

#include <libgtftk.h>

extern GTF_READER *get_gtf_reader(char *query);
extern void make_columns();
extern char *get_next_gtf_line(GTF_READER *gr, char *buffer);
extern int split_ip(char ***tab, char *s, char *delim);

extern COLUMN **column;
extern int nb_column;

__attribute__ ((visibility ("default")))
GTF_DATA *load_GTF(char *input) {
	char *buffer = (char *)calloc(10000, sizeof(char));
	char **token;
	ATTRIBUTES *a;
	int i, l, found, k, nb_row;

	GTF_READER *gr = get_gtf_reader(input);
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));
	GTF_ROW *row;

	make_columns();
	nb_row = 0;
	while (get_next_gtf_line(gr, buffer) != NULL) {
		if (*buffer != '#') {
			*(buffer + strlen(buffer) - 1) = 0;
			ret->data = (GTF_ROW **)realloc(ret->data, (nb_row + 1) * sizeof(GTF_ROW *));
			row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
			ret->data[nb_row] = row;

			if (split_ip(&token, buffer, "\t") != nb_column) {
				if (!strcmp(gr->filename, "-"))
					fprintf(stderr, "ERROR : standard input is not a valid GTF stream\n");
				else
					fprintf(stderr, "ERROR : GTF file %s is not valid\n", gr->filename);
				exit(0);
			}
			row->data = (void *)calloc(nb_column, sizeof(void *));
			for (i = 0; i < (nb_column - 1); i++) {
				row->data[i] = column[i]->convert(token[i], column[i]->default_value, row);
				column[i]->make_index(nb_row, row, column[i]);
			}
			row->data[8] = column[8]->convert(token[8], column[8]->default_value, row);
			a = (ATTRIBUTES *)(row->data[8]);
			for (l = 0; l < a->nb; l++) {
				found = 0;
				for (k = 0; k < column[8]->nb_index; k++)
					if (!strcmp(a->attr[l]->key, column[8]->index[k]->key)) {
						found = 1;
						break;
					}
				if (!found) {
					column[8]->nb_index++;
					column[8]->index = realloc(column[8]->index, column[8]->nb_index * sizeof(INDEX *));
					column[8]->index[column[8]->nb_index - 1] = (INDEX *)calloc(1, sizeof(INDEX));
					column[8]->index[column[8]->nb_index - 1]->key = strdup(a->attr[l]->key);
					column[8]->make_index(nb_row, row, column[8]);
				}
			}

			nb_row++;
			free(token);
		}
	}
	free(buffer);
	ret->size = nb_row;
	return ret;
}
