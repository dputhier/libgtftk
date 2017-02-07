/*
 * load_gtf.c
 *
 *  Created on: Jan 12, 2017
 *      Author: fafa
 *
 * This source file contains functions to load a GTF file into memory, to get
 * a GTF_DATA structure and to index the whole thing with a given column name
 * or attribute name. The GTF_DATA structure is the object that must be used
 * as input with most of the functions of the library.
 */

#include "libgtftk.h"
#include <search.h>
/*
 * external functions declaration
 */
extern GTF_READER *get_gtf_reader(char *query);
extern void make_columns();
extern char *get_next_gtf_line(GTF_READER *gr, char *buffer);
extern int split_ip(char ***tab, char *s, char *delim);

/*
 * global variables declaration
 */
extern COLUMN **column;
extern int nb_column;

/*
 * Reserve memory for a new attribute index and add it at the end of the
 * attributes column. Used by index_gtf in this source file.
 *
 * Parameters:
 * 		key:	an attribute name
 *
 * Returns:		the rank of the index in the index table of the column
 */
int add_index(char *key) {
	column[8]->nb_index++;
	column[8]->index = realloc(column[8]->index, column[8]->nb_index * sizeof(INDEX *));
	column[8]->index[column[8]->nb_index - 1] = (INDEX *)calloc(1, sizeof(INDEX));
	column[8]->index[column[8]->nb_index - 1]->key = strdup(key);
	return column[8]->nb_index - 1;
}

int compatt(const void *s1, const void *s2) {
	char *cs1 = *(char **)s1;
	char *cs2 = *(char **)s2;
	return strcmp(cs1, cs2);
}

/*
 * This function read GTF data from a GTF file (gzipped or not) or standard
 * input and makes a GTF_DATA object that contains all the GTF file data.
 * This GTF_DATA is intended to be eventually indexed and used as input of one
 * of the library functions. This function is visible from outside the library.
 *
 * Parameters:
 * 		input:	a GTF file name (gzipped or not) or "-" for standard input
 *
 * Returns:		a GTF_DATA structure filled with GTF data
 */
__attribute__ ((visibility ("default")))
GTF_DATA *load_GTF(char *input) {
	// a buffer to store rows of GTF file as they are read
	char *buffer = (char *)calloc(10000, sizeof(char));

	// a string table to split GTF rows into fields
	char **token;

	int i, nb_row;

	// creates a GTF_READER to read from input
	GTF_READER *gr = get_gtf_reader(input);
	if (gr == NULL) return NULL;

	// reserve memory for the GTF_DATA structure to return
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));

	GTF_ROW *row;

	// make the GTF column model
	make_columns();

	// a counter for the rows
	nb_row = 0;

	// loop on the GTF file rows
	while (get_next_gtf_line(gr, buffer) != NULL) {
		if (*buffer != '#') {
			*(buffer + strlen(buffer) - 1) = 0;

			// reserve memory for a new row in the row table
			ret->data = (GTF_ROW **)realloc(ret->data, (nb_row + 1) * sizeof(GTF_ROW *));

			/* reserve memory for the new row; this row variable is used to
			 * improve readability of the code
			 */
			row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));

			// put the new row at the end of the row table
			ret->data[nb_row] = row;

			// split the row and check the number of fields (should be 9)
			if (split_ip(&token, buffer, "\t") != nb_column) {
				if (!strcmp(gr->filename, "-"))
					fprintf(stderr, "ERROR : standard input is not a valid GTF stream\n");
				else
					fprintf(stderr, "ERROR : GTF file %s is not valid\n", gr->filename);
				exit(0);
			}

			// reserve memory for the 9 fields
			row->data = (void *)calloc(nb_column, sizeof(void *));

			// store the 9 fields after conversion in their appropriate type
			for (i = 0; i < nb_column; i++)
				row->data[i] = column[i]->convert(token[i], column[i]->default_value);

			// keep the rank number of the row
			row->rank = nb_row;

			// one more row has been read
			nb_row++;

			// free memory used to split the row
			free(token);
		}
	}

	// free the buffer used to read GTF file
	free(buffer);

	// store the final number of rows
	ret->size = nb_row;

	return ret;
}

STRING_LIST *get_attribute_list(GTF_DATA *gtf_data) {
	int i, j;
	STRING_LIST *sl = (STRING_LIST *)calloc(1, sizeof(STRING_LIST));
	char **pkey;
	GTF_ROW *row;

	for (j = 0; j < gtf_data->size; j++) {
		row = gtf_data->data[j];
		for (i = 0; i < ((ATTRIBUTES *)row->data[8])->nb; i++) {
			pkey = &((ATTRIBUTES *)row->data[8])->attr[i]->key;
			if (!lfind(pkey, sl->list, (size_t *)(&sl->nb), sizeof(char *), compatt)) {
				sl->list = (char **)realloc(sl->list, (sl->nb + 1) * sizeof(char *));
				lsearch(pkey, sl->list, (size_t *)(&sl->nb), sizeof(char *), compatt);
			}
		}
	}
	return sl;
}
void action_destroy(void *nodep) {

}
/*
 * Indexes a GTF_DATA with a column name or an attribute name. The return value
 * is the rank of the column for a column name index or the rank of the
 * attribute index in the index table of the column + 8. So, a return value
 * greater than 7 means that data has been indexed on an attribute name.
 *
 * Parameters:
 * 		gtf_data:	the GTF_DATA to index
 * 		key:		the column name or attribute name
 *
 * Returns:			the column rank (0-7) or the attribute index rank + 8
 */
int index_gtf(GTF_DATA *gtf_data, char *key) {
	int i, j, k, found;

	found = 0;

	// look for key in the 8 first column names and index the column if found
	for (i = 0; i < (nb_column - 1); i++)
		if (!strcmp(column[i]->name, key)) {
			/*
			 * if an index on this column exists, we need to destroy it
			 */
			if (column[i]->index[0]->data != NULL) {
				tdestroy(column[i]->index[0]->data, action_destroy);
				column[i]->index[0]->data = NULL;
			}
			for (k = 0; k < gtf_data->size; k++)
				column[i]->index_row(k, column[i]->convert_to_string(gtf_data->data[k]->data[i], column[i]->default_value), column[i]->index[0]);
			found = 1;
			break;
		}
	if (!found) {
		/* key is not a column name so look for key in the attribute index
		 * table
		 */
		for (i = 0; i < column[8]->nb_index; i++)
			if (!strcmp(column[8]->index[i]->key, key)) {
				found = 1;

				/*
				 * if an index on this attribute exists, we need to destroy it
				 */
				if (column[8]->index[i]->data != NULL) {
					tdestroy(column[8]->index[i]->data, action_destroy);
					column[8]->index[i]->data = NULL;
				}
				break;
			}

		/* if there is no index for key, we need to make one and add it in the
		 * table
		 */
		if (!found) i = add_index(key);

		/*
		 * Now we have an index for key and i is his rank. We can parse the
		 * GTF_DATA, look for key attribute in each row, and if found, add the
		 * row in the attribute key index
		 */
		for (k = 0; k < gtf_data->size; k++)
			for (j = 0; j < ((ATTRIBUTES *)gtf_data->data[k]->data[8])->nb; j++)
				if (!strcmp(key, ((ATTRIBUTES *)gtf_data->data[k]->data[8])->attr[j]->key)) {
					column[8]->index_row(k, ((ATTRIBUTES *)gtf_data->data[k]->data[8])->attr[j]->value, column[8]->index[i]);
					break;
				}

		// add 8 to rank for an attribute name index
		i += 8;
	}
	return i;
}
