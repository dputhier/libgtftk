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
extern void split_key_value(char *s, char **key, char **value);
extern void index_row(int row_nb, char *value, INDEX *index);
extern INDEX *make_index(INDEX_ID **index_id, char *key);

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
int add_index(char *key, GTF_DATA *gtf_data) {
	column[8]->nb_index++;
	column[8]->index = realloc(column[8]->index, column[8]->nb_index * sizeof(INDEX));
	column[8]->index[column[8]->nb_index - 1].key = strdup(key);
	column[8]->index[column[8]->nb_index - 1].data = NULL;
	column[8]->index[column[8]->nb_index - 1].gtf_data = gtf_data;
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
	char **token, **attr;

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
			ret->data = (GTF_ROW *)realloc(ret->data, (nb_row + 1) * sizeof(GTF_ROW));

			/*
			 * reserve memory for the new row; this row variable is used to
			 * improve readability of the code
			 */
			row = ret->data + nb_row;

			// split the row and check the number of fields (should be 9)
			if (split_ip(&token, buffer, "\t") != nb_column) {
				if (!strcmp(gr->filename, "-"))
					fprintf(stderr, "ERROR : standard input is not a valid GTF stream\n");
				else
					fprintf(stderr, "ERROR : GTF file %s is not valid\n", gr->filename);
				exit(0);
			}

			// reserve memory for the 8 first fields
			row->field = (char **)calloc(8, sizeof(char *));

			// store the 8 first fields
			for (i = 0; i < 8; i++) row->field[i] = strdup(token[i]);

			// split the attributes
			row->nb_attributes = split_ip(&attr, token[8], ";");

			// reserve memory for the attributes
			row->key = (char **)calloc(row->nb_attributes, sizeof(char *));
			row->value = (char **)calloc(row->nb_attributes, sizeof(char *));

			// store the attributes
			for (i = 0; i < row->nb_attributes; i++)
				split_key_value(attr[i], &row->key[i], &row->value[i]);

			// keep the rank number of the row
			row->rank = nb_row;

			// one more row has been read
			nb_row++;

			// free memory used to split the row and the attributes
			free(token);
			free(attr);
		}
	}

	// free the buffer used to read GTF file
	free(buffer);

	// store the final number of rows
	ret->size = nb_row;

	return ret;
}

STRING_LIST *get_all_attributes(GTF_DATA *gtf_data) {
	int i, j;
	STRING_LIST *sl = (STRING_LIST *)calloc(1, sizeof(STRING_LIST));
	char **pkey;
	GTF_ROW *row;

	for (j = 0; j < gtf_data->size; j++) {
		row = &gtf_data->data[j];
		for (i = 0; i < row->nb_attributes; i++) {
			pkey = &row->key[i];
			if (!lfind(pkey, sl->list, (size_t *)(&sl->nb), sizeof(char *), compatt)) {
				sl->list = (char **)realloc(sl->list, (sl->nb + 1) * sizeof(char *));
				lsearch(pkey, sl->list, (size_t *)(&sl->nb), sizeof(char *), compatt);
			}
		}
	}
	return sl;
}

/*
 * an empty function used to destroy the indexes when needed (tdestroy)
 */
void action_destroy(void *nodep) {

}

INDEX_ID *get_index(GTF_DATA *gtf_data, char *key) {
	int c, k;
	INDEX_ID *ret = calloc(1, sizeof(INDEX_ID));

	ret->column = ret->index_rank = -1;

	// look for key in the 8 first column names and index the column if found
	for (c = 0; c < (nb_column - 1); c++)
		if (!strcmp(column[c]->name, key)) {
			/*
			 * key is the name of a column
			 * if an index on this column exists and correspond to the gtf_data
			 * parameter, we return it
			 */
			ret->column = c;
			for (k = 0; k < column[c]->nb_index; k++) {
				if ((column[c]->index[k].data != NULL) && (column[c]->index[k].gtf_data == gtf_data)) {
					ret->index_rank = k;
					//ret = &(column[c]->index[k]);
				}
				break;
			}
			break;
		}

	if (ret->column == -1) {
		/*
		 * key is not a column name so look for key in the attribute index
		 * table of the 8th column
		 */
		ret->column = 8;
		for (c = 0; c < column[8]->nb_index; c++)
			/*
			 * if an index on this parameter exists and correspond to the
			 * gtf_data parameter, we return it
			 */
			if (!strcmp(column[8]->index[c].key, key)) {
				/*
				 * key is the name of a parameter
				 * if an index on this parameter exists and correspond to the
				 * gtf_data parameter, we return it
				 */
				if ((column[8]->index[c].data != NULL) && (column[8]->index[c].gtf_data == gtf_data)) {
					ret->index_rank = c;
					//ret = &(column[8]->index[c]);
				}
				break;
			}
	}
	return ret;
}

INDEX_ID *index_gtf(GTF_DATA *gtf_data, char *key, int *index_rank) {
	int j, k;

	INDEX_ID *index_id = get_index(gtf_data, key);

	if (index_id->index_rank == -1) {
		//fprintf(stderr, "NULL index for %s\n", key);
		make_index(&index_id, key);
		//fprintf(stderr, "coucoucou : %s %d\n", column[index_id->column]->index[index_id->index_rank]->key, column[index_id->column]->nb_index);
		if (index_id->column != 8) {
			for (k = 0; k < gtf_data->size; k++) {
				//fprintf(stderr, "indexing : %d %s ... ", k, gtf_data->data[k].field[i]);
				index_row(k, gtf_data->data[k].field[index_id->column], column[index_id->column]->index + index_id->index_rank);
				//fprintf(stderr, "done.\n");
			}
			column[index_id->column]->index[index_id->index_rank].gtf_data = gtf_data;
		}
		else {
			//i = add_index(key, gtf_data);
			for (k = 0; k < gtf_data->size; k++)
				for (j = 0; j < gtf_data->data[k].nb_attributes; j++)
					if (!strcmp(key, gtf_data->data[k].key[j])) {
						index_row(k, gtf_data->data[k].value[j], column[index_id->column]->index + index_id->index_rank);
						break;
					}
			//i += 8;
			column[index_id->column]->index[index_id->index_rank].gtf_data = gtf_data;
		}
	}

	//fprintf(stderr, "i = %d\n", i);
	return index_id;
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
int index_gtf_old(GTF_DATA *gtf_data, char *key) {
	int i, j, k, found;

	found = 0;

	/*
	 * look for key in the 8 first column names and index the column if found
	 */
	for (i = 0; i < (nb_column - 1); i++)
		if (!strcmp(column[i]->name, key)) {
			/*
			 * if an index on this column doesn't exist, we create it by
			 * indexing each row of the GTF data
			 */
			if (column[i]->index->data == NULL)
				for (k = 0; k < gtf_data->size; k++)
					index_row(k, gtf_data->data[k].field[i], column[i]->index);
			found = 1;
			break;
		}
	if (!found) {
		/*
		 * key is not a column name so look for key in the attribute index
		 * table
		 */
		for (i = 0; i < column[8]->nb_index; i++)
			if (!strcmp(column[8]->index[i].key, key)) {
				found = 1;
				break;
			}

		/*
		 * if there is no index for key, we need to make one and add it in the
		 * table
		 */
		if (!found) {
			i = add_index(key, gtf_data);

			/*
			 * Now we have an index for key and i is his rank. We can parse the
			 * GTF_DATA, look for key attribute in each row, and if found, add the
			 * row in the attribute key index
			 */
			for (k = 0; k < gtf_data->size; k++)
				for (j = 0; j < gtf_data->data[k].nb_attributes; j++)
					if (!strcmp(key, gtf_data->data[k].key[j])) {
						index_row(k, gtf_data->data[k].value[j], column[8]->index + i);
						break;
					}

		}
		/*
		 * add 8 to rank for an attribute name index
		 */
		i += 8;
	}
	return i;
}
