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
 * In a GTF_DATA, the rows are stored in two manners: a linked list and a table.
 * The linked list is usefull for adding rows (ie convert_to_ensembl function)
 * while the table is quicker to address a particular row.
 * After some operations, the rows can be stored only in the linked list (a
 * pointer on the first row is stored in gtf_data->data[0]). This function
 * rebuild the corresponding table in gtf_data->data from the linked list.
 *
 * Parameters:
 * 		gtf_data:	the input GTF data
 *
 * Returns:			0 if success
 */
int update_row_table(GTF_DATA *gtf_data) {
	GTF_ROW *row = gtf_data->data[0];
	int i;

	gtf_data->data = (GTF_ROW **)realloc(gtf_data->data, gtf_data->size * sizeof(GTF_ROW *));
	for (i = 0; i < gtf_data->size; i++) {
		gtf_data->data[i] = row;
		row = row->next;
	}
	return 0;
}

/*
 * The symetric function of update_row_table. This function rebuild the linked
 * list from the table in gtf_data->data.
 *
 * Parameters:
 * 		gtf_data:	the input GTF data
 *
 * Returns:			0 if success
 */
int update_linked_list(GTF_DATA *gtf_data) {
	int i, j;
	GTF_ROW *row;

	for (i = 0; i < gtf_data->size - 1; i++) {
		if ((gtf_data->data[i]->next != NULL) && (gtf_data->data[i]->next != gtf_data->data[i + 1])) {
			row = gtf_data->data[i]->next;
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
		gtf_data->data[i]->next = gtf_data->data[i + 1];
	}
	return 0;
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

	// reserve memory for one row in the table.
	ret->data = (GTF_ROW **)calloc(1, sizeof(GTF_ROW *));

	GTF_ROW *row, *previous_row;

	// make the GTF column model
	make_columns();

	// a counter for the rows
	nb_row = 0;

	// loop on the GTF file rows
	while (get_next_gtf_line(gr, buffer) != NULL) {
		if (*buffer != '#') {
			*(buffer + strlen(buffer) - 1) = 0;

			// reserve memory for a new row in the row list
			row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
			if (nb_row == 0) ret->data[0] = row;
			//ret->data = (GTF_ROW *)realloc(ret->data, (nb_row + 1) * sizeof(GTF_ROW));

			/*
			 * reserve memory for the new row; this row variable is used to
			 * improve readability of the code
			 */
			//row = ret->data + nb_row;

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
			row->attributes.nb = split_ip(&attr, token[8], ";");

			// reserve memory for the attributes
			row->attributes.attr = (ATTRIBUTE **)calloc(row->attributes.nb, sizeof(ATTRIBUTE *));

			// store the attributes
			for (i = 0; i < row->attributes.nb; i++) {
				row->attributes.attr[i] = (ATTRIBUTE *)calloc(1, sizeof(ATTRIBUTE));
				split_key_value(attr[i], &(row->attributes.attr[i]->key), &(row->attributes.attr[i]->value));
			}

			// keep the rank number of the row
			row->rank = nb_row;

			// insert the new row in the list
			//if (nb_row > 0) (ret->data + nb_row - 1)->next = row;
			if (nb_row > 0) previous_row->next = row;
			previous_row = row;

			// one more row has been read
			nb_row++;

			// free memory used to split the row and the attributes
			free(token);
			free(attr);
		}
	}

	// store the final number of rows
	ret->size = nb_row;

	// update the row tlist
	update_row_table(ret);

	// free the buffer used to read GTF file
	free(buffer);

	return ret;
}

STRING_LIST *get_all_attributes(GTF_DATA *gtf_data) {
	int i, j;
	STRING_LIST *sl = (STRING_LIST *)calloc(1, sizeof(STRING_LIST));
	char **pkey;
	GTF_ROW *row;

	for (j = 0; j < gtf_data->size; j++) {
		row = gtf_data->data[j];
		for (i = 0; i < row->attributes.nb; i++) {
			pkey = &(row->attributes.attr[i]->key);
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
		make_index(&index_id, key);
		if (index_id->column != 8) {
			for (k = 0; k < gtf_data->size; k++)
				index_row(k, gtf_data->data[k]->field[index_id->column], column[index_id->column]->index + index_id->index_rank);
			column[index_id->column]->index[index_id->index_rank].gtf_data = gtf_data;
		}
		else {
			for (k = 0; k < gtf_data->size; k++)
				for (j = 0; j < gtf_data->data[k]->attributes.nb; j++)
					if (!strcmp(key, gtf_data->data[k]->attributes.attr[j]->key)) {
						index_row(k, gtf_data->data[k]->attributes.attr[j]->value, column[index_id->column]->index + index_id->index_rank);
						break;
					}
			column[index_id->column]->index[index_id->index_rank].gtf_data = gtf_data;
		}
	}

	return index_id;
}
