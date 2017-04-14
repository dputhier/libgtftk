/*
 * column.c
 *
 *  Created on: Jan 9, 2017
 *      Author: fafa
 *
 *  This source file contains all that is related to column management :
 *  	- all COLUMN structure functions, for each kind of column (string, int,
 *  	  float, char, attributes): convert, convert_to_string, index_row and
 *  	  print
 *		- print_row to print a whole row with a given delimiter
 *		- index_row to add a row in an index corresponding to the given value
 *		- make_columns/make_column to build the GTF column model
 */

#include "libgtftk.h"

/*
 * external functions declaration
 */
extern int compare_row_list(const void *p1, const void *p2);

/*
 * global variables declaration
 */
extern int nb_column;
extern COLUMN **column;

/*
 * These functions prints in output a column value (a feature, a start, an
 * attributes list ...) with delimiter. Default values are printed as ".".
 * These functions are the "print" function of the COLUMN structure.
 *
 * Parameters:
 * 		token:	the typed value to print
 * 		output:	where to print (a file or stdout)
 * 		col:	the corresponding column (NULL if attributes column)
 * 		delim:	the delimiter character
 */

// print for the 8 first columns
void print_string(void *token, FILE *output, void *col, char delim) {
	if (((COLUMN *)col)->default_value != NULL)
		if (!strcmp((char *)(token), (char *)((COLUMN *)col)->default_value))
			delim != 0 ? fprintf(output, ".%c", delim) : fprintf(output, ".");
		else
			delim != 0 ? fprintf(output, "%s%c", (char *)(token), delim) : fprintf(output, "%s", (char *)(token));
	else
		delim != 0 ? fprintf(output, "%s%c", (char *)(token), delim) : fprintf(output, "%s", (char *)(token));
}

// print for attributes column
void print_attributes(void *token, FILE *output, void *col, char delim) {
	int k;
	GTF_ROW *row = (GTF_ROW *)token;
	if (row->nb_attributes != -1) {
		fprintf(output, "%s \"%s\";", row->key[0], row->value[0]);
		for (k = 1; k < row->nb_attributes; k++)
			fprintf(output, " %s \"%s\";", row->key[k], row->value[k]);
	}
}

void make_index(INDEX_ID **p_index_id, char *key) {
	column[(*p_index_id)->column]->nb_index++;
	column[(*p_index_id)->column]->index = (INDEX *)realloc(column[(*p_index_id)->column]->index, column[(*p_index_id)->column]->nb_index * sizeof(INDEX));
	INDEX *index = column[(*p_index_id)->column]->index + column[(*p_index_id)->column]->nb_index - 1;
	index->data = NULL;
	index->gtf_data = NULL;
	index->key = strdup(key);
	(*p_index_id)->index_rank = column[(*p_index_id)->column]->nb_index - 1;
}

/*
 * This function adds a row (row_nb) into an index.
 *
 * Parameters:
 * 		row_nb:		the row number to add
 * 		value:		the value with which to row is associated
 * 		index:		the index of ROW_LIST elements
 */
void index_row(int row_nb, char *value, INDEX *index) {
	ROW_LIST *test_row_list, *row_list, *find_row_list;

	if (index != NULL) {
		// build a ROW_LIST to check if value is already indexed

		test_row_list = calloc(1, sizeof(ROW_LIST));
		test_row_list->token = value;
		//fprintf(stderr, "ho\n");
		find_row_list = tfind(test_row_list, &(index->data), compare_row_list);
		if (find_row_list == NULL) {
			/*
			 * value is not in the index so we reserve a ROW_LIST, initialize it
			 * with one row (row_nb) and put it in the index (tsearch)
			 */

			//fprintf(stderr, "creating index element %s\n", value);
			row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
			row_list->token = value;
			row_list->nb_row = 1;
			row_list->row = (int *)calloc(1, sizeof(int));
			row_list->row[row_list->nb_row - 1] = row_nb;
			tsearch(row_list, &(index->data), compare_row_list);
		}
		else {
			/*
			 * value is already in the index so we just have to add row_nb into
			 * the ROW_LIST element found
			 */
			//fprintf(stderr, "adding row %d in element %s\n", row_nb, value);
			row_list = *((ROW_LIST **)find_row_list);
			row_list->nb_row++;
			row_list->row = (int *)realloc(row_list->row, row_list->nb_row * sizeof(int));
			row_list->row[row_list->nb_row - 1] = row_nb;
		}
		free(test_row_list);
	}
}

/*
 * Prints a GTF row in output with the given character delimiter
 *
 * Parameters:
 * 		output:		where to print
 * 		r:			the row to print
 * 		delim:		the delimiter character
 */
void print_row(FILE *output, GTF_ROW *r, char delim) {
	int i;

	for (i = 0; i < 8; i++) print_string(r->field[i], output, column[i], delim);
	print_attributes(r, output, NULL, 0);
	fprintf(output, "\n");
}

/*
 * Reserve memory for a COLUMN structure and fill it whith parameters.
 *
 * Parameters:
 * 		type:	the column type ('S' for string and 'A' for attributes)
 * 		i:		the number of the column (0 to 8)
 * 		dv:		the default value
 * 		name:	the column name
 *
 * Returns:		the initialized COLUMN structure
 */
COLUMN *make_column(char type, int i, void *dv, char *name) {
	COLUMN *column = (COLUMN *)calloc(1, sizeof(COLUMN));
	column->num = i;
	column->name = strdup(name);
	column->index = NULL;
	column->nb_index = 0;
	if (dv != NULL) column->default_value = strdup((char *)dv);

	return column;
}

/*
 * Creates and initialize a column model suitable for GTF data.
 * This function creates an index for each column except the last one
 * (attributes). The indexing of a column depends on the kind of operation to
 * perform on GTF data. For example, "select_by_key gene_biotype:protein_coding"
 * needs to index only "gene_biotype" attribute. This function is called by
 * loadGTF function, the first function a client should call to use the library.
 */
void make_columns() {
	int i;

	nb_column = 9;
	column = (COLUMN **)calloc(nb_column, sizeof(COLUMN *));
	column[0] = make_column('S', 0, ".", "seqid");
	column[1] = make_column('S', 1, ".", "source");
	column[2] = make_column('S', 2, ".", "feature");
	column[3] = make_column('S', 3, ".", "start");
	column[4] = make_column('S', 4, ".", "end");
	column[5] = make_column('S', 5, ".", "score");
	column[6] = make_column('S', 6, ".", "strand");
	column[7] = make_column('S', 7, ".", "phase");
	column[8] = make_column('A', 8, ".", "attributes");

	for (i = 0; i < (nb_column - 1); i++) {
		//column[i]->index = (INDEX *)calloc(1, sizeof(INDEX));
		//column[i]->index->key = column[i]->name;
		//column[i]->nb_index = 1;
	}
}
