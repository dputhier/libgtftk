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
extern int is_in_attrs(GTF_ROW *row, char *at);

/*
 * global variables declaration
 */
extern COLUMN **column;
extern int nb_column;

/*
 * Look for a column in the column names.
 *
 * Parameters:
 * 		col:	the column name to look for
 *
 * Returns:		the rank of the found column (or -1 if not found)
 */
int is_in_columns(char *col) {
	int ret = -1, i;
	for (i = 0; i < nb_column; i++)
		if (!strcmp(column[i]->name, col)) {
			ret = i;
			break;
		}
	return ret;
}

/*
 * This function extracts data from a GTF_DATA structure and pack it into
 * a RAW_DATA structure, that is a matrix of strings.
 *
 * Parameters:
 * 		gtf_data:	the GTF source data
 * 		key:		the list of column names and attribute names that are to be
 * 					extracted (separated by "," character)
 *
 * Returns:			a matrix of strings, packed into a RAW_DATA structure
 */
__attribute__ ((visibility ("default")))
RAW_DATA *extract_data(GTF_DATA *gtf_data, char *key, int base, int uniq) {
	/*
	 * Declaration and reservation of the structure to be returned
	 */
	RAW_DATA *ret = (RAW_DATA *)calloc(1, sizeof(RAW_DATA));

	/*
	 * Some convenient local variables
	 */
	int i, k, n, j;

	/*
	 * The list of all attributes in the given GTF_DATA
	 */
	STRING_LIST *attributes;

	/*
	 * if key is "all", we extract ALL the data in the GTF_DATA
	 */
	if (!strcmp(key, "all")) {
		/*
		 * get the list of all attributes in the GTF_DATA
		 */
		attributes = get_attribute_list(gtf_data);

		/*
		 * reserve the memory for the list of column names, in the results
		 * (ROW_DATA)
		 */
		ret->column_name = (char **)malloc((8 + attributes->nb) * sizeof(char *));

		/*
		 * Fill the list of result columns with the column names and the
		 * attribute names of the GTF_DATA
		 */
		for (i = 0; i < 8; i++) ret->column_name[ret->nb_columns++] = column[i]->name;
		for (i = 0; i < attributes->nb; i++) ret->column_name[ret->nb_columns++] = strdup(attributes->list[i]);
		free(attributes->list);
	}
	else {
		/*
		 * if key contains a list of column and attribute names, split this
		 * list right in the RAW_DATA structure
		 */
		//fprintf(stderr, "coucou : %s %d %lu\n", key, ret->nb_columns, &(ret->column_name));
		ret->nb_columns = split_ip(&(ret->column_name), key, ",");
		//fprintf(stderr, "coucou : %d\n", ret->nb_columns);
	}
	/*
	 * reserve the memory for the string matrix and setup the number of rows
	 * in the matrix (always the same as in the GTF_DATA)
	 */
	ret->data = (char ***)calloc(gtf_data->size, sizeof(char **));
	ret->nb_rows = gtf_data->size;

	/*
	 * The loop on the rows of GTF_DATA
	 */
	for (k = 0; k < gtf_data->size; k++) {
		/*
		 * reserve the memory for the corresponding result row
		 */
		ret->data[k] = (char **)calloc(ret->nb_columns, sizeof(char *));

		/*
		 * The loop on the columns and attributes to extract
		 */
		for (i = 0; i < ret->nb_columns; i++) {

			/*
			 * column extraction
			 */
			if ((n = is_in_columns(ret->column_name[i])) != -1) {
				/*
				 * do a copy of the field, as it may be altered
				 */
				ret->data[k][i] = strdup(gtf_data->data[k].field[n]);

				/*
				 * if we want data 0-based, we must remove 1 from the start
				 * value
				 */
				if (!strcmp(ret->column_name[i], "start") && !base) {
					j = atoi(ret->data[k][i]) - 1;
					sprintf(ret->data[k][i], "%d", j);
				}
			}

			/*
			 * attribute extraction
			 */
			else if ((n = is_in_attrs(&gtf_data->data[k], ret->column_name[i])) != -1)
				ret->data[k][i] = strdup(gtf_data->data[k].value[n]);

			/*
			 * not a column and not an attribute ! (some attributes are not
			 * present in all the rows of a GTF_DATA)
			 */
			else
				ret->data[k][i] = strdup(".");
		}
	}

	/*
	 * Returns the extracted data
	 */
	return ret;
}
