/*
 * select_by_key.c
 *
 *  Created on: Jan 25, 2017
 *      Author: fafa
 *
 * This is the implementation of select_by_key function.
 * This function selects rows in GTF data that contains (or not) some values
 * for a column or for an attribute.
 */

#include "libgtftk.h"

/*
 * external functions declaration
 */
extern int comprow(const void *m1, const void *m2);
extern int compare_row_list(const void *p1, const void *p2);
extern int add_row_list(ROW_LIST *src, ROW_LIST *dst);
extern int index_gtf(GTF_DATA *gtf_data, char *key);
extern int split_ip(char ***tab, char *s, char *delim);

/*
 * global variables declaration
 */
extern COLUMN **column;

/*
 * select_by_key function selects rows in GTF_DATA that contains some given
 * value(s) for a column or an attribute. For instance:
 *  - to get all "gene" and "transcript" rows:
 *  	 select_by_key(gtf_data, "feature", "gene,transcript", 0)
 * 	- to get all the rows concerning the genes ESR1 and ESR2:
 * 		select_by_key(gtf_data, "gene_name", "ESR1,ESR2", 0)
 * 	- to search for all lincRNA genes:
 * 		select_by_key(gtf_data, "gene_biotype", "lincRNA", 0)
 * 	- to get all non coding genes:
 * 		select_by_key(gtf_data, "gene_biotype", "protein_coding", 1)
 * As output from select_by_key can be used as input, one can chain several
 * function calls. For instance, to get the first exon af all transcripts of
 * gene ESR2:
 * 	select_by_key(select_by_key(select_by_key(gtf_data,
 * 		"gene_name", "ESR2", 0),
 * 		"feature", "exon", 0),
 * 		"exon_number", "1", 0)
 *
 * Parameters:
 * 		gtf_data:	a GTF_DATA structure obtained by a call to loadGTF function
 * 		key:		a column name or an attribute name
 * 		value:		a set of values (separated by "," characters) for key to
 * 					look for in the GTF_DATA
 * 		not:		1 to get the rows that don't contains the requested values
 *
 * Returns:		a GTF_DATA structure that contains the result of the query
 */
__attribute__ ((visibility ("default")))
GTF_DATA *select_by_key(GTF_DATA *gtf_data, char *key, char *value, int not) {
	int i, j, k, p, n = 0;

	// reserve memory for the GTF_DATA structure to return
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));

	/*
	 * Some ROW_LIST variables needed to perform different tasks:
	 * 	- test_row_list will be used to look for a ROW_LIST element in an index
	 * 	- row_list is used to build the list of rows to keep in result
	 * 	- find_row_list contains the found row list for each value, that will
	 * 	  be merged into row_list
	 */
	ROW_LIST *test_row_list = calloc(1, sizeof(ROW_LIST)), *row_list = NULL, **find_row_list = NULL;

	// splitting the given values with the "," character
	char **values;
	int nb_value = split_ip(&values, value, ",");

	/*
	 * indexing the GTF_DATA with the given key (column name or attribute name)
	 * and get the rank of the concerned index (i)
	 */
	i = index_gtf(gtf_data, key);
	// kept just to remind me how to parse a binary tree with twalk !
	/*for (k = 0; k < column[8]->nb_index; k++) {
		N = 0;
		twalk(column[8]->index[k]->data, action_nb);
		fprintf(stderr, "%s : %d\n", column[8]->index[k]->key, N);
	}*/

	// reserve memory for the final ROW_LIST
	row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));

	/*
	 * this test on i is intended to setup the values of i and k.
	 * i must contain the rank of the concerned column (0-8) and k must
	 * contain the value of the index rank in the column (0 for the 8 first
	 * columns that contains only one index, and another value for the 9th
	 * column)
	 */
	if (i < 8) {
		/*
		 * the given key is a column name so i is the column rank and doesn't
		 * need to be updated. k is setup to 0 because the 8 first columns
		 * contain only one index
		 */
		k = 0;
	}
	else {
		/*
		 * the given key is an attribute name so k is setup to i - 8 to
		 * retrieve the value of the concerned attribute index and i is setup
		 * to 8 (the rank of the 9th column)
		 */
		k = i - 8;
		i = 8;
	}

	for (p = 0; p < nb_value; p++) {
		/*
		 * for each value, we look for it in the kth index of the ith column
		 * and if we find it, we merge (add_row_list) the returned ROW_LIST
		 * (find_row_list) with the final one (row_list). Duplicates are
		 * suppressed by add_row_list function.
		 */
		test_row_list->token = values[p];
		find_row_list = (ROW_LIST **)tfind(test_row_list, &(column[i]->index[k].data), compare_row_list);
		if (find_row_list != NULL) add_row_list(*find_row_list, row_list);
	}

	/*
	 * as searches were made in the order of the given values, row_list may
	 * contains unsorted rows, so we need to sort them in ascending order for
	 * the output not to be mixed
	 */
	qsort(row_list->row, row_list->nb_row, sizeof(int), comprow);

	/*
	 * now it's time to build the GTF_DATA result
	 */
	if (!not) {
		/*
		 * if the query was positive, the size of the result is the number of
		 * rows in row_list
		 */
		ret->size = row_list->nb_row;

		/*
		 * we reserve memory for the table of rows
		 */
		ret->data = (GTF_ROW *)calloc(ret->size, sizeof(GTF_ROW));

		/*
		 * each row in row_list is a number used to get the real GTF_ROW in the
		 * whole GTF data
		 */
		for (j = 0; j < ret->size; j++) {
			ret->data[j].rank = gtf_data->data[row_list->row[j]].rank;
			ret->data[j].nb_attributes = gtf_data->data[row_list->row[j]].nb_attributes;

			ret->data[j].field = (char **)calloc(8, sizeof(char*));
			for (i = 0; i < 8; i++) ret->data[j].field[i] = strdup(gtf_data->data[row_list->row[j]].field[i]);
			//ret->data[j].field = gtf_data->data[row_list->row[j]].field;

			ret->data[j].value = (char **)calloc(gtf_data->data[row_list->row[j]].nb_attributes, sizeof(char*));
			for (i = 0; i < gtf_data->data[row_list->row[j]].nb_attributes; i++)
				ret->data[j].value[i] = strdup(gtf_data->data[row_list->row[j]].value[i]);
			//ret->data[j].value = gtf_data->data[row_list->row[j]].value;

			ret->data[j].key = (char **)calloc(gtf_data->data[row_list->row[j]].nb_attributes, sizeof(char*));
			for (i = 0; i < gtf_data->data[row_list->row[j]].nb_attributes; i++)
				ret->data[j].key[i] = strdup(gtf_data->data[row_list->row[j]].key[i]);
			//ret->data[j].key = gtf_data->data[row_list->row[j]].key;
		}
	}
	else {
		/*
		 * if the query was negative, we want to keep all the rows that are NOT
		 * in row_list, so the size of the result is the difference between the
		 * total number of rows and the number of rows in row_list
		 */
		ret->size = gtf_data->size - row_list->nb_row;

		/*
		 * we reserve memory for the table of rows
		 */
		ret->data = (GTF_ROW *)calloc(ret->size, sizeof(GTF_ROW));

		/*
		 * an ugly code to get the "complement" rows in gtf_data
		 */
		j = 0;
		for (k = 0; k < gtf_data->size; k++) {
			if (k < row_list->row[j]) {
				ret->data[n].field = gtf_data->data[k].field;
				ret->data[n].key = gtf_data->data[k].key;
				ret->data[n].value = gtf_data->data[k].value;
				ret->data[n].nb_attributes = gtf_data->data[k].nb_attributes;
				ret->data[n].rank = gtf_data->data[k].rank;
				n++;
			}
			else if (k == row_list->row[j])
				j++;
		}
		if (n != ret->size) {
			for (k = row_list->row[row_list->nb_row - 1] + 1; k < gtf_data->size; k++) {
				ret->data[n].field = gtf_data->data[k].field;
				ret->data[n].key = gtf_data->data[k].key;
				ret->data[n].value = gtf_data->data[k].value;
				ret->data[n].nb_attributes = gtf_data->data[k].nb_attributes;
				ret->data[n].rank = gtf_data->data[k].rank;
				n++;
			}
		}
	}

	free(values);
	free(test_row_list);
	free(row_list);

	return ret;
}

