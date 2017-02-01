/*
 * select_by_genomic_location.c
 *
 *  Created on: Mar 30, 2016
 *      Author: fafa
 */

#include "libgtftk.h"

/*
 * external functions declaration
 */
extern int index_gtf(GTF_DATA *gtf_data, char *key);
extern int compare_row_list(const void *p1, const void *p2);

/*
 * global variables declaration
 */
extern COLUMN **column;
extern int nb_column;

int i, start, end;
ROW_LIST **find_row_list, *row_list, *test_row_list;

COLUMN *get_column_by_name(char *name, COLUMN **column, int nb_column) {
	int i;

	for (i = 0; i < nb_column; i++)
		if (!strcmp(column[i]->name, name))
			return column[i];
	return NULL;
}

__attribute__ ((visibility ("default")))
GTF_DATA *select_by_genomic_location(GTF_DATA *gtf_data, char *chr, int begin_gl, int end_gl) {
	int i;

	/*
	 * reserve memory for the GTF_DATA structure to return
	 */
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));

	/*
	 * indexes the GTF_DATA with seqid column to create an index containing,
	 * for each chromosome, the list of his related rows.
	 */
	i = index_gtf(gtf_data, "seqid");

	/*
	 * Finds all rows related to chromosome "chr"
	 */
	test_row_list = calloc(1, sizeof(ROW_LIST));
	test_row_list->token = chr;
	find_row_list = (ROW_LIST **)tfind(test_row_list, &(column[0]->index[0]->data), compare_row_list);

	/*
	 * Get the number of columns start and end (3 and 4 in GTF)
	 */
	int start_col = get_column_by_name("start", column, nb_column)->num;
	int end_col = get_column_by_name("end", column, nb_column)->num;

	/*
	 * If "chr" was in GTF file
	 */
	if (find_row_list != NULL) {
		/*
		 * Dereferencing found row list
		 */
		row_list = *find_row_list;

		/*
		 * now we fill the resulting GTF_DATA with the found rows and return it
		 */
		ret->data = NULL;
		for (i = 0; i < row_list->nb_row; i++) {
			/*
			 * For each row, get the start and end values
			 */
			start = *((int *)gtf_data->data[row_list->row[i]]->data[start_col]);
			end = *((int *)gtf_data->data[row_list->row[i]]->data[end_col]);

			/*
			 * If start and end values of the row match with the given interval,
			 * add a GTF_ROW in the results
			 */
			if ((begin_gl >= start && begin_gl <= end) || (end_gl >= start && end_gl <= end) || (begin_gl <= start && end_gl >= end)) {
				ret->data = (GTF_ROW **)realloc(ret->data, (ret->size + 1) * sizeof(GTF_ROW *));
				ret->data[ret->size] = gtf_data->data[row_list->row[i]];
				ret->size++;
			}
		}
	}
	return ret;
}
