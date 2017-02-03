/*
 * select_by_genomic_location.c
 *
 *  Created on: Mar 30, 2016
 *      Author: fafa
 *
 * This function selects all rows in GTF data with coordinates that intersect
 * a given genomic location.
 */

#include "libgtftk.h"

/*
 * global variables declaration
 */
extern COLUMN **column;

/*
 * external functions declaration
 */
extern int index_gtf(GTF_DATA *gtf_data, char *key);
extern int compare_row_list(const void *p1, const void *p2);

/*
 * select_by_genomic_location function selects rows in GTF_DATA that match with
 * the given genomic location (chr:start-end).
 *
 * Parameters:
 * 		gtf_data:	a GTF_DATA structure
 * 		chr:		the chromosome
 * 		begin_gl:	the start position on chromosome
 * 		end_gl:		the end position on chromosome
 *
 * Returns:			a GTF_DATA structure that contains the result of the query
 */
__attribute__ ((visibility ("default")))
GTF_DATA *select_by_genomic_location(GTF_DATA *gtf_data, char *chr, int begin_gl, int end_gl) {
	int i, start, end;
	ROW_LIST **find_row_list, *row_list, *test_row_list;

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
			 * For each row, get the start and end values (start and end
			 * columns are 3rd and 4th columns in GTF format)
			 */
			start = *((int *)gtf_data->data[row_list->row[i]]->data[3]);
			end = *((int *)gtf_data->data[row_list->row[i]]->data[4]);

			/*
			 * If start and end values of the row match with the given
			 * interval, add a GTF_ROW in the results
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
