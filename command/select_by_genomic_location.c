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
extern int add_row(int row, ROW_LIST *dst);

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
GTF_DATA *select_by_genomic_location(GTF_DATA *gtf_data, int nb_loc, char **chr, int *begin_gl, int *end_gl) {
	int i, j, start, end, k;
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

	// reserve memory for the final ROW_LIST
	row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));

	// reserve memory for the ROW_LIST used to search for chromosomes in index
	test_row_list = calloc(1, sizeof(ROW_LIST));

	/*
	 * Loop on the number of locations
	 */
	for (k = 0; k < nb_loc; k++) {

		/*
		 * Find all rows related to the given chromosome
		 */
		test_row_list->token = chr[k];
		find_row_list = (ROW_LIST **)tfind(test_row_list, &(column[0]->index[0]->data), compare_row_list);
		if (find_row_list != NULL) {
			for (j = 0; j < (*find_row_list)->nb_row; j++) {
				/*
				 * For each row, get the start and end values (start and end
				 * columns are 3rd and 4th columns in GTF format)
				 */
				start = *((int *)gtf_data->data[(*find_row_list)->row[j]]->data[3]);
				end = *((int *)gtf_data->data[(*find_row_list)->row[j]]->data[4]);

				/*
				 * If start and end values of the row match with the given
				 * interval, add a GTF_ROW in the results
				 */
				if ((begin_gl[k] >= start && begin_gl[k] <= end) ||
						(end_gl[k] >= start && end_gl[k] <= end) ||
						(begin_gl[k] <= start && end_gl[k] >= end))
					add_row((*find_row_list)->row[j], row_list);
			}
		}
	}

	/*
	 * now we fill the resulting GTF_DATA with the found rows and return it
	 */
	ret->data = (GTF_ROW **)calloc(row_list->nb_row, sizeof(GTF_ROW *));
	for (i = 0; i < row_list->nb_row; i++) ret->data[i] = gtf_data->data[row_list->row[i]];
	ret->size = row_list->nb_row;
	return ret;
}
