/*
 * select_by_number_of_exon.c
 *
 *  Created on: Mar 30, 2016
 *      Author: fafa
 *
 * Implementation of select_by_number_of_exons function.
 * This function select rows in GTF_DATA that correspond to transcripts with a
 * number of exons between the given min and max values.
 */

#include "libgtftk.h"

/*
 * external functions declaration
 */
extern int index_gtf(GTF_DATA *gtf_data, char *key);
extern int comprow(const void *m1, const void *m2);
extern int add_row_list(ROW_LIST *src, ROW_LIST *dst);

/*
 * global variables declaration
 */
extern COLUMN **column;

/*
 * We need some local variables because the research is made with the twalk
 * mechanism (tree browsing) in a separate function (action_sbnoe) with
 * restricted arguments.
 * 	row_list:			a ROW_LIST to aggregate all the selected rows (their rank)
 * 	gtf_d:				a local copy of the GTF_DATA to process
 * 	min_noe, max_noe:	the local copies of min and max number of exons values
 */
ROW_LIST *row_list;
GTF_DATA *gtf_d;
int min_noe, max_noe;

/*
 * The comparison function used by twalk.
 * This function is used to browse an index on "transcript_id" attribute
 * containing ROW_LIST elements. For each transcript, we compute the number of
 * exons and if this number is between min and max values, the associated
 * ROW_LIST element is merged with the local row_list.
 * For information about the parameters, see man pages of twalk.
 */
static void action_sbnoe(const void *nodep, const VISIT which, const int depth) {
	ROW_LIST *datap;
	GTF_ROW *row;
	int i, noe;

	switch (which) {
		case preorder:
			break;

		case leaf:
		case postorder:
			datap = *((ROW_LIST **)nodep);
			noe = 0;
			for (i = 0; i < datap->nb_row; i++) {
				row = &gtf_d->data[datap->row[i]];
				if (!strcmp(row->field[2], "exon")) noe++;
			}
			if ((noe >= min_noe) && (noe <= max_noe)) add_row_list(datap, row_list);
			break;

		case endorder:
			break;
	}
}

/*
 * select_by_number_of_exon function selects rows in GTF_DATA that correspond
 * to transcripts that contains a number of exons between the given min and max
 * values.
 *
 * Parameters:
 * 		gtf_data:	a GTF_DATA structure
 * 		min:		the minimum number of exons
 * 		max:		the maximum number of exons
 *
 * Returns:			a GTF_DATA structure that contains the result of the query
 */
__attribute__ ((visibility ("default")))
GTF_DATA *select_by_number_of_exon(GTF_DATA *gtf_data, int min, int max) {
	int i;

	/*
	 * reserve memory for the GTF_DATA structure to return
	 */
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));

	/*
	 * reserve memory for the local ROW_LIST
	 */
	row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));

	/*
	 * setup local variables to allow action_sbnoe function to access these
	 * values
	 */
	gtf_d = gtf_data;
	min_noe = min;
	max_noe = max;

	/*
	 * indexes the GTF_DATA with transcript_id attribute to create an index
	 * containing, for each transcript, the list of his related rows. The real
	 * rank of the index (- 8) in the attributes column index table is stored
	 * in i.
	 */
	i = index_gtf(gtf_data, "transcript_id") - 8;

	// tree browsing of the transcript_id index
	twalk(column[8]->index[i]->data, action_sbnoe);

	/*
	 * we sort the resulting row list to be sure to respect the original order
	 * of the transcripts
	 */
	qsort(row_list->row, row_list->nb_row, sizeof(int), comprow);

	/*
	 * now we fill the resulting GTF_DATA with the found rows and return it
	 */
	ret->data = (GTF_ROW *)calloc(row_list->nb_row, sizeof(GTF_ROW));
	for (i = 0; i < row_list->nb_row; i++) {
		ret->data[i].field = gtf_data->data[row_list->row[i]].field;
		ret->data[i].key = gtf_data->data[row_list->row[i]].key;
		ret->data[i].value = gtf_data->data[row_list->row[i]].value;
		ret->data[i].nb_attributes = gtf_data->data[row_list->row[i]].nb_attributes;
		ret->data[i].rank = gtf_data->data[row_list->row[i]].rank;
	}
	ret->size = row_list->nb_row;

	return ret;
}
