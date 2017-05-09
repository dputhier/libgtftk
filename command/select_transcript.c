/*
 * select_transcript.c
 *
 *  Created on: Feb 8, 2017
 *      Author: fafa
 *
 * The select_transcript function selects the transcripts based on several
 * criteria:
 * - the length of the mature transcript (shortest or longest)
 * - the most 5' transcript
 * The criteria is specified in the second argument of the function:
 * - 0 : shortest
 * - 1 : longest
 * - 2 : most 5'
 */

#include "libgtftk.h"

/*
 * external functions declaration
 */
extern INDEX_ID *index_gtf(GTF_DATA *gtf_data, char *key);
extern int comprow(const void *m1, const void *m2);
extern int compare_row_list(const void *p1, const void *p2);
extern int add_row_list(ROW_LIST *src, ROW_LIST *dst);
extern int add_row(int row, ROW_LIST *dst);

/*
 * global variables declaration
 */
extern COLUMN **column;

/*
 * We need some local variables because the research is made with the twalk
 * mechanism (tree browsing) in a separate function (action_sst) with
 * restricted arguments.
 * 	row_list:		aggregates all the selected rows (their rank)
 * 	test_row_list:	used to look for a transcript in the 2nd index of the
 * 					attributes column
 * 	find_row_list:	used to store the data related to the shortest transcript
 * 					of each gene
 * 	gtf_d:			a local copy of the GTF_DATA to process
 * 	tr_type:		the kind of transcript we want to get (shortest or longest)
 * 	tid_index:		the rank of the transcript_id index
 */
ROW_LIST *row_list, *test_row_list, **find_row_list;
GTF_DATA *gtf_d;
int tr_type;
INDEX_ID *tid_index;

/*
 * The function used by twalk on each node of the index tree.
 * This function is used to browse an index on "gene_id" attribute containing
 * ROW_LIST elements. For each gene, we compute the lengthes of his transcripts
 * and we keep the shortest and the longest ones.
 * For information about the parameters, see man pages of twalk.
 */
static void action_st(const void *nodep, const VISIT which, const int depth) {

	/*
	 * the feature value of a row
	 */
	char *feature;

	/*
	 * Variables used to find the shortest, the longest and the most 5'
	 * transcript of a gene.
	 * min_i, max_i:	the ranks of the first row of the shortest and the
	 * 					longest transcripts
	 * min_j, max_j:	the ranks of the transcript_id attribute in the min_i
	 * 					and max_i rows
	 * min_trsize:		the size of the shortest transcript
	 * max_trsize:		the size of the longest transcript
	 * most_5p:			the position of the most 5' transcript
	 * most_5p_i:		the rank of the first row of the most 5' transcript
	 * most_5p_j:		the ranks of the transcript_id attribute in the
	 * 					most_5p_i row
	 */
	int min_i = 0, min_j = 0, min_trsize;
	int max_i = 0, max_j = 0, max_trsize;
	int most_5p, most_5p_i = 0, most_5p_j = 0;

	// the size of the current transcript
	int trsize;

	// the coordinates of the current exon or transcript
	int start, end;

	// the "gene" row rank
	int r = -1;

	// loop variables
	int i, j, k;

	// the row list of the current gene
	ROW_LIST *datap = *((ROW_LIST **)nodep);

	// the current row
	GTF_ROW *row;

	switch (which) {
		case preorder:
			break;

		/*
		 * The operations are made on internal nodes and leaves of the tree.
		 */
		case leaf :
		case postorder:
			trsize = 0;
			min_trsize = 10000000;
			max_trsize = 0;
			most_5p = 0;

			/*
			 * loop on the rows of a gene
			 */
			for (i = 0; i < datap->nb_row; i++) {

				// the current row
				row = &gtf_d->data[datap->row[i]];

				// the current feature value
				feature = row->field[2];
				if (!strcmp(feature, "gene")) {
					/*
					 * if we parse a "gene" row, we have to keep in mind the
					 * rank of this row to print this "gene" row in the results
					 */
					r = datap->row[i];

					if (*(row->field[6]) == '+') most_5p = 300000000;

				}
				else if (!strcmp(feature, "transcript")) {
					/*
					 * if we parse a "transcript" row, we have to compute his
					 * length
					 */

					/*
					 * loop on the attributes of the current transcript, just to
					 * find the transcript_id value
					 */
					for (j = 0; j < row->nb_attributes; j++) {
						if (!strcmp(row->key[j], "transcript_id")) {
							/*
							 * once we know the transcript_id value, we use it
							 * to get all the rows related to it
							 */
							row_list->token = row->value[j];
							find_row_list = (ROW_LIST **)tfind(row_list, &(column[8]->index[tid_index->index_rank].data), compare_row_list);
							if (find_row_list != NULL) {
								trsize = 0;

								/*
								 * we loop on all exon rows of the transcript
								 * to compute his length
								 */
								for (k = 0; k < (*find_row_list)->nb_row; k++) {
									if (!strcmp(gtf_d->data[(*find_row_list)->row[k]].field[2], "exon")) {
										start = atoi(gtf_d->data[(*find_row_list)->row[k]].field[3]);
										end = atoi(gtf_d->data[(*find_row_list)->row[k]].field[4]);
										trsize += (end - start + 1);
									}
								}
								/*
								 * now we check if it is the current shortest
								 * transcript and, if it is the case, we store:
								 * 	the length
								 * 	the rank of the first row
								 * 	the rank of the transcript_id attribute
								 */
								if (trsize < min_trsize) {
									min_trsize = trsize;
									min_i = i;
									min_j = j;
								}

								/*
								 * the same for the longest transcript
								 */
								if (trsize > max_trsize) {
									max_trsize = trsize;
									max_i = i;
									max_j = j;
								}
							}

							/*
							 * now we check if it is the most 5' transcript and,
							 * if it is the case, we store:
							 * the position (most_5p)
							 * the rank of the first row (most_5p_i)
							 * the rank of the transcript_id attribute (most_5p_j)
							 */
							start = atoi(row->field[3]);
							end = atoi(row->field[4]);
							if (*(row->field[6]) == '+') {
								if (most_5p > start) {
									most_5p = start;
									most_5p_i = i;
									most_5p_j = j;
								}
							}
							else {
								if (most_5p < end) {
									most_5p = end;
									most_5p_i = i;
									most_5p_j = j;
								}
							}

							/*
							 * once we have found the transcript_id attribute,
							 * we don't need to check the other attributes
							 */
							break;
						}
					}
				}
			}
			/*
			 * after parsing the gene, min_i and max_i are the ranks of the
			 * first rows ("transcript" rows) of the shortest and longest
			 * transcripts of the gene and min_j and max_j are the ranks of the
			 * "transcript_id" attribute of min_i and max_i rows
			 */

			if (tr_type == SHORTEST_TRANSCRIPT) {
				/*
				 * we look for the shortest transcript
				 */
				test_row_list->token = gtf_d->data[datap->row[min_i]].value[min_j];
			}
			else if (tr_type == LONGEST_TRANSCRIPT) {
				/*
				 * we look for the longest transcript
				 */
				test_row_list->token = gtf_d->data[datap->row[max_i]].value[max_j];
			}
			else if (tr_type == MOST5P_TRANSCRIPT) {
				/*
				 * we look for the longest transcript
				 */
				test_row_list->token = gtf_d->data[datap->row[most_5p_i]].value[most_5p_j];
			}
			/*
			 * get transcript rows
			 */
			find_row_list = tfind(test_row_list, &(column[8]->index[tid_index->index_rank].data), compare_row_list);

			/*
			 * we add the gene and transcript rows in the row_list structure
			 */
			add_row_list(*find_row_list, row_list);
			add_row(r, row_list);
			break;

		case endorder:
			break;
	}
}

__attribute__ ((visibility ("default")))
GTF_DATA *select_transcript(GTF_DATA *gtf_data, int type) {
	int i, j;
	INDEX_ID *gid_index;

	/*
	 * we save the transcript type to allow action_st to access it
	 */
	tr_type = type;

	/*
	 * reserve memory for the GTF_DATA structure to return
	 */
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));

	/*
	 * indexes the GTF_DATA with gene_id and transcript_id attributes. The
	 * first index contains, for each gene, the list of his related rows.
	 * The rank of this index is 0. The second index contains, for each
	 * transcript, the list of his related rows. The rank of this index is 1.
	 */
	gid_index = index_gtf(gtf_data, "gene_id");
	tid_index = index_gtf(gtf_data, "transcript_id");

	/*
	 * setup local variables to allow action_sst function to access these
	 * values
	 */
	gtf_d = gtf_data;

	/*
	 * reserve memory for the local ROW_LIST
	 */
	row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
	test_row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));

	// tree browsing of the gene_id index (rank 0)
	twalk(column[8]->index[gid_index->index_rank].data, action_st);

	/*
	 * we sort the resulting row list to respect the original order of the
	 * transcripts
	 */
	qsort(row_list->row, row_list->nb_row, sizeof(int), comprow);

	/*
	 * now we fill the resulting GTF_DATA with the found rows and return it
	 */
	ret->data = (GTF_ROW *)calloc(row_list->nb_row, sizeof(GTF_ROW));
	for (i = 0; i < row_list->nb_row; i++) {
		ret->data[i].field = (char **)calloc(8, sizeof(char *));
		for (j = 0; j < 8; j++)
			ret->data[i].field[j] = strdup(gtf_data->data[row_list->row[i]].field[j]);

		ret->data[i].key = (char **)calloc(gtf_data->data[row_list->row[i]].nb_attributes, sizeof(char *));
		ret->data[i].value = (char **)calloc(gtf_data->data[row_list->row[i]].nb_attributes, sizeof(char *));
		for (j = 0; j < gtf_data->data[row_list->row[i]].nb_attributes; j++) {
			ret->data[i].key[j] = strdup(gtf_data->data[row_list->row[i]].key[j]);
			ret->data[i].value[j] = strdup(gtf_data->data[row_list->row[i]].value[j]);
		}
		//ret->data[i].key = gtf_data->data[row_list->row[i]].key;
		//ret->data[i].value = gtf_data->data[row_list->row[i]].value;
		//ret->data[i].field = gtf_data->data[row_list->row[i]].field;
		ret->data[i].nb_attributes = gtf_data->data[row_list->row[i]].nb_attributes;
		ret->data[i].rank = gtf_data->data[row_list->row[i]].rank;
	}
	ret->size = row_list->nb_row;
	return ret;
}
