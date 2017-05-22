/*
 * convert_to_ensembl.c
 *
 *  Created on: May 18, 2017
 *      Author: fafa
 *
 *  The convert_to_ensembl function add "gene" and "transcript" rows if not
 *  present to conform to Ensembl GTF format
 */

#include "libgtftk.h"

/*
 * external functions declaration
 */
extern INDEX_ID *index_gtf(GTF_DATA *gtf_data, char *key);
extern char *get_attribute_value(GTF_ROW *row, char *attr);
extern int update_row_table(GTF_DATA *gtf_data);
//extern void print_row(FILE *output, GTF_ROW *r, char delim);

/*
 * global variables declaration
 */
extern COLUMN **column;
INDEX_ID *tid_index, *gid_index;
int n, nbrow;
GTF_DATA *gtf_d;

static void action_transcript(const void *nodep, const VISIT which, const int depth) {
	int i, ok, start, end, k, n;
	char *feature;
	GTF_ROW *row, *tr_row;;

	// the row list of the current gene
	ROW_LIST *datap = *((ROW_LIST **)nodep);

	switch (which) {
		case preorder:
			break;
		/*
		 * The operations are made on internal nodes and leaves of the tree.
		 */
		case leaf :
		case postorder:
			ok = 0;
			row = NULL;
			for (i = 0; i < datap->nb_row; i++) {
				// the current row
				row = gtf_d->data[datap->row[i]];

				// the current feature value
				feature = row->field[2];
				if (!strcmp(feature, "transcript")) {
					ok = 1;
					break;
				}
			}
			if (!ok) {
				fprintf(stderr, "missing transcript row for %s\n", get_attribute_value(row, "transcript_id"));
				tr_row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
				tr_row->rank = -1;
				tr_row->field = (char **)calloc(8, sizeof(char *));
				start = INT_MAX;
				end = 0;
				for (i = 0; i < datap->nb_row; i++) {
					// the current row
					row = gtf_d->data[datap->row[i]];

					// the current feature value
					feature = row->field[2];
					//fprintf(stderr, "feature = %s\n", feature);

					if (!strcmp(feature, "exon")) {
						k = atoi(row->field[3]);
						if (k < start) start = k;
						k = atoi(row->field[4]);
						if (k > end) end = k;
						//fprintf(stderr, "start = %d ; end = %d\n", start, end);
						if (!ok) {
							for (k = 0; k < row->nb_attributes; k++)
								if (strncmp(row->key[k], "exon", 4))
									tr_row->nb_attributes++;
							tr_row->key = (char **)calloc(tr_row->nb_attributes, sizeof(char *));
							tr_row->value = (char **)calloc(tr_row->nb_attributes, sizeof(char *));
							n = 0;
							for (k = 0; k < row->nb_attributes; k++) {
								if (strncmp(row->key[k], "exon", 4)) {
									tr_row->key[n] = strdup(row->key[k]);
									tr_row->value[n] = strdup(row->value[k]);
									n++;
								}
							}
							tr_row->field[0] = strdup(row->field[0]);
							tr_row->field[1] = strdup(row->field[1]);
							tr_row->field[2] = strdup("transcript");
							asprintf(&(tr_row->field[3]), "%d", start);
							asprintf(&(tr_row->field[4]), "%d", end);
							tr_row->field[5] = strdup(row->field[5]);
							tr_row->field[6] = strdup(row->field[6]);
							tr_row->field[7] = strdup(row->field[7]);
							nbrow++;
							ok = 1;
						}
					}
				}
				if (datap->row[0] != 0)	gtf_d->data[datap->row[0] - 1]->next = tr_row;
				tr_row->next = gtf_d->data[datap->row[0]];
				if (datap->row[0] == 0) gtf_d->data[0] = tr_row;
				//print_row(stderr, tr_row, '\t');
			}
			break;
		case endorder:
			break;
	}
}

static void action_gene(const void *nodep, const VISIT which, const int depth) {
	switch (which) {
		case preorder:
			break;
		/*
		 * The operations are made on internal nodes and leaves of the tree.
		 */
		case leaf :
		case postorder:
			n++;
			break;
		case endorder:
			break;
	}
}

__attribute__ ((visibility ("default")))
GTF_DATA *convert_to_ensembl(GTF_DATA *gtf_data) {
	int j, r;
	/*
	 * reserve memory for the GTF_DATA structure to return
	 */
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));

	/*
	 * indexing the gtf with transcript_id and gene_id
	 */
	//fprintf(stderr, "indexing gtf on transcript ids\n");
	tid_index = index_gtf(gtf_data, "transcript_id");
	//fprintf(stderr, "indexing gtf on gene ids\n");
	gid_index = index_gtf(gtf_data, "gene_id");
	//fprintf(stderr, "OK\n");

	gtf_d = gtf_data;

	/*
	 * tree browsing of the trancript_id index (rank 0)
	 */
	n = nbrow = 0;
	twalk(column[tid_index->column]->index[tid_index->index_rank].data, action_transcript);
	fprintf(stderr, "nb transcript : %d\n", n);

	/*
	 * tree browsing of the gene_id index (rank 0)
	 */
	n = 0;
	twalk(column[gid_index->column]->index[gid_index->index_rank].data, action_gene);
	fprintf(stderr, "nb gene : %d\n", n);

	ret->data = (GTF_ROW **)calloc(1, sizeof(GTF_ROW *));
	GTF_ROW *row = gtf_d->data[0], *new_row, *previous_row = NULL;
	r = 0;
	while (row != NULL) {
		new_row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
		if (row == gtf_d->data[0]) ret->data[0] = new_row;
		new_row->field = (char **)calloc(8, sizeof(char *));
		for (j = 0; j < 8; j++) new_row->field[j] = strdup(row->field[j]);
		new_row->key = (char **)calloc(row->nb_attributes, sizeof(char *));
		new_row->value = (char **)calloc(row->nb_attributes, sizeof(char *));
		for (j = 0; j < row->nb_attributes; j++) {
			new_row->key[j] = strdup(row->key[j]);
			new_row->value[j] = strdup(row->value[j]);
		}
		new_row->nb_attributes = row->nb_attributes;
		new_row->rank = r++; //row->rank;
		if (row != gtf_d->data[0]) previous_row->next = new_row;
		previous_row = new_row;
		row = row->next;
	}
	ret->size = nbrow + gtf_d->size;
	update_row_table(ret);
	return ret;
}
