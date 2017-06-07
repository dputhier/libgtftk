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
extern int update_linked_list(GTF_DATA *gtf_data);
//extern void print_row(FILE *output, GTF_ROW *r, char delim);

/*
 * global variables declaration
 */
extern COLUMN **column;
INDEX_ID *tid_index, *gid_index;
int n, nbrow;
GTF_DATA *gtf_d;

static void action_transcript(const void *nodep, const VISIT which, const int depth) {
	int i, ok, start, end, k, na;
	char *feature;
	GTF_ROW *row, *tr_row;

	// the row list of the current transcript
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

					k = atoi(row->field[3]);
					if (k < start) start = k;
					k = atoi(row->field[4]);
					if (k > end) end = k;

					if (!ok)
						if (!strcmp(feature, "exon")) {
							for (k = 0; k < row->nb_attributes; k++)
								if (strncmp(row->key[k], "exon", 4))
									tr_row->nb_attributes++;
							tr_row->key = (char **)calloc(tr_row->nb_attributes, sizeof(char *));
							tr_row->value = (char **)calloc(tr_row->nb_attributes, sizeof(char *));
							na = 0;
							for (k = 0; k < row->nb_attributes; k++)
								if (strncmp(row->key[k], "exon", 4)) {
									tr_row->key[na] = strdup(row->key[k]);
									tr_row->value[na] = strdup(row->value[k]);
									na++;
								}
							tr_row->field[0] = strdup(row->field[0]);
							tr_row->field[1] = get_attribute_value(row, "transcript_source");
							if (tr_row->field[1] == NULL) tr_row->field[1] = strdup(row->field[1]);
							tr_row->field[2] = strdup("transcript");
							tr_row->field[5] = strdup(row->field[5]);
							tr_row->field[6] = strdup(row->field[6]);
							tr_row->field[7] = strdup(row->field[7]);
							ok = 1;
							nbrow++;
						}
				}
				asprintf(&(tr_row->field[3]), "%d", start);
				asprintf(&(tr_row->field[4]), "%d", end);

				if (!strcmp(gtf_d->data[datap->row[0]]->field[2], "gene")) {
					tr_row->next = gtf_d->data[datap->row[1]];
					gtf_d->data[datap->row[0]]->next = tr_row;
				}
				else {
					tr_row->next = gtf_d->data[datap->row[0]];
					if (datap->row[0] != 0)
						gtf_d->data[datap->row[0] - 1]->next = tr_row;
					else
						gtf_d->data[0] = tr_row;
				}
			}
			n++;
			break;
		case endorder:
			break;
	}
}

static void action_gene(const void *nodep, const VISIT which, const int depth) {
	int i, ok, start, end, k, na;
	char *feature;
	GTF_ROW *row, *g_row;

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
				if (!strcmp(feature, "gene")) {
					ok = 1;
					break;
				}
			}
			if (!ok) {
				g_row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
				g_row->rank = -1;
				g_row->field = (char **)calloc(8, sizeof(char *));
				start = INT_MAX;
				end = 0;
				for (i = 0; i < datap->nb_row; i++) {
					// the current row
					row = gtf_d->data[datap->row[i]];

					// the current feature value
					feature = row->field[2];

					k = atoi(row->field[3]);
					if (k < start) start = k;
					k = atoi(row->field[4]);
					if (k > end) end = k;

					if (!strcmp(feature, "exon") || !strcmp(feature, "transcript")) {
						if (!ok) {
							for (k = 0; k < row->nb_attributes; k++)
								if (!strncmp(row->key[k], "gene", 4) || strstr(row->key[k], "_gene_") ||
									!strncmp(row->key[k] + strlen(row->key[k]) - 5, "_gene", 5))
									g_row->nb_attributes++;
							g_row->key = (char **)calloc(g_row->nb_attributes, sizeof(char *));
							g_row->value = (char **)calloc(g_row->nb_attributes, sizeof(char *));
							na = 0;
							for (k = 0; k < row->nb_attributes; k++) {
								if (!strncmp(row->key[k], "gene", 4) || strstr(row->key[k], "_gene_") ||
									!strncmp(row->key[k] + strlen(row->key[k]) - 5, "_gene", 5)) {
									g_row->key[na] = strdup(row->key[k]);
									g_row->value[na] = strdup(row->value[k]);
									na++;
								}
							}
							g_row->field[0] = strdup(row->field[0]);
							g_row->field[1] = get_attribute_value(row, "gene_source");
							if (g_row->field[1] == NULL) g_row->field[1] = strdup(row->field[1]);
							g_row->field[2] = strdup("gene");
							g_row->field[5] = strdup(row->field[5]);
							g_row->field[6] = strdup(row->field[6]);
							g_row->field[7] = strdup(row->field[7]);
							ok = 1;
							nbrow++;
						}
					}
				}
				asprintf(&(g_row->field[3]), "%d", start);
				asprintf(&(g_row->field[4]), "%d", end);

				if (datap->row[0] != 0)	gtf_d->data[datap->row[0] - 1]->next = g_row;
				g_row->next = gtf_d->data[datap->row[0]];
				if (datap->row[0] == 0)	gtf_d->data[0] = g_row;
			}
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
	GTF_DATA *inter = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));

	/*
	 * indexing the gtf with transcript_id
	 */
	tid_index = index_gtf(gtf_data, "transcript_id");

	/*
	 * tree browsing of the transcript_id index (rank 0)
	 */
	gtf_d = gtf_data;
	n = nbrow = 0;
	twalk(column[tid_index->column]->index[tid_index->index_rank].data, action_transcript);

	inter->data = (GTF_ROW **)calloc(1, sizeof(GTF_ROW *));
	GTF_ROW *row = gtf_d->data[0], *new_row, *previous_row = NULL;
	r = 0;
	while (row != NULL) {
		new_row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
		if (previous_row == NULL) inter->data[0] = new_row;
		new_row->field = (char **)calloc(8, sizeof(char *));
		for (j = 0; j < 8; j++) new_row->field[j] = strdup(row->field[j]);
		new_row->key = (char **)calloc(row->nb_attributes, sizeof(char *));
		new_row->value = (char **)calloc(row->nb_attributes, sizeof(char *));
		for (j = 0; j < row->nb_attributes; j++) {
			new_row->key[j] = strdup(row->key[j]);
			new_row->value[j] = strdup(row->value[j]);
		}
		new_row->nb_attributes = row->nb_attributes;
		new_row->rank = r++;
		if (previous_row != NULL) previous_row->next = new_row;
		previous_row = new_row;
		row = row->next;
	}
	inter->size = nbrow + gtf_d->size;
	update_row_table(inter);
	update_linked_list(gtf_data);

	/*
	 * indexing the gtf with gene_id
	 */
	gid_index = index_gtf(inter, "gene_id");

	/*
	 * tree browsing of the gene_id index (rank 1)
	 */
	gtf_d = inter;
	n = nbrow = 0;
	twalk(column[gid_index->column]->index[gid_index->index_rank].data, action_gene);

	ret->data = (GTF_ROW **)calloc(1, sizeof(GTF_ROW *));
	row = inter->data[0];
	previous_row = NULL;
	r = 0;
	while (row != NULL) {
		new_row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
		if (previous_row == NULL) ret->data[0] = new_row;
		new_row->field = (char **)calloc(8, sizeof(char *));
		for (j = 0; j < 8; j++) new_row->field[j] = strdup(row->field[j]);
		new_row->key = (char **)calloc(row->nb_attributes, sizeof(char *));
		new_row->value = (char **)calloc(row->nb_attributes, sizeof(char *));
		for (j = 0; j < row->nb_attributes; j++) {
			new_row->key[j] = strdup(row->key[j]);
			new_row->value[j] = strdup(row->value[j]);
		}
		new_row->nb_attributes = row->nb_attributes;
		new_row->rank = r++;
		if (previous_row != NULL) previous_row->next = new_row;
		previous_row = new_row;
		row = row->next;
	}
	ret->size = nbrow + inter->size;
	update_row_table(ret);
	update_linked_list(inter);

	return ret;
}
