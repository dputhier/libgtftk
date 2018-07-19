/*
 * load_blast.c
 *
 *  Created on: Jul 17, 2018
 *      Author: fafa
 */

#include "libgtftk.h"

extern GTF_READER *get_blast_reader(char *query);
extern char *get_blast_header(GTF_READER *, BLAST_HEADER *);
extern char *get_next_blast_hsp(GTF_READER *, BLAST_HSP *, char *);
extern void make_columns(void);
extern int update_row_table(GTF_DATA *gtf_data);

void new_attribute(char *key, char *value, GTF_ROW *row, BLAST_HSP *hsp) {
	row->attributes.nb++;
	row->attributes.attr = (ATTRIBUTE **)realloc(row->attributes.attr, row->attributes.nb * sizeof(ATTRIBUTE *));
	row->attributes.attr[row->attributes.nb - 1] = (ATTRIBUTE *)calloc(1, sizeof(ATTRIBUTE));
	row->attributes.attr[row->attributes.nb - 1]->key = strdup(key);
	row->attributes.attr[row->attributes.nb - 1]->value = strdup(value);
}

GTF_ROW *make_row(BLAST_HSP *hsp, GTF_DATA *gtf_data, int rank) {
	char *buffer = (char *)calloc(10000, sizeof(char));

	GTF_ROW *row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
	if (rank == 0) gtf_data->data[0] = row;
	row->field = (char **)calloc(8, sizeof(char *));
	row->field[0] = strdup(hsp->bs.subject_name);
	row->field[1] = strdup(hsp->bh.program_name);
	row->field[2] = strdup("HSP");
	sprintf(buffer, "%d", hsp->subject_start);
	row->field[3] = strdup(buffer);
	sprintf(buffer, "%d", hsp->subject_end);
	row->field[4] = strdup(buffer);
	sprintf(buffer, "%f", hsp->score);
	row->field[5] = strdup(buffer);
	row->field[6] = (char *)calloc(2, sizeof(char));
	*(row->field[6]) = hsp->strand_subject;
	row->field[7] = (char *)calloc(2, sizeof(char));
	*(row->field[7]) = '.';
	row->attributes.nb = 0;
	row->attributes.attr = NULL;
	new_attribute("database_name", hsp->bh.database_name, row, hsp);
	sprintf(buffer, "%u", hsp->bh.database_length);
	new_attribute("database_length", buffer, row, hsp);
	sprintf(buffer, "%d", hsp->bh.database_nb_sequences);
	new_attribute("database_nb_sequences", buffer, row, hsp);
	new_attribute("query_name", hsp->bq.query_name, row, hsp);
	sprintf(buffer, "%d", hsp->bq.query_length);
	new_attribute("query_length", buffer, row, hsp);
	sprintf(buffer, "%d", hsp->bs.subject_length);
	new_attribute("subject_length", buffer, row, hsp);
	sprintf(buffer, "%g", hsp->expect);
	new_attribute("expect", buffer, row, hsp);
	new_attribute("identities", hsp->identities, row, hsp);
	sprintf(buffer, "%d", hsp->identities_percent);
	new_attribute("identities_percent", buffer, row, hsp);
	if (hsp->gaps != NULL) {
		new_attribute("gaps", hsp->gaps, row, hsp);
		sprintf(buffer, "%d", hsp->gap_percent);
		new_attribute("gaps_percent", buffer, row, hsp);
	}
	row->rank = rank;

	free(buffer);
	return row;
}

__attribute__ ((visibility ("default")))
GTF_DATA *load_blast(char *input) {
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));
	ret->data = (GTF_ROW **)calloc(1, sizeof(GTF_ROW *));

	make_columns();

	GTF_READER *gr = get_blast_reader(input);
	BLAST_HSP *hsp = (BLAST_HSP *)calloc(1, sizeof(BLAST_HSP));
	char *r = get_blast_header(gr, &(hsp->bh));
	GTF_ROW *row, *previous_row;

	int n = 0;
	while ((r = get_next_blast_hsp(gr, hsp, r)) != NULL) {
		row = make_row(hsp, ret, n);
		if (n > 0) previous_row->next = row;
		previous_row = row;
		n++;
	}
	row = make_row(hsp, ret, n);
	if (n > 0) previous_row->next = row;
	previous_row = row;
	n++;
	free(r);
	ret->size = n;
	update_row_table(ret);
	return ret;
}
