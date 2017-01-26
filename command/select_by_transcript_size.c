/****
 * select_by_transcript_size.c
 *
 *  Created on: Jan 12, 2017
 *      Author: fafa
 */

#include "libgtftk.h"

GTF_DATA *gtf_d;
int min_ts, max_ts;

extern COLUMN **column;
extern void print_row(FILE *output, GTF_ROW *r, char delim);

static void action_sbts(const void *nodep, const VISIT which, const int depth) {
	ROW_LIST *datap;
	int i, trsize;

	switch (which) {
		case preorder:
			break;
		case postorder:
		case leaf:
			datap = *((ROW_LIST **)nodep);
			trsize = 0;
			for (i = 0; i < datap->nb_row; i++)
				if (!strcmp((char *)(gtf_d->data[datap->row[i]]->data[2]), "exon"))
					trsize += (*(int *)(gtf_d->data[datap->row[i]]->data[4]) - *(int *)(gtf_d->data[datap->row[i]]->data[3]) + 1);
			if ((trsize >= min_ts) && (trsize <= max_ts))
				for (i = 0; i < datap->nb_row; i++)
					print_row(stdout, gtf_d->data[datap->row[i]], '\t');
			break;
		case endorder:
			break;
	}
}

__attribute__ ((visibility ("default")))
GTF_DATA *select_by_transcript_size(GTF_DATA *gtf_data, int min, int max) {
	int i, found;
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));

	gtf_d = gtf_data;
	min_ts = min;
	max_ts = max;
	fprintf(stderr, "min = %d ; max = %d\n", min_ts, max_ts);
	found = 0;
	for (i = 0; i < column[8]->nb_index; i++)
		if (!strcmp(column[8]->index[i]->key, "transcript_id")) {
			found = 1;
			break;
		}
	fprintf(stderr, "found = %d ; i = %d ; id = %s\n", found, i, column[8]->index[i]->key);
	if (found)
		twalk(column[8]->index[i]->data, action_sbts);
	return ret;
}
