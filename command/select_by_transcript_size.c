/****
 * select_by_transcript_size.c
 *
 *  Created on: Jan 12, 2017
 *      Author: fafa
 */

#include <libgtftk.h>

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
			/*for (i = 0; i < datap->nb_row; i++)
				if (!strcmp((char *)(data[datap->row[i]]->data[2]), "exon"))
					trsize += (*(int *)(data[datap->row[i]]->data[4]) - *(int *)(data[datap->row[i]]->data[3]) + 1);
			if ((trsize >= min_ts) && (trsize <= max_ts))
				for (i = 0; i < datap->nb_row; i++)
					print_row(stdout, data[datap->row[i]], '\t', &column, nb_column);*/
			break;
		case endorder:
			break;
	}
}

__attribute__ ((visibility ("default")))
char **select_by_transcript_size(GTF_DATA *gtf_data, int min, int max) {
	/*twalk(column[8]->index[1]->data, action_sbts);*/
	return NULL;
}
