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

/*
 * global variables declaration
 */
extern COLUMN **column;
INDEX_ID *tid_index, *gid_index;
int n;

static void action_transcript(const void *nodep, const VISIT which, const int depth) {
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

	/*
	 * tree browsing of the trancript_id index (rank 0)
	 */
	n = 0;
	twalk(column[tid_index->column]->index[tid_index->index_rank].data, action_transcript);
	fprintf(stderr, "nb transcript : %d\n", n);

	/*
	 * tree browsing of the gene_id index (rank 0)
	 */
	n = 0;
	twalk(column[gid_index->column]->index[gid_index->index_rank].data, action_gene);
	fprintf(stderr, "nb gene : %d\n", n);

	return ret;
}
