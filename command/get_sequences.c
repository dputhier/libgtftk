/*
 * get_fasta.c
 *
 *  Created on: Mar 30, 2016
 *      Author: fafa
 */

#include  "libgtftk.h"

extern void print_attribute(void *token, char *attr, char *output, char delim);
extern int split_ip(char ***tab, char *s, char *delim);
extern char *get_attribute_value(void *token, char *attr, void *r);
extern int compare_row_list(const void *p1, const void *p2);
extern int index_gtf(GTF_DATA *gtf_data, char *key);

//extern int color;
extern COLUMN **column;

void revcomp(char *s, int n) {
	int i;
	char c;

	for (i = 0; i < (n + 1) / 2; i++) {
		c = s[i];
		s[i] = COMPLEMENT(s[n - i - 1]);
		s[n - i - 1] = COMPLEMENT(c);
	}
}

void get_chunk(char *ret, FILE *fasta_file, long seqpos, int L, int N, int p, char strand) {
	int sr, sc, er, ec, reste_row = N, toget, reste_row_file, eof;

	sr = (p - 1) / L;
	sc = p - sr * L - 1;
	er = (p + N - 2) / L;
	ec = p + N - 2 - er * L;
	fseek(fasta_file, seqpos, SEEK_SET);
	if (strand == '+') {
		fseek(fasta_file, sr * (L + 1) + sc, SEEK_CUR);
		reste_row_file = L - sc;
		do {
			toget = MIN(reste_row, reste_row_file);
			eof = (fgets(ret + N - reste_row, toget + 1, fasta_file) == NULL);
			if (*(ret + strlen(ret) - 1) == '\n') *(ret + strlen(ret) - 1) = 0;
			reste_row -= toget;
			reste_row_file -= toget;
			if (!reste_row_file) {
				fgetc(fasta_file);
				reste_row_file = L;
			}
		} while (reste_row && !eof);
	}
	else {
		fseek(fasta_file, er * (L + 1) + ec, SEEK_CUR);
		reste_row_file = ec + 1;
		do {
			toget = MIN(reste_row, reste_row_file);
			fseek(fasta_file, 1 - toget, SEEK_CUR);
			fgets(ret + N - reste_row, toget + 1, fasta_file);
			revcomp(ret + N - reste_row, toget);
			reste_row -= toget;
			reste_row_file -= toget;
			fseek(fasta_file, -toget - 1, SEEK_CUR);
			if (!reste_row_file) {
				fseek(fasta_file, -1, SEEK_CUR);
				reste_row_file = L;
			}
		} while (reste_row);
	}
}

SEQFRAG *make_seqfrag(char *seqid, int start, int end, char strand, char *style, char *color) {
	SEQFRAG * sf = (SEQFRAG *)calloc(1, sizeof(SEQFRAG));

	sf->start = start;
	sf->end = end;
	sf->strand = strand;

	return sf;
}

static int compare_seqfrag(const void *p1, const void *p2) {
	SEQFRAG *e1 = (SEQFRAG *)p1;
	SEQFRAG *e2 = (SEQFRAG *)p2;
	return e1->start - e2->start;
}

static int compare_seqfrag_rc(const void *p1, const void *p2) {
	SEQFRAG *e1 = (SEQFRAG *)p1;
	SEQFRAG *e2 = (SEQFRAG *)p2;
	return e2->start - e1->start;
}

static int compare_feature_p(const void *p1, const void *p2) {
	FEATURE *e1 = *(FEATURE **)p1;
	FEATURE *e2 = *(FEATURE **)p2;
	if (e1->start != e2->start)
		return e1->start - e2->start;
	if (!strcmp(e1->name, "exon"))
		return -1;
	return 1;
}

static int compare_feature_m(const void *p1, const void *p2) {
	FEATURE *e1 = *(FEATURE **)p1;
	FEATURE *e2 = *(FEATURE **)p2;
	if (e1->end != e2->end)
		return e2->end - e1->end;
	if (!strcmp(e1->name, "exon"))
		return -1;
	return 1;
}

static int compare_feature_tr(const void *p1, const void *p2) {
	FEATURE *e1 = *(FEATURE **)p1;
	FEATURE *e2 = *(FEATURE **)p2;
	if (e1->tr_start != e2->tr_start)
		return e1->tr_start - e2->tr_start;
	if (!strcmp(e1->name, "exon"))
		return -1;
	return 1;
}

/*void print_header(FILE *out, int row, GTF_ROW **data, int intron, int rc) {
	fprintf(out, ">");
	print_attribute(data[row]->data[8], "gene_id", out, '_');
	print_attribute(data[row]->data[8], "gene_name", out, '_');
	print_attribute(data[row]->data[8], "transcript_id", out, '_');
	column[0]->print(data[row]->data[0], out, column[0], ':');
	column[3]->print(data[row]->data[3], out, column[3], '-');
	column[4]->print(data[row]->data[4], out, column[4], '_');
	column[6]->print(data[row]->data[6], out, column[6], 0);
	if (rc && *(char *)(data[row]->data[6]) == '-') fprintf(out, "_RC");
	if (!intron) fprintf(out, "_mRNA");
	fprintf(out, "\n");
}*/

char *make_header(int row, GTF_ROW **data, int intron, int rc) {
	char *buffer = (char *)calloc(1000, sizeof(char));
	char *tmp;

	strcat(buffer, ">");
	print_attribute(data[row]->data[8], "gene_id", buffer + strlen(buffer), '`');
	print_attribute(data[row]->data[8], "gene_name", buffer + strlen(buffer), '`');
	print_attribute(data[row]->data[8], "transcript_id", buffer + strlen(buffer), '`');
	tmp = column[0]->convert_to_string(data[row]->data[0], column[0]->default_value);
	strcat(buffer, tmp);
	strcat(buffer, ":");
	free(tmp);
	tmp = column[3]->convert_to_string(data[row]->data[3], column[3]->default_value);
	strcat(buffer, tmp);
	strcat(buffer, "-");
	free(tmp);
	tmp = column[4]->convert_to_string(data[row]->data[4], column[4]->default_value);
	strcat(buffer, tmp);
	strcat(buffer, "`");
	free(tmp);
	tmp = column[6]->convert_to_string(data[row]->data[6], column[6]->default_value);
	strcat(buffer, tmp);
	free(tmp);
	if (rc && *(char *)(data[row]->data[6]) == '-') strcat(buffer, "`RC");
	if (!intron) strcat(buffer, "`mRNA");
	buffer = (char *)realloc(buffer, (strlen(buffer) + 1) * sizeof(char));
	return buffer;
}

FILE *get_fasta_file(char *fasta_file, char *buffer) {
	FILE *ff = NULL;
	char *fname;

	if (!access(fasta_file, F_OK)) {
		ff = fopen(fasta_file, "ro");
		fname = strrchr(fasta_file, '/');
		if (fname != NULL)
			fname++;
		else
			fname = fasta_file;
		strcpy(buffer, getenv("HOME"));
		strcat(buffer, "/.gtftk/");
		strcat(buffer, fname);
	}
	else {
		strcpy(buffer, getenv("HOME"));
		strcat(buffer, "/.gtftk/");
		strcat(buffer, fasta_file);
		if (!access(buffer, F_OK)) ff = fopen(buffer, "ro");
	}
	return ff;
}

FILE *get_fasta_file_index(FILE *fasta_file, char *index) {
	FILE *ffi = NULL;
	long pfasta;
	char *buffer = (char *)calloc(1000, sizeof(char));
	int maxLineSize = 0;

	if (access(index, F_OK)) {
		ffi = fopen(index, "w+");
		pfasta = ftell(fasta_file);
		while (fgets(buffer, 999, fasta_file) != NULL) {
			if (*buffer == '>') {
				*(buffer + strlen(buffer) - 1) = 0;
				fprintf(ffi, "%s\t%ld\t%ld\n", buffer + 1, pfasta, ftell(fasta_file));
			}
			else
				maxLineSize = MAX(maxLineSize, strlen(buffer));
			pfasta = ftell(fasta_file);
		}
		fprintf(ffi, "%d\n", maxLineSize - 1);
		rewind(fasta_file);
		fflush(ffi);
		rewind(ffi);
	}
	else
		ffi = fopen(index, "r");
	return ffi;
}

char *get_attribute_value(void *token, char *attr, void *r) {
	int k;
	if (((ATTRIBUTES *)token)->nb > 0)
		for (k = 0; k < ((ATTRIBUTES *)token)->nb; k++)
			if (!strcmp(attr, ((ATTRIBUTES *)token)->attr[k]->key))
				return ((ATTRIBUTES *)token)->attr[k]->value;
	return NULL;
}

void print_fasta_sequence(SEQUENCE *seq) {
	int k;
	FEATURE *feat;

	fprintf(stdout, "%s\n", seq->header);
	for (k = 0; k < strlen(seq->sequence); k += 60)
		fprintf(stdout, "%.60s\n", seq->sequence + k);
	for (k = 0; k < seq->features->nb; k++) {
		feat = seq->features->feature[k];
		fprintf(stdout, "  %s : %d-%d (%d-%d)\n", feat->name, feat->start, feat->end, feat->tr_start, feat->tr_end);
	}
}

/*
 * This function retrieve the nucleotide sequences of the transcripts
 * described by GTF data. For each transcript, the function provides also
 * all the features coordinates, in the genome and in the transcript.
 *
 * Parameters:
 * 		gtf_data:		the GTF data describing the transcripts
 * 		genome_file:	the fasta file containing the genome from which to
 * 						extract the transcript sequences
 * 		intron:			1 to keep the introns on the transcript sequences
 * 		rc:				1 to get reverse-complemented sequences for transcripts
 * 						located on minus strand
 *
 * Return:
 * 		a pointer on a SEQUENCES structure
 */
__attribute__ ((visibility ("default")))
SEQUENCES *get_sequences(GTF_DATA *gtf_data, char *genome_file, int intron, int rc) {
	SEQUENCES *ret = (SEQUENCES *)calloc(1, sizeof(SEQUENCES));
	SEQUENCE *sequence = NULL;
	FEATURE *feat, *new_feat;
	FILE *ff = NULL, *ffi;
	char *buffer = (char *)calloc(10000, sizeof(char));
	int j, k, pb, pe, pf, start, end, tr_start, tr_end;
	ENTRY item, *e;
	SEQFRAG *seqfrag;

	char **token, *feature;
	int i, n, nb_exon = 0, tr_len, maxLineSize = 0, pcdna;
	ROW_LIST *test_row_list = calloc(1, sizeof(ROW_LIST)), **find_row_list;
	GTF_ROW *row;

	/*
	 * Indexes GTF data with transcripts IDs
	 */
	index_gtf(gtf_data, "transcript_id");

	/*
	 * Test if genome file exists
	 */
	ff = get_fasta_file(genome_file, buffer);
	if (ff != NULL) {
		strcat(buffer, ".gtftk");

		/*
		 * gets the index file of the genome file. If it doesn't exist, creates
		 * one in ~/.gtftk directory with the same name as genome file plus
		 * ".gtftk" at the end.
		 */
		ffi = get_fasta_file_index(ff, buffer);

		/*
		 * From index file, creates a hastable with :
		 * 	key = chromosome name
		 * 	value = file position of the start of the sequence
		 */
		hdestroy();
		hcreate(100);
		while (fgets(buffer, 999, ffi) != NULL) {
			n = split_ip(&token, buffer, " \t\n");
			if (n > 1) {
				item.key = strdup(token[0]);
				item.data = (long *)malloc(sizeof(long));
				*(long *)item.data = atol(token[n - 1]);
				hsearch(item, ENTER);
			}
			else
				maxLineSize = atoi(token[0]);
			free(token);
		}
		fclose(ffi);

		/*
		 * The main loop on GTF rows
		 */
		for (k = 0; k < gtf_data->size; k++) {
			row = gtf_data->data[k];
			/*
			 * get only transcripts and non coding RNAs (transcripted sequences)
			 */
			if (!strcmp(row->data[2], "transcript") || !strcmp(row->data[2], "ncRNA")) {
				/*
				 * Create a new SEQUENCE and add it in the results table
				 */
				sequence = (SEQUENCE *)calloc(1, sizeof(SEQUENCE));
				ret->sequence = (SEQUENCE **)realloc(ret->sequence, (ret->nb + 1) * sizeof(SEQUENCE *));
				ret->sequence[ret->nb] = sequence;
				ret->nb++;

				/*
				 * Make the header of the new sequence
				 */
				sequence->header = make_header(k, gtf_data->data, intron, rc);

				/*
				 * save the strand
				 */
				sequence->strand = *(char *)(row->data[6]);

				/*
				 * save the start of the transcript
				 */
				sequence->start = *(int *)(row->data[3]);
				sequence->end = *(int *)(row->data[4]);

				/*
				 * save the sequence id (chromosome)
				 */
				sequence->seqid = strdup((char *)(row->data[0]));

				/*
				 * Look for chromosome in the hashtable made from genome index
				 * file
				 */
				item.key = row->data[0];
				e = hsearch(item, FIND);
				if (e != NULL) {
					/*
					 * Search for transcript ID in the index
					 */
					test_row_list->token = get_attribute_value(row->data[8], "transcript_id", column[8]);
					find_row_list = (ROW_LIST **)tfind(test_row_list, &(column[8]->index[0]->data), compare_row_list);

					tr_len = 0;
					if (find_row_list != NULL) {
						/*
						 * find the number of exons of the current transcript
						 */
						nb_exon = 0;
						for (i = 0; i < (*find_row_list)->nb_row; i++)
							if (!strcmp((char *)(gtf_data->data[(*find_row_list)->row[i]]->data[2]), "exon"))
								nb_exon++;

						/*
						 * reserve memory for the table of exons
						 */
						seqfrag = (SEQFRAG *)calloc(nb_exon, sizeof(SEQFRAG));

						/*
						 * fill the table of exons and compute the transcript
						 * size (sum of exon sizes)
						 */
						nb_exon = 0;
						for (i = 0; i < (*find_row_list)->nb_row; i++)
							if (!strcmp(gtf_data->data[(*find_row_list)->row[i]]->data[2], "exon")) {
								seqfrag[nb_exon].start = *(int *)(gtf_data->data[(*find_row_list)->row[i]]->data[3]);
								seqfrag[nb_exon].end = *(int *)(gtf_data->data[(*find_row_list)->row[i]]->data[4]);
								seqfrag[nb_exon].strand = *(char *)(gtf_data->data[(*find_row_list)->row[i]]->data[6]);
								tr_len += seqfrag[nb_exon].end - seqfrag[nb_exon].start + 1;
								nb_exon++;
							}
						/*
						 * sort exons by location on chromosome, in decreasing
						 * order if transcript is on minus strand and user
						 * asked for reverse-complemented sequence (mRNA). This
						 * is necessary for the next step : getting sequence
						 * chunks in right order
						 */
						qsort(seqfrag, nb_exon, sizeof(SEQFRAG),
								rc && seqfrag[0].strand == '-' ? compare_seqfrag_rc : compare_seqfrag);

						/*
						 * If introns are to be kept in the sequence, the
						 * sequence length is computed from the start and end
						 * fields of the "transcript" row. We just have to call
						 * get_chunk one time to get the sequence from the fasta
						 * indexed file.
						 */
						if (intron) {
							tr_len = *(int *)(row->data[4]) - *(int *)(row->data[3]) + 1;
							sequence->sequence = (char *)calloc(tr_len + 1, sizeof(char));
							get_chunk(sequence->sequence, ff, *(long *)(e->data), maxLineSize, tr_len, *(int *)(row->data[3]), rc ? *(char *)(row->data[6]) : '+');
						}

						/*
						 * If we want only exons (to get messenger sequence),
						 * we use the transcript size computed before, we call
						 * get_chunk for each exon and we concatenate the exons
						 * by adding pcdna (the actual length of the sequence)
						 * in the first argument of get_chunk.
						 */
						else {
							sequence->sequence = (char *)calloc(tr_len + 1, sizeof(char));
							pcdna = 0;
							for (i = 0; i < nb_exon; i++) {
								get_chunk(sequence->sequence + pcdna, ff, *(long *)(e->data), maxLineSize, seqfrag[i].end - seqfrag[i].start + 1, seqfrag[i].start, rc ? seqfrag[i].strand : '+');
								pcdna += seqfrag[i].end - seqfrag[i].start + 1;
							}
						}

						/*
						 * Make the features
						 */
						sequence->features = (FEATURES *)calloc(1, sizeof(FEATURES));
						for (i = 0; i < (*find_row_list)->nb_row; i++) {
							feature = gtf_data->data[(*find_row_list)->row[i]]->data[2];
							if (strcmp(feature, "transcript") && strcmp(feature, "ncRNA")) {
								sequence->features->feature = (FEATURE **)realloc(sequence->features->feature, (sequence->features->nb + 1) * sizeof(FEATURE *));
								feat = sequence->features->feature[sequence->features->nb] = (FEATURE *)calloc(1, sizeof(FEATURE));
								sequence->features->nb++;
								feat->name = strdup(feature);
								feat->start = *(int *)(gtf_data->data[(*find_row_list)->row[i]]->data[3]);
								feat->end = *(int *)(gtf_data->data[(*find_row_list)->row[i]]->data[4]);
							}
						}

						/*
						 * sort the features by genomic location, in decreasing
						 * order if transcript is on minus strand. This second
						 * sorting is made for adding intron features.
						 */
						qsort(sequence->features->feature, sequence->features->nb, sizeof(FEATURE *),
								seqfrag[0].strand == '-' ? compare_feature_m : compare_feature_p);

						/*
						 * Add introns in features if needed
						 */
						if (intron) {
							if (sequence->strand == '-') {
								pb = 0;
								for (j = 0; j < sequence->features->nb; j++) {
									feat = sequence->features->feature[j];
									if (!strcmp(feat->name, "exon")) {
										if (pb == 0)
											pb = feat->start;
										else {
											sequence->features->feature = (FEATURE **)realloc(sequence->features->feature, (sequence->features->nb + 1) * sizeof(FEATURE *));
											new_feat = sequence->features->feature[sequence->features->nb] = (FEATURE *)calloc(1, sizeof(FEATURE));
											sequence->features->nb++;
											new_feat->name = strdup("intron");
											new_feat->start = feat->end + 1;
											new_feat->end = pb - 1;
											pb = feat->start;
										}
									}
								}
							}
							else {
								pe = 0;
								for (j = 0; j < sequence->features->nb; j++) {
									feat = sequence->features->feature[j];
									if (!strcmp(feat->name, "exon")) {
										if (pe == 0)
											pe = feat->end;
										else {
											sequence->features->feature = (FEATURE **)realloc(sequence->features->feature, (sequence->features->nb + 1) * sizeof(FEATURE *));
											new_feat = sequence->features->feature[sequence->features->nb] = (FEATURE *)calloc(1, sizeof(FEATURE));
											sequence->features->nb++;
											new_feat->name = strdup("intron");
											new_feat->start = pe + 1;
											new_feat->end = feat->start - 1;
											pe = feat->end;
										}
									}
								}
							}
						}

						/*
						 * Now that we have all features, we can compute their
						 * local coordinates, i.e. the coordinates in the
						 * transcript (tr_start and tr_end).
						 */
						if (intron) {
							/*
							 * If we keep introns, we just have to re-adjust the
							 * coordinates relatively to the end or the start of
							 * the transcript, depending of the strand (minus
							 * or plus).
							 */
							qsort(sequence->features->feature, sequence->features->nb, sizeof(FEATURE *),
									seqfrag[0].strand == '-' ? compare_feature_m : compare_feature_p);
							if (rc && seqfrag[0].strand == '-')
								for (j = 0; j < sequence->features->nb; j++) {
									feat = sequence->features->feature[j];
									feat->tr_start = sequence->end - feat->end;
									feat->tr_end = sequence->end - feat->start;
								}
							else
								for (j = 0; j < sequence->features->nb; j++) {
									feat = sequence->features->feature[j];
									feat->tr_start = feat->start - sequence->start;
									feat->tr_end = feat->end - sequence->start;
								}
						}
						else {
							/*
							 * If we don't keep introns, the operation is more
							 * complicated because we have to base the
							 * computation on exons to determine the correct
							 * coordinates of all other features.
							 */
							qsort(sequence->features->feature, sequence->features->nb, sizeof(FEATURE *),
									sequence->strand == '-' && rc ? compare_feature_m : compare_feature_p);
							pf = start = end = tr_start = tr_end = 0;
							for (j = 0; j < sequence->features->nb; j++) {
								feat = sequence->features->feature[j];
								if (!strcmp(feat->name, "exon")) {
									tr_start = feat->tr_start = pf;
									tr_end = feat->tr_end = pf + feat->end - feat->start;
									pf += feat->end - feat->start + 1;
									start = feat->start;
									end = feat->end;
								}
								else {
									if (sequence->strand == '-' && rc) {
										if ((feat->start >= start) && (feat->start <= end))
											feat->tr_start = end - feat->end + tr_start;
										else
											feat->tr_start = -1;
										if ((feat->end >= start) && (feat->end <= end))
											feat->tr_end = end - feat->start + tr_start;
										else
											feat->tr_end = -1;
									}
									else {
										if ((feat->start >= start) && (feat->start <= end))
											feat->tr_start = feat->start - start + tr_start;
										else
											feat->tr_start = -1;
										if ((feat->end >= start) && (feat->end <= end))
											feat->tr_end = feat->end - start + tr_start;
										else
											feat->tr_end = -1;
									}
								}
							}
						}
					}
				}

				/*
				 * The last sort operation on features, to ensure that they are
				 * in ascending order, based on their local start coordinate in
				 * the transcript.
				 */
				qsort(sequence->features->feature, sequence->features->nb, sizeof(FEATURE *), compare_feature_tr);
			}
		}
	}
	return ret;
}
