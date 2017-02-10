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

SEQFRAG *insert_seqfrag(SEQFRAG *root, SEQFRAG *src, SEQFRAG *dest) {
	SEQFRAG *sf;

	//fprintf(stderr, "inserting (%d,%d) in (%d,%d)\n", src->start, src->end, dest->start, dest->end);

	if ((src->start - dest->start) > 0) {	//create seqfrag1
		//fprintf(stderr, "doing frag1 ...\n");
		sf = (SEQFRAG *)calloc(1, sizeof(SEQFRAG));
		sf->seqid = src->seqid;
		sf->start = dest->start;
		sf->end = src->start - 1;
		sf->strand = src->strand;
		sf->style = strdup(dest->style);
		sf->next = src;
		sf->previous = dest->previous;
		if (dest->previous != NULL) dest->previous->next = sf;
		src->previous = sf;
		if (root == dest) root = sf;
	}
	else {
		//fprintf(stderr, "no frag1\n");
		if (dest->previous != NULL) dest->previous->next = src;
		src->previous = dest->previous;
		if (root == dest) root = src;
	}
	//fprintf(stderr, "continue ...\n");
	if ((dest->end - src->end) > 0) {	//create seqfrag3
		//fprintf(stderr, "doing frag3 ...\n");
		sf = (SEQFRAG *)calloc(1, sizeof(SEQFRAG));
		sf->seqid = src->seqid;
		sf->start = src->end + 1;
		sf->end = dest->end;
		sf->strand = src->strand;
		sf->style = strdup(dest->style);
		sf->next = dest->next;
		sf->previous = src;
		src->next = sf;
		if (dest->next != NULL) dest->next->previous = sf;
	}
	else {
		//fprintf(stderr, "no frag3\n");
		src->next = dest->next;
		if (dest->next != NULL) dest->next->previous = src;
	}
	//fprintf(stderr, "finito !\n");
	return root;
}

SEQFRAG *make_seqfrag(char *seqid, int start, int end, char strand, char *style, char *color) {
	SEQFRAG * sf = (SEQFRAG *)calloc(1, sizeof(SEQFRAG));

	sf->seqid = seqid;
	sf->start = start;
	sf->end = end;
	sf->strand = strand;
	sf->style = (char *)calloc(strlen(style) + strlen(color) + 1, sizeof(char));
	strcpy(sf->style, style);
	strcat(sf->style, color);

	return sf;
}

static int compare_seqfrag(const void *p1, const void *p2) {
	SEQFRAG *e1 = (SEQFRAG *)p1;
	SEQFRAG *e2 = (SEQFRAG *)p2;
	if (e1->strand == '-')
		return e2->start - e1->start;
	return e1->start - e2->start;
}

static int compare_seqfrag_norc(const void *p1, const void *p2) {
	SEQFRAG *e1 = (SEQFRAG *)p1;
	SEQFRAG *e2 = (SEQFRAG *)p2;
	return e1->start - e2->start;
}

static int compare_seqchunk(const void *p1, const void *p2) {
	SEQ_CHUNK *e1 = (SEQ_CHUNK *)p1;
	SEQ_CHUNK *e2 = (SEQ_CHUNK *)p2;
	return e1->start - e2->start;
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
	print_attribute(data[row]->data[8], "gene_id", buffer + strlen(buffer), '_');
	print_attribute(data[row]->data[8], "gene_name", buffer + strlen(buffer), '_');
	print_attribute(data[row]->data[8], "transcript_id", buffer + strlen(buffer), '_');
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
	strcat(buffer, "_");
	free(tmp);
	tmp = column[6]->convert_to_string(data[row]->data[6], column[6]->default_value);
	strcat(buffer, tmp);
	free(tmp);
	if (rc && *(char *)(data[row]->data[6]) == '-') strcat(buffer, "_RC");
	if (!intron) strcat(buffer, "_mRNA");
	buffer = (char *)realloc(buffer, strlen(buffer) * sizeof(char));
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

/*__attribute__ ((visibility ("default")))
int get_fasta(FILE *output, GTF_DATA *gtf_data, char *genome_file, int intron, int rc, int color) {
	FILE *ff = NULL, *ffi;
	char **token, *seq = NULL, format[100], *seqid, strand;
	char *buffer = (char *)calloc(10000, sizeof(char));
	int i, k, maxLineSize = 0, n, nb_exon = 0, pfasta, tmp, start,end, pcdna, tr_len, nb;
	int f_toprint, f_i, f_reste, f_min;
	ENTRY item, *e;
	ROW_LIST *test_row_list = calloc(1, sizeof(ROW_LIST)), **find_row_list;
	SEQFRAG *seqfrag = NULL, *sf, *p_sf, *p_tmp;
	GTF_ROW *row;

	index_gtf(gtf_data, "transcript_id");

	nb = 0;
	ff = get_fasta_file(genome_file, buffer);
	//fprintf(stderr, "%s\n", buffer);
	if (ff != NULL) {
		ffi = get_fasta_file_index(ff, buffer);
		hdestroy();fprintf(out, ">");
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
		hcreate(100);
		while (fgets(buffer, 999, ffi) != NULL) {
			//fprintf(stderr, "%s\n", buffer);
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
		//fprintf(stderr, "COUCOU %d\n", gtf_data->size);
		for (k = 0; k < gtf_data->size; k++) {
			//fprintf(stderr, "coucou : %s %d %d\n", data[k]->data[2], (*(int *)(data[k]->data[3])), *(int *)(data[k]->data[4]));
			row = gtf_data->data[k];
			if (!strcmp(row->data[2], "transcript") || !strcmp(row->data[2], "ncRNA")) {
				nb++;
				print_header(output, k, gtf_data->data, intron, rc);
				item.key = row->data[0];
				e = hsearch(item, FIND);
				if (e != NULL) {
					tr_len = 0;
					test_row_list->token = get_attribute_value(row->data[8], "transcript_id", column[8]);
					//fprintf(stderr, "coucou : %s\n", test_row_list->token);
					find_row_list = (ROW_LIST **)tfind(test_row_list, &(column[8]->index[0]->data), compare_row_list);
					if (find_row_list != NULL) {
						nb_exon = 0;
						for (i = 0; i < (*find_row_list)->nb_row; i++)
							if (!strcmp((char *)(gtf_data->data[(*find_row_list)->row[i]]->data[2]), "exon"))
								nb_exon++;
						seqfrag = (SEQFRAG *)calloc(nb_exon, sizeof(SEQFRAG));
						nb_exon = 0;
						for (i = 0; i < (*find_row_list)->nb_row; i++)
							if (!strcmp(gtf_data->data[(*find_row_list)->row[i]]->data[2], "exon")) {
								seqfrag[nb_exon].start = *(int *)(gtf_data->data[(*find_row_list)->row[i]]->data[3]);
								seqfrag[nb_exon].end = *(int *)(gtf_data->data[(*find_row_list)->row[i]]->data[4]);
								seqfrag[nb_exon].seqid = (char *)(gtf_data->data[(*find_row_list)->row[i]]->data[0]);
								seqfrag[nb_exon].strand = *(char *)(gtf_data->data[(*find_row_list)->row[i]]->data[6]);
								seqfrag[nb_exon].style = (char *)calloc(strlen(plain) + strlen(black) + 1, sizeof(char));
								strcpy(seqfrag[nb_exon].style, plain);
								strcat(seqfrag[nb_exon].style, black);
								tr_len += seqfrag[nb_exon].end - seqfrag[nb_exon].start + 1;
								nb_exon++;
							}
						qsort(seqfrag, nb_exon, sizeof(SEQFRAG), rc && seqfrag[0].strand == '-' ? compare_seqfrag : compare_seqfrag_norc);
						if (intron) {
							tr_len = *(int *)(row->data[4]) - *(int *)(row->data[3]) + 1;
							seq = (char *)calloc(tr_len + 1, sizeof(char));
							get_chunk(seq, ff, *(long *)(e->data), maxLineSize, tr_len, *(int *)(row->data[3]), rc ? *(char *)(row->data[6]) : '+');
						}
						else {
							seq = (char *)calloc(tr_len + 1, sizeof(char));
							pcdna = 0;
							for (i = 0; i < nb_exon; i++) {
								get_chunk(seq + pcdna, ff, *(long *)(e->data), maxLineSize, seqfrag[i].end - seqfrag[i].start + 1, seqfrag[i].start, rc ? seqfrag[i].strand : '+');
								pcdna += seqfrag[i].end - seqfrag[i].start + 1;
							}
						}
					}
					if (!color)
						if (output == stdout)
							for (i = 0; i < tr_len; i += 60)
								fprintf(output, "%.60s\n", seq + i);
						else
							fprintf(output, "%s\n", seq);
					else {
						qsort(seqfrag, nb_exon, sizeof(SEQFRAG), compare_seqfrag_norc);
						for (i = 0; i < nb_exon - 1; i++) {
							if ((seqfrag[i + 1].start - seqfrag[i].end - 1) > 0 && intron) {
								sf = make_seqfrag(seqfrag[i].seqid, seqfrag[i].end + 1, seqfrag[i + 1].start - 1, seqfrag[i].strand, plain, gray);
								seqfrag[i].next = sf;
								sf->previous = &(seqfrag[i]);
								sf->next = &(seqfrag[i + 1]);
								seqfrag[i + 1].previous = sf;
							}
							else {
								seqfrag[i].next = &(seqfrag[i + 1]);
								seqfrag[i + 1].previous = &(seqfrag[i]);
							}
						}
						int tr_start = seqfrag[0].start;
						p_sf = seqfrag;
						while (p_sf != NULL) {
							p_sf->start -= tr_start;
							p_sf->end -= tr_start;
							if (p_sf->strand == '-' && rc) {
								tmp = p_sf->start;
								p_sf->start = tr_len - p_sf->end - 1;
								p_sf->end = tr_len - tmp - 1;
							}
							p_sf = p_sf->next;
						}
						if (rc && seqfrag[0].strand == '-') {
							p_sf = seqfrag;
							while (p_sf->next != NULL) p_sf = p_sf->next;
							seqfrag = p_sf;
							while (p_sf != NULL) {
								p_tmp = p_sf->previous;
								p_sf->previous = p_sf->next;
								p_sf->next = p_tmp;
								p_sf = p_sf->next;
							}
						}
						for (i = 0; i < (*find_row_list)->nb_row; i++) {
							seqid = (char *)(gtf_data->data[(*find_row_list)->row[i]]->data[0]);
							start = *(int *)(gtf_data->data[(*find_row_list)->row[i]]->data[3]) - tr_start;
							end = *(int *)(gtf_data->data[(*find_row_list)->row[i]]->data[4]) - tr_start;
							strand = *(char *)(gtf_data->data[(*find_row_list)->row[i]]->data[6]);
							if (strand == '-' && rc) {
								tmp = start;
								start = tr_len - end - 1;
								end = tr_len - tmp - 1;
							}
							sf = NULL;
							if (!strcmp(gtf_data->data[(*find_row_list)->row[i]]->data[2], "start_codon"))
								sf = make_seqfrag(seqid, start, end, strand, bold, green);
							else if (!strcmp(gtf_data->data[(*find_row_list)->row[i]]->data[2], "stop_codon"))
								sf = make_seqfrag(seqid, start, end, strand, bold, red);
							else if (!strcmp(gtf_data->data[(*find_row_list)->row[i]]->data[2], "five_prime_utr"))
								sf = make_seqfrag(seqid, start, end, strand, plain, yellow);
							else if (!strcmp(gtf_data->data[(*find_row_list)->row[i]]->data[2], "three_prime_utr"))
								sf = make_seqfrag(seqid, start, end, strand, plain, yellow);
							if (sf != NULL) {
								p_sf = seqfrag;
								while (p_sf != NULL)
									if ((sf->start >= p_sf->start) && (sf->end <= p_sf->end)) {
										seqfrag = insert_seqfrag(seqfrag, sf, p_sf);
										break;
									}
									else
										p_sf = p_sf->next;
							}
						}
						pfasta = 0;
						p_sf = seqfrag;
						if (!intron)
							while (p_sf != NULL) {
								tmp = p_sf->start;
								p_sf->start = pfasta;
								p_sf->end = pfasta + p_sf->end - tmp;
								pfasta += p_sf->end - p_sf->start + 1;
								p_sf = p_sf->next;
							}

						format[0] = '%';
						format[1] = '.';
						p_sf = seqfrag;
						f_reste = 0;
						while (p_sf != NULL) {
							f_i = 0;
							f_toprint = p_sf->end - p_sf->start + 1;
							fprintf(output, "%s", p_sf->style);
							if (f_reste > 0) {
								f_min = MIN(f_reste, f_toprint);
								sprintf(format + 2, "%d", f_min);
								strcat(format, "s");
								fprintf(output, format, seq + p_sf->start);
								f_reste -= f_min;
								f_toprint -= f_min;
								if (f_reste == 0) fprintf(output, "\n");
								f_i += f_min;
							}
							while (f_toprint > 0) {
								f_min = MIN(60, f_toprint);
								sprintf(format + 2, "%d", f_min);
								strcat(format, "s");
								fprintf(output, format, seq + p_sf->start + f_i);
								if (f_min == f_toprint) f_reste = 60 - f_toprint;
								f_toprint -= f_min;
								if (f_reste == 0) fprintf(output, "\n");
								f_i += f_min;
							}
							p_sf = p_sf->next;
						}
						fprintf(output, "\n");
						if (k != (gtf_data->size - 1)) fprintf(output, "%s%s", plain, black);
					}
				}
				free(seq);
			}
		}
	}
	return nb;
}*/

__attribute__ ((visibility ("default")))
SEQUENCES *get_fasta(GTF_DATA *gtf_data, char *genome_file, int intron, int rc) {
	SEQUENCES *ret = (SEQUENCES *)calloc(1, sizeof(SEQUENCES));
	SEQUENCE *sequence = NULL;
	SEQ_CHUNK *pchunk;
	FILE *ff = NULL, *ffi;
	char *buffer = (char *)calloc(10000, sizeof(char));
	int k;
	ENTRY item, *e;

	char **token, *seq = NULL, format[100], *seqid, strand;
	int i, maxLineSize = 0, n, nb_exon = 0, pfasta, tmp, start,end, pcdna, tr_len;
	int f_toprint, f_i, f_reste, f_min;
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
		ret->file_name = strdup(buffer);
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
				 * save the chromosome of the new sequence
				 */
				sequence->seqid = row->data[0];

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
						 * create the chunks for exons
						 */
						nb_exon = 0;
						pchunk = NULL;
						for (i = 0; i < (*find_row_list)->nb_row; i++)
							if (!strcmp((char *)(gtf_data->data[(*find_row_list)->row[i]]->data[2]), "exon")) {
								if (sequence->chunk == NULL)
									pchunk = sequence->chunk = (SEQ_CHUNK *)calloc(1, sizeof(SEQ_CHUNK));
								else {
									pchunk->next = (SEQ_CHUNK *)calloc(1, sizeof(SEQ_CHUNK));
									pchunk->next->previous = pchunk;
									pchunk = pchunk->next;
								}
								pchunk->start = *(int *)(gtf_data->data[(*find_row_list)->row[i]]->data[3]);
								pchunk->end = *(int *)(gtf_data->data[(*find_row_list)->row[i]]->data[4]);
								pchunk->type = EXON;
								tr_len += pchunk->end - pchunk->start + 1;
								nb_exon++;
							}
						//qsort(sequence->chunk, nb_exon, sizeof(SEQ_CHUNK), compare_seqchunk);

						/*
						 * add chunks for introns ... stay zen
						 */
						pchunk = sequence->chunk;
						while (pchunk->next != NULL) {
							pchunk->next->previous = (SEQ_CHUNK *)calloc(1, sizeof(SEQ_CHUNK));
							pchunk->next->previous->next = pchunk->next;
							pchunk->next->previous->previous = pchunk;
							pchunk->next = pchunk->next->previous;
							pchunk = pchunk->next->next;
							pchunk->previous->type = INTRON;
							pchunk->previous->start = pchunk->previous->previous->end + 1;
							pchunk->previous->end = pchunk->start - 1;;
						}

						/*int tr_start = seqfrag[0].start;
						p_sf = seqfrag;
						while (p_sf != NULL) {
							p_sf->start -= tr_start;
							p_sf->end -= tr_start;
							if (p_sf->strand == '-' && rc) {
								tmp = p_sf->start;
								p_sf->start = tr_len - p_sf->end - 1;
								p_sf->end = tr_len - tmp - 1;
							}
							p_sf = p_sf->next;
						}
						if (rc && seqfrag[0].strand == '-') {
							p_sf = seqfrag;
							while (p_sf->next != NULL) p_sf = p_sf->next;
							seqfrag = p_sf;
							while (p_sf != NULL) {
								p_tmp = p_sf->previous;
								p_sf->previous = p_sf->next;
								p_sf->next = p_tmp;
								p_sf = p_sf->next;
							}
						}
						for (i = 0; i < (*find_row_list)->nb_row; i++) {
							seqid = (char *)(gtf_data->data[(*find_row_list)->row[i]]->data[0]);
							start = *(int *)(gtf_data->data[(*find_row_list)->row[i]]->data[3]) - tr_start;
							end = *(int *)(gtf_data->data[(*find_row_list)->row[i]]->data[4]) - tr_start;
							strand = *(char *)(gtf_data->data[(*find_row_list)->row[i]]->data[6]);
							if (strand == '-' && rc) {
								tmp = start;
								start = tr_len - end - 1;
								end = tr_len - tmp - 1;
							}
							sf = NULL;
							if (!strcmp(gtf_data->data[(*find_row_list)->row[i]]->data[2], "start_codon"))
								sf = make_seqfrag(seqid, start, end, strand, bold, green);
							else if (!strcmp(gtf_data->data[(*find_row_list)->row[i]]->data[2], "stop_codon"))
								sf = make_seqfrag(seqid, start, end, strand, bold, red);
							else if (!strcmp(gtf_data->data[(*find_row_list)->row[i]]->data[2], "five_prime_utr"))
								sf = make_seqfrag(seqid, start, end, strand, plain, yellow);
							else if (!strcmp(gtf_data->data[(*find_row_list)->row[i]]->data[2], "three_prime_utr"))
								sf = make_seqfrag(seqid, start, end, strand, plain, yellow);
							if (sf != NULL) {
								p_sf = seqfrag;
								while (p_sf != NULL)
									if ((sf->start >= p_sf->start) && (sf->end <= p_sf->end)) {
										seqfrag = insert_seqfrag(seqfrag, sf, p_sf);
										break;
									}
									else
										p_sf = p_sf->next;
							}
						}
						pfasta = 0;
						p_sf = seqfrag;
						if (!intron)
							while (p_sf != NULL) {
								tmp = p_sf->start;
								p_sf->start = pfasta;
								p_sf->end = pfasta + p_sf->end - tmp;
								pfasta += p_sf->end - p_sf->start + 1;
								p_sf = p_sf->next;
							}
						format[0] = '%';
						format[1] = '.';
						p_sf = seqfrag;
						f_reste = 0;
						while (p_sf != NULL) {
							f_i = 0;
							f_toprint = p_sf->end - p_sf->start + 1;
							fprintf(output, "%s", p_sf->style);
							if (f_reste > 0) {
								f_min = MIN(f_reste, f_toprint);
								sprintf(format + 2, "%d", f_min);
								strcat(format, "s");
								fprintf(output, format, seq + p_sf->start);
								f_reste -= f_min;
								f_toprint -= f_min;
								if (f_reste == 0) fprintf(output, "\n");
								f_i += f_min;
							}
							while (f_toprint > 0) {
								f_min = MIN(60, f_toprint);
								sprintf(format + 2, "%d", f_min);
								strcat(format, "s");
								fprintf(output, format, seq + p_sf->start + f_i);
								if (f_min == f_toprint) f_reste = 60 - f_toprint;
								f_toprint -= f_min;
								if (f_reste == 0) fprintf(output, "\n");
								f_i += f_min;
							}
							p_sf = p_sf->next;
						}
						fprintf(output, "\n");
						if (k != (gtf_data->size - 1)) fprintf(output, "%s%s", plain, black);*/
					}
				}
				free(seq);
			}
		}
	}
	return ret;
}
