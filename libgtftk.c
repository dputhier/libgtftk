/*
 * gtftk_lib.c
 *
 *  Created on: Apr 1, 2016
 *      Author: fafa
 */

#include "libgtftk.h"

extern void print_row(FILE *output, GTF_ROW *r, char delim);
//extern int extract_data(FILE *output, char *gtf_filename, char *key, COLUMN **column, ROW **data, int nb_row, int nb_column);
//extern int select_by_genomic_location(FILE *output, char *gtf_filename, char *chr, int start, int end, COLUMN **column, ROW **data, int nb_column);
//extern int get_fasta(FILE *output, char *gtf_filename, char *fasta_file, int intron, int rc, COLUMN **column, ROW **data, int nb_row);
//extern void init_gtf(char *gtf_filename, char *key, COLUMN ***column, int *nb_column, ROW ***p_data, int *nb_row);
//extern void init_gtf_from_gtfresult(GTF_RESULT *gtf_result, char *key, COLUMN ***column, int *nb_column, ROW ***p_data, int *nb_row);
//extern int get_feature_list(FILE *output);
//extern int get_attribute_list(FILE *output);
//extern int get_seq_list(FILE *output);
//extern void init_ftp_data();
//extern char **select_by_transcript_size(char *gtf_filename, char *key, char *value, int not, COLUMN **column, ROW **data);
//extern FILE *open_memstream(char **ptr, size_t *sizeloc);
//extern FILE *fmemopen(void *buf, size_t size, const char *mode);

COLUMN **column;
int nb_column;

__attribute__ ((visibility ("default")))
void print_gtf_data(GTF_DATA *gtf_data) {
	int i;
	for (i = 0; i < gtf_data->size; i++) print_row(stdout, gtf_data->data[i], '\t');
}

/*__attribute__ ((visibility ("default")))
TAB_RESULT *extract_data_lib(char *gtf_filename, char *key) {
	char *buffer;
	int i, n;
	size_t size = 0;
	FILE *output = open_memstream(&buffer, &size);
	TAB_RESULT *ret = (TAB_RESULT *)calloc(1, sizeof(TAB_RESULT));

	column = NULL;
	data = NULL;
	nb_column = nb_row = 0;
	init_ftp_data();
	init_gtf(gtf_filename, key, &column, &nb_column, &data, &nb_row);
	n = extract_data(output, gtf_filename, key, column, data, nb_row, nb_column);
	fflush(output);
	fclose(output);

	ret->data = (char **)calloc(n, sizeof(char *));
	FILE *input = fmemopen(buffer, size, "r");
	char *line = (char *)calloc(10000, sizeof(char));
	i = 0;
	while (fgets(line, 9999, input) != NULL) {
		*(line + strlen(line) - 1) = 0;
		ret->data[i] = strdup(line);
		i++;
	}
	ret->size = n;
	fclose(input);
	free(line);
	if (buffer != NULL) free(buffer);
	return ret;
}

__attribute__ ((visibility ("default")))
TAB_RESULT *get_feature_list_lib(char *gtf_filename) {
	char *buffer;
	size_t size = 0;
	FILE *output = open_memstream(&buffer, &size);
	TAB_RESULT *ret = (TAB_RESULT *)calloc(1, sizeof(TAB_RESULT));
	char *line = (char *)calloc(10000, sizeof(char)), *ptr;
	int n, i;

	column = NULL;
	data = NULL;
	nb_column = nb_row = 0;
	init_ftp_data();
	init_gtf(gtf_filename, "feature", &column, &nb_column, &data, &nb_row);
	n = get_feature_list(output);
	fflush(output);
	fclose(output);
	FILE *input = fmemopen(buffer, size, "r");
	ret->data = (char **)calloc(n, sizeof(char *));
	for (i = 0; i < n; i++) {
		fgets(line, 9999, input);
		*(line + strlen(line) - 1) = 0;
		ret->data[i] = strdup(line);
		ptr = strchr(ret->data[i], (int)':');
		if (ptr != NULL) *ptr = '\t';
	}
	ret->size = n;
	if (buffer != NULL) free(buffer);
	return ret;
}

__attribute__ ((visibility ("default")))
TAB_RESULT *get_attribute_list_lib(char *gtf_filename) {
	char *buffer;
	size_t size = 0;
	FILE *output = open_memstream(&buffer, &size);
	TAB_RESULT *ret = (TAB_RESULT *)calloc(1, sizeof(TAB_RESULT));
	char *line = (char *)calloc(10000, sizeof(char)), *ptr;
	int n, i;

	column = NULL;
	data = NULL;
	nb_column = nb_row = 0;
	init_ftp_data();
	init_gtf(gtf_filename, NULL, &column, &nb_column, &data, &nb_row);
	n = get_attribute_list(output);
	fflush(output);
	fclose(output);
	FILE *input = fmemopen(buffer, size, "r");
	ret->data = (char **)calloc(n, sizeof(char *));
	for (i = 0; i < n; i++) {
		fgets(line, 9999, input);
		*(line + strlen(line) - 1) = 0;
		ret->data[i] = strdup(line);
		ptr = strchr(ret->data[i], (int)':');
		if (ptr != NULL) *ptr = '\t';
	}
	ret->size = n;
	if (buffer != NULL) free(buffer);
	return ret;
}

__attribute__ ((visibility ("default")))
TAB_RESULT *get_seq_list_lib(char *gtf_filename) {
	char *buffer;
	size_t size = 0;
	FILE *output = open_memstream(&buffer, &size);
	TAB_RESULT *ret = (TAB_RESULT *)calloc(1, sizeof(TAB_RESULT));
	char *line = (char *)calloc(10000, sizeof(char)), *ptr;
	int n, i;

	column = NULL;
	data = NULL;
	nb_column = nb_row = 0;
	init_ftp_data();
	init_gtf(gtf_filename, "seqid", &column, &nb_column, &data, &nb_row);
	n = get_seq_list(output);
	fflush(output);
	fclose(output);
	FILE *input = fmemopen(buffer, size, "r");
	ret->data = (char **)calloc(n, sizeof(char *));
	for (i = 0; i < n; i++) {
		fgets(line, 9999, input);
		*(line + strlen(line) - 1) = 0;
		ret->data[i] = strdup(line);
		ptr = strchr(ret->data[i], (int)':');
		if (ptr != NULL) *ptr = '\t';
	}
	ret->size = n;
	if (buffer != NULL) free(buffer);
	return ret;
}

__attribute__ ((visibility ("default")))
GTF_RESULT *select_by_genomic_location_lib(char *gtf_filename, char *chr, int start, int end) {
	GTF_RESULT *ret = (GTF_RESULT *)calloc(1, sizeof(GTF_RESULT));
	int i, j, k, k0, n;
	size_t size = 0;
	char *buffer = NULL, **token, **attr;
	char *line = (char *)calloc(10000, sizeof(char));
	FILE *output = open_memstream(&buffer, &size);

	column = NULL;
	data = NULL;
	nb_column = nb_row = 0;
	init_ftp_data();
	init_gtf(gtf_filename, "seqid", &column, &nb_column, &data, &nb_row);

	ret->size = select_by_genomic_location(output, gtf_filename, chr, start, end, column, data, nb_column);

	fflush(output);
	fclose(output);
	ret->data = (GTF_ROW **)calloc(ret->size, sizeof(GTF_ROW *));
	i = 0;
	FILE *input = fmemopen(buffer, size, "r");
	while (fgets(line, 9999, input) != NULL) {
		ret->data[i] = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
		ret->data[i]->field = (char **)calloc(8, sizeof(char *));
		n = split_ip(&token, line, "\t\n");
		for (j = 0; j < 8; j++) ret->data[i]->field[j] = strdup(token[j]);
		if (n == 9) {
			ret->data[i]->attribute = (GTF_ROW_ATTRIBUTES *)calloc(1, sizeof(GTF_ROW_ATTRIBUTES));
			n = split_ip(&attr, token[8], ";\n");
			ret->data[i]->attribute->size = n;
			ret->data[i]->attribute->key = (char **)calloc(n, sizeof(char *));
			ret->data[i]->attribute->value = (char **)calloc(n, sizeof(char *));
			for (j = 0; j < n; j++) {
				k = 0;
				while (*(attr[j] + k) == ' ') k++;
				k0 = k;
				while (*(attr[j] + k) != ' ' && *(attr[j] + k) != '\n') k++;
				if (*(attr[j] + k) != '\n') {
					ret->data[i]->attribute->key[j] = strndup(attr[j] + k0, k - k0);
					if (*(attr[j] + k + 1) == '"') k++;
					if (*(attr[j] + k + 1 + strlen(attr[j] + k + 1) - 1) == '"') *(attr[j] + k + 1 + strlen(attr[j] + k + 1) - 1) = 0;
					ret->data[i]->attribute->value[j] = strdup(attr[j] + k + 1);
				}
			}
			free(attr);
		}
		free(token);
		i++;
	}
	fclose(input);

	if (buffer != NULL) free(buffer);
	free(line);
	return ret;
}

__attribute__ ((visibility ("default")))
FASTA_RESULT *get_fasta_lib(char *gtf_filename, char *genome_file, int intron, int rc) {
	char *buffer = NULL;
	size_t size = 0;
	FILE *output = open_memstream(&buffer, &size);
	FASTA_RESULT *ret = (FASTA_RESULT *)calloc(1, sizeof(FASTA_RESULT));

	column = NULL;
	data = NULL;
	nb_column = nb_row = 0;
	init_ftp_data();
	init_gtf(gtf_filename, NULL, &column, &nb_column, &data, &nb_row);

	int n = get_fasta(output, gtf_filename, genome_file, intron, rc, column, data, nb_row);
	fflush(output);
	fclose(output);

	FILE *input = fmemopen(buffer, size, "r");
	ret->size = n;
	ret->data = (FASTA **)calloc(n, sizeof(FASTA *));
	int i;
	for (i = 0; i < n; i++) {
		ret->data[i] = (FASTA *)calloc(1, sizeof(FASTA));
		size = 0;
		getline(&(ret->data[i]->header), &size, input);
		if (*(ret->data[i]->header + strlen(ret->data[i]->header) - 1) == '\n') *(ret->data[i]->header + strlen(ret->data[i]->header) - 1) = 0;
		getline(&(ret->data[i]->sequence), &size, input);
		if (*(ret->data[i]->sequence + strlen(ret->data[i]->sequence) - 1) == '\n') *(ret->data[i]->sequence + strlen(ret->data[i]->sequence) - 1) = 0;
	}
	fclose(input);
	if (buffer != NULL) free(buffer);
	return ret;
}*/
