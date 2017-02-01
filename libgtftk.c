/*
 * gtftk_lib.c
 *
 *  Created on: Apr 1, 2016
 *      Author: fafa
 *
 *  The aim of this library is to allow clients written in other languages
 *  than C to take advantage of the speed of C language to access and
 *  manipulate data stored in GTF files. Thus, libgtftk is a library that
 *  provides a set of basic functions to parse, index and extract data from a
 *  GTF file. Those functions can be called from a client written in any
 *  language (C, C++, Java, Python ...) provided that a dedicated C interface
 *  is available (ex : JNI for Java or ctypes for Python). Obviously, C and
 *  C++ clients doesn't need such an interface. Each GTF related function is
 *  implemented in a separate source file. So, this file is just here to hold
 *  common utility functions and most of them are not callable from the shared
 *  library.
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

/*
 * This is a pointer on the column model of a GTF file. It is initialized by
 * the make_columns() function in column.c source file. It is accessed by most
 * of the functions of the library.
 */
COLUMN **column;

/*
 * The number of columns, that is set to 9 for GTF files in column.c source
 * file.
 */
int nb_column;

/*
 * A variable used to count the elements in a index with the twalk C function.
 * Used in action_nb function in this file.
 */
int N;

/*
 * This function splits a character string (s) into a table of words (*tab),
 * according to a set of delimiters (delim). Each character of delim is a
 * delimiter. The string is splitted in place and the resulting word table
 * just contains pointers to the words.
 * Parameters:
 * 		tab:	the address of a pointer on a table of string characters
 * 				the table must NOT be reserved before calling this function
 * 		s:		the character string to be splitted
 * 		delim:	the set of charater delimiters
 *
 * Return:		the number of words that were found
 */
int split_ip(char ***tab, char *s, char *delim) {
	int i, n, k, in_token, l;

	in_token = n = k = 0;
	l = strlen(s);
	for (i = 0; i < l; i++)
		if (strchr(delim, (int)(*(s + i))) != NULL) {
			*(s + i) = 0;
			in_token = 0;
		}
		else if (!in_token) {
			in_token = 1;
			n++;
		}
	*tab = (char **)calloc(n, sizeof(char *));
	for (i = 0; i < l; i++)
		if (*(s + i ) != 0) {
			(*tab)[k++] = s + i;
			i += strlen(s + i);
		}
	return n;
}

/*
 * This function trim the extra space characters at the beginning and at the
 * end of the string s. The trimming is performed in place, so no extra memory
 * is reserved to store the trimmed string. The returned value is a pointer on
 * the first non space character of s, or NULL if s contains only space
 * characters.
 *
 * Parameters:
 * 		s:		the string to be trimmed
 *
 * Return:		a pointer on the trimmed string
 */
char *trim_ip(char *s) {
	int b, e, l;

	l = strlen(s);
	for (b = 0; b < l; b++)
		if (*(s + b) != ' ')
			break;
	for (e = l - 1; e > 0; e--)
		if (*(s + e) == ' ')
			*(s + e) = 0;
		else
			break;
	return s + b;
}

/*
 * This function split a key/value pair (as specified in the GTF format for
 * the attributes column) into 2 character strings. value can be delimitted by
 * double quote characters like in Ensembl format or not. The strings pointed
 * by key and value are allocated in this function.
 *
 * Parameters:
 * 		s:		the key/value pair from a GTF attribute
 * 		key:	a pointer to store the address of the key
 * 		value:	a pointer to store the address of the value
 */
void split_key_value(char *s, char **key, char **value) {
	int k = 0;

	while (*s == ' ') s++;
	while (*(s + k) != ' ') k++;
	*(s + k) = 0;
	*key = strdup(s);
	s += k + 1;
	while ((*s == ' ') || (*s == '"')) s++;
	k = 0;
	while ((*(s + k) != '"') && (*(s + k) != ' ') && (*(s + k) != 0)) k++;
	*(s + k) = 0;
	*value = strdup(s);
}

/*
 * This function is used by the C tsearch and tfind functions to search and
 * add ROW_LIST elements in all the indexes on the columns of a GTF file. For
 * more information, see the man pages of tsearch.
 */
int compare_row_list(const void *p1, const void *p2) {
	ROW_LIST *rl1 = ((ROW_LIST *)p1);
	ROW_LIST *rl2 = ((ROW_LIST *)p2);
	return strcmp(rl1->token, rl2->token);
}

/*
 * This function is used by the C qsort function to sort a table of integers
 * representing rows in a GTF file. As some operations can modify the order of
 * the rows in the results, this function is used to be sure to output the
 * rows in the original order.
 */
int comprow(const void *m1, const void *m2) {
	 int *r1 = (int *)m1;
	 int *r2 = (int *)m2;
	 return *r1 - *r2;
}

/*
 * These functions merge two lists of rows (ROW_LIST structures) with
 * suppression of duplications. The resulting list is not sorted.
 *
 * Parameters:
 * 		scr: 	pointer on the source list
 * 		dsr:	pointer on the destination list
 *
 * Return:		the number of rows in dest list, after merging
 */
int add_row(int src, ROW_LIST *dst) {
	if (dst == NULL) {
		dst = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
		dst->row = (int *)calloc(1, sizeof(int));
		dst->row[dst->nb_row] = src;
		dst->nb_row++;
	}
	else if (bsearch(&src, dst->row, dst->nb_row, sizeof(int), comprow) == NULL) {
		dst->row = (int *)realloc(dst->row, (dst->nb_row + 1) * sizeof(int));
		dst->row[dst->nb_row] = src;
		dst->nb_row++;
	}
	return dst->nb_row;
}

int add_row_list(ROW_LIST *src, ROW_LIST *dst) {
	int i;
	for (i = 0; i < src->nb_row; i++) add_row(src->row[i], dst);
	return dst->nb_row;
}

/*
 * This function prints the content of a GTF_DATA on the standard output. It
 * can be called from the library (default visibility).
 *
 * Parameters:
 * 		gtf_data:	a pointer on the GTF data to be printed
 */
__attribute__ ((visibility ("default")))
void print_gtf_data(GTF_DATA *gtf_data) {
	int i;
	for (i = 0; i < gtf_data->size; i++) print_row(stdout, gtf_data->data[i], '\t');
}

__attribute__ ((visibility ("default")))
GTF_ROW_CHAR *gtf_row_to_char(GTF_ROW *row) {
	int i;
	GTF_ROW_CHAR *ret = (GTF_ROW_CHAR *)calloc(1, sizeof(GTF_ROW_CHAR));
	ret->rank = row->rank;
	ret->data = (char **)calloc(9, sizeof(char *));
	for (i = 0; i < 9; i++)
		ret->data[i] = column[i]->convert_to_string(row->data[i], column[i]->default_value);
	return ret;
}

/*
 * This function is intended to be used with the C twalk function. It is used
 * to evaluate the number of elements in an index. the N variable is declared
 * at the beginning of this file. For more information about twalk, see the
 * man pages.
 */
static void action_nb(const void *nodep, const VISIT which, const int depth) {
	switch (which) {
		case preorder:
			break;
		case postorder:
		case leaf:
			N++;
			break;
		case endorder:
			break;
	}
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
