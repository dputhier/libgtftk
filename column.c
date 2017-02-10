/*
 * column.c
 *
 *  Created on: Jan 9, 2017
 *      Author: fafa
 *
 *  This source file contains all that is related to column management :
 *  	- all COLUMN structure functions, for each kind of column (string, int,
 *  	  float, char, attributes): convert, convert_to_string, index_row and
 *  	  print
 *		- print_row to print a whole row with a given delimiter
 *		- index_row to add a row in an index corresponding to the given value
 *		- make_columns/make_column to build the GTF column model
 */

#include "libgtftk.h"

/*
 * external functions declaration
 */
extern int split_ip(char ***tab, char *s, char *delim);
extern void split_key_value(char *s, char **key, char **value);
extern int compare_row_list(const void *p1, const void *p2);

/*
 * global variables declaration
 */
extern int nb_column;
extern COLUMN **column;

/*
 * Conversion functions between string and column data types. Default values
 * are taken into account. These functions are the "convert" function of the
 * COLUMN structure
 *
 * Parameters:
 * 		token:	the value to convert
 * 		def:	the default value in the appropriate type
 *
 * Returns:		a void pointer to the converted value
 */

// String to string conversion
void *convert_string(char *token, void *def) {
	char *s;

	if (!strcmp(token, "."))
		s = strdup(def);
	else
		s = strdup(token);
	return s;
}

// String to int conversion
void *convert_int(char *token, void *def) {
	int *i = (int *)calloc(1, sizeof(int));

	if (!strcmp(token, "."))
		*i = *(int *)def;
	else
		*i = atoi(token);
	return i;
}

// String to float conversion
void *convert_float(char *token, void *def) {
	float *f = (float *)calloc(1, sizeof(float));

	if (!strcmp(token, "."))
		*f = *(float *)def;
	else
		*f = atof(token);
	return f;
}

// String to char conversion
void *convert_char(char *token, void *def) {
	char *c = (char *)calloc(1, sizeof(char));

	if (!strcmp(token, "."))
		*c = *(char *)def;
	else
		*c = *token;
	return c;
}

// String to ATTRIBUTES conversion
void *convert_attributes(char *token, void *def) {
	char **attribute;
	int j, na, l;
	ATTRIBUTES *attributes;

	// remove extra spaces at the end of the attributes list
	l = strlen(token);
	while (*(token + l - 1) == ' ') l--;
	*(token + l) = 0;

	// split attributes (separated by ";" characters)
	na = split_ip(&attribute, token, ";");

	// reserve and initialize the ATTRIBUTES structure
	attributes = (ATTRIBUTES *)calloc(1, sizeof(ATTRIBUTES));
	attributes->attr = (ATTRIBUTE **)calloc(na, sizeof(ATTRIBUTE *));
	attributes->nb = na;

	/* For each attribute, reserve the ATTRIBUTE structure and call the
	 * split_key_value function to fill it
	 */
	for (j = 0; j < na; j++) {
		attributes->attr[j] = (ATTRIBUTE *)calloc(1, sizeof(ATTRIBUTE));
		split_key_value(attribute[j], &(attributes->attr[j]->key), &(attributes->attr[j]->value));
	}

	/* free the splitted attributes (they have been duplicated by
	 * split_key_value function) and return
	 */
	free(attribute);
	return attributes;
}

/*
 * Conversion functions between column data types and string (the opposite of
 * the previous functions). Default values are taken into account. These
 * functions are the "convert_to_string" function of the COLUMN structure
 *
 * Parameters:
 * 		s, i, f, c, attributes:		a pointer to the typed value
 * 		def:						a pointer to the typed default value
 *
 * Returns:							the string value
 */

// String to string conversion
char *convert_from_string(void *s, void *def) {
	char *ret = (char *)calloc(100, sizeof(char));

	if (!strcmp((char *)s, (char *)def))
			ret = strdup(".");
		else {
			snprintf(ret, 100, "%s", (char *)s);
			ret = realloc(ret, strlen(ret) + 1);
		}
	return ret;
}

// int to string conversion
char *convert_from_int(void *i, void *def) {
	char *ret = (char *)calloc(100, sizeof(char));

	if (*(int *)i == *(int *)def)
		ret = strdup(".");
	else {
		snprintf(ret, 100, "%d", *(int *)i);
		ret = realloc(ret, strlen(ret) + 1);
	}
	return ret;
}

// float to string conversion
char *convert_from_float(void *f, void *def) {
	char *ret = (char *)calloc(100, sizeof(char));

	if (*(float *)f == *(float *)def)
		ret = strdup(".");
	else {
		snprintf(ret, 100, "%f", *(float *)f);
		ret = realloc(ret, strlen(ret) + 1);
	}
	return ret;
}

// char to string conversion
char *convert_from_char(void *c, void *def) {
	char *ret = (char *)calloc(100, sizeof(char));

	if (*(char *)c == *(char *)def)
			ret = strdup(".");
		else {
			snprintf(ret, 100, "%c", *(char *)c);
			ret = realloc(ret, strlen(ret) + 1);
		}
	return ret;
}

/*
 * ATTRIBUTES to string conversion
 * prints each ATTRIBUTE as: key "value";
 * double quote characters are added around value, even if they were not
 * present in the GTF file
 */
char *convert_from_attributes(void *attributes, void *def) {
	char *ret = (char *)calloc(1000, sizeof(char));
	int j, na;

	na = ((ATTRIBUTES *)attributes)->nb;
	for (j = 0; j < na; j++) {
		strcat(ret, ((ATTRIBUTES *)attributes)->attr[j]->key);
		strcat(ret, " \"");
		strcat(ret, ((ATTRIBUTES *)attributes)->attr[j]->value);
		strcat(ret, "\"; ");
	}
	*(ret + strlen(ret) - 1) = 0;
	return ret;
}

/*
 * These functions prints in output a column value (a feature, a start, an
 * attributes list ...) with delimiter. Default values are printed as ".".
 * These functions are the "print" function of the COLUMN structure.
 *
 * Parameters:
 * 		token:	the typed value to print
 * 		output:	where to print (a file or stdout)
 * 		col:	the corresponding column (NULL if attributes column)
 * 		delim:	the delimiter character
 */

// print for integer columns (start, end and phase)
void print_int(void *token, FILE *output, void *col, char delim) {
	if (((COLUMN *)col)->default_value != NULL)
		if (*(int *)(token) == *(int *)((COLUMN *)col)->default_value)
			delim != 0 ? fprintf(output, ".%c", delim) : fprintf(output, ".");
		else
			delim != 0 ? fprintf(output, "%d%c", *(int *)(token), delim) : fprintf(output, "%d", *(int *)(token));
	else
		delim != 0 ? fprintf(output, "%d%c", *(int *)(token), delim) : fprintf(output, "%d", *(int *)(token));
}

// print for float columns (score)
void print_float(void *token, FILE *output, void *col, char delim) {
	if (((COLUMN *)col)->default_value != NULL)
		if (*(float *)(token) == *(float *)((COLUMN *)col)->default_value)
			delim != 0 ? fprintf(output, ".%c", delim) : fprintf(output, ".");
		else
			delim != 0 ? fprintf(output, "%f%c", *(float *)(token), delim) : fprintf(output, "%f", *(float *)(token));
	else
		delim != 0 ? fprintf(output, "%f%c", *(float *)(token), delim) : fprintf(output, "%f", *(float *)(token));
}

// print for string columns (seqid, source and feature)
void print_string(void *token, FILE *output, void *col, char delim) {
	if (((COLUMN *)col)->default_value != NULL)
		if (!strcmp((char *)(token), (char *)((COLUMN *)col)->default_value))
			delim != 0 ? fprintf(output, ".%c", delim) : fprintf(output, ".");
		else
			delim != 0 ? fprintf(output, "%s%c", (char *)(token), delim) : fprintf(output, "%s", (char *)(token));
	else
		delim != 0 ? fprintf(output, "%s%c", (char *)(token), delim) : fprintf(output, "%s", (char *)(token));
}

// print for char columns (strand)
void print_char(void *token, FILE *output, void *col, char delim) {
	if (((COLUMN *)col)->default_value != NULL)
		if (*(char *)(token) == *(char *)((COLUMN *)col)->default_value)
			delim != 0 ? fprintf(output, ".%c", delim) : fprintf(output, ".");
		else
			delim != 0 ? fprintf(output, "%c%c", *(char *)(token), delim) : fprintf(output, "%c", *(char *)(token));
	else
		delim != 0 ? fprintf(output, "%c%c", *(char *)(token), delim) : fprintf(output, "%c", *(char *)(token));
}

// print for attributes column
void print_attributes(void *token, FILE *output, void *col, char delim) {
	int k;
	if (((ATTRIBUTES *)token)->nb > 0) {
		fprintf(output, "%s \"%s\";", ((ATTRIBUTES *)token)->attr[0]->key, ((ATTRIBUTES *)token)->attr[0]->value);
		for (k = 1; k < ((ATTRIBUTES *)token)->nb; k++)
			fprintf(output, " %s \"%s\";", ((ATTRIBUTES *)token)->attr[k]->key, ((ATTRIBUTES *)token)->attr[k]->value);
	}
}

/*
 * This function adds a row (row_nb) into an index.
 *
 * Parameters:
 * 		row_nb:		the row number to add
 * 		value:		the value with which to row is associated
 * 		index:		the index of ROW_LIST elements
 */
void index_row(int row_nb, char *value, INDEX *index) {
	ROW_LIST *test_row_list, *row_list, *find_row_list;

	if (index != NULL) {
		// build a ROW_LIST to check if value is already indexed
		test_row_list = calloc(1, sizeof(ROW_LIST));
		test_row_list->token = value;
		find_row_list = tfind(test_row_list, &(index->data), compare_row_list);

		if (find_row_list == NULL) {
			/*
			 * value is not in the index so we reserve a ROW_LIST, initialize it
			 * with one row (row_nb) and put it in the index (tsearch)
			 */
			row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
			row_list->token = value;
			row_list->nb_row = 1;
			row_list->row = (int *)calloc(1, sizeof(int));
			row_list->row[row_list->nb_row - 1] = row_nb;
			tsearch(row_list, &(index->data), compare_row_list);
		}
		else {
			/*
			 * value is already in the index so we just have to add row_nb into
			 * the ROW_LIST element found
			 */
			row_list = *((ROW_LIST **)find_row_list);
			row_list->nb_row++;
			row_list->row = (int *)realloc(row_list->row, row_list->nb_row * sizeof(int));
			row_list->row[row_list->nb_row - 1] = row_nb;
		}
		free(test_row_list);
	}
}

/*
 * Prints a GTF row in output with the given character delimiter
 *
 * Parameters:
 * 		output:		where to print
 * 		r:			the row to print
 * 		delim:		the delimiter character
 */
void print_row(FILE *output, GTF_ROW *r, char delim) {
	int i;

	for (i = 0; i < nb_column - 1; i++) column[i]->print(r->data[i], output, column[i], delim);
	column[i]->print(r->data[i], output, NULL, 0);
	fprintf(output, "\n");
}

/*
 * Prints the value of an attribute (attr) from token (ATTRIBUTES) in output,
 * whith the given delimiter. If the attribute attr is not in token, prints
 * "NA". This function is used in get_fasta.c to print the header of fasta
 * sequences.
 *
 * Parameters:
 * 		token:		the attributes in which to search
 * 		attr:		the attribute to print (value)
 * 		output:		where to print
 * 		delim:		the delimiter character
 */
/*void print_attribute(void *token, char *attr, FILE *output, char delim) {
	int k;
	if (((ATTRIBUTES *)token)->nb > 0) {
		for (k = 0; k < ((ATTRIBUTES *)token)->nb; k++) {
			if (!strcmp(attr, ((ATTRIBUTES *)token)->attr[k]->key)) {
				delim != 0 ? fprintf(output, "%s%c", ((ATTRIBUTES *)token)->attr[k]->value, delim) : fprintf(output, "%s", ((ATTRIBUTES *)token)->attr[k]->value);
				break;
			}
		}
		if (k == ((ATTRIBUTES *)token)->nb) delim != 0 ? fprintf(output, "NA%c", delim) : fprintf(output, "NA");
	}
}*/
void print_attribute(void *token, char *attr, char *output, char delim) {
	int k;
	if (((ATTRIBUTES *)token)->nb > 0) {
		for (k = 0; k < ((ATTRIBUTES *)token)->nb; k++) {
			if (!strcmp(attr, ((ATTRIBUTES *)token)->attr[k]->key)) {
				delim != 0 ? fprintf(output, "%s%c", ((ATTRIBUTES *)token)->attr[k]->value, delim) : fprintf(output, "%s", ((ATTRIBUTES *)token)->attr[k]->value);
				break;
			}
		}
		if (k == ((ATTRIBUTES *)token)->nb) delim != 0 ? fprintf(output, "NA%c", delim) : fprintf(output, "NA");
	}
}

/*
 * Reserve memory for a COLUMN structure and fill it whith parameters.
 * Initialize the function pointers to point on the good functions, depending
 * on the columnn type.
 *
 * Parameters:
 * 		type:	the column type ('S' for string, 'I' for integer, 'F' for
 * 				float, 'C' for char and 'A' for attributes)
 * 		i:		the number of the column (0 to 8)
 * 		dv:		the typed default value
 * 		name:	the column name
 *
 * Returns:		the initialized COLUMN structure
 */
COLUMN *make_column(char type, int i, void *dv, char *name) {
	COLUMN *column = (COLUMN *)calloc(1, sizeof(COLUMN));
	column->type = type;
	column->num = i;
	column->name = strdup(name);
	column->index = NULL;
	column->nb_index = 0;
	column->index_row = index_row;
	if (type == 'S') {
		if (dv != NULL) column->default_value = strdup((char *)dv);
		column->convert = convert_string;
		column->convert_to_string = convert_from_string;
		column->print = print_string;
	}
	else if (type == 'I') {
		column->default_value = (int *)calloc(1, sizeof(int));
		if (dv != NULL) *((int *)(column->default_value)) = *(int *)dv;
		column->convert = convert_int;
		column->convert_to_string = convert_from_int;
		column->print = print_int;
	}
	else if (type == 'F') {
		column->default_value = (float *)calloc(1, sizeof(float));
		if (dv != NULL) (*((float *)(column->default_value))) = *(float *)dv;
		column->convert = convert_float;
		column->convert_to_string = convert_from_float;
		column->print = print_float;
	}
	else if (type == 'C') {
		column->default_value = (char *)malloc(sizeof(char));
		if (dv != NULL) (*((char *)(column->default_value))) = *(char *)dv;
		column->convert = convert_char;
		column->convert_to_string = convert_from_char;
		column->print = print_char;
	}
	else if (type == 'A') {
		if (dv != NULL) column->default_value = strdup((char *)dv);
		column->convert = convert_attributes;
		column->convert_to_string = convert_from_attributes;
		column->print = print_attributes;
	}
	return column;
}

/*
 * Creates and initialize a column model suitable for GTF data.
 * This function creates an index for each column except the last one
 * (attributes). By default, the 8 first columns are indexed but not the last
 * one (it takes too much time). The indexing of the attributes column depends
 * on the kind of operation to perform on GTF data. For example,
 * "select_by_key gene_biotype:protein_coding" needs to index only
 * "gene_biotype" attribute. This function is called by loadGTF function, the
 * first function a client should call to use the library.
 */
void make_columns() {
	int default_value_int = -1;
	float default_value_float = -1;
	int i;

	nb_column = 9;
	column = (COLUMN **)calloc(nb_column, sizeof(COLUMN *));
	column[0] = make_column('S', 0, ".", "seqid");
	column[1] = make_column('S', 1, ".", "source");
	column[2] = make_column('S', 2, ".", "feature");
	column[3] = make_column('I', 3, &default_value_int, "start");
	column[4] = make_column('I', 4, &default_value_int, "end");
	column[5] = make_column('F', 5, &default_value_float, "score");
	column[6] = make_column('C', 6, ".", "strand");
	column[7] = make_column('I', 7, &default_value_int, "phase");
	column[8] = make_column('A', 8, ".", "attributes");

	for (i = 0; i < (nb_column - 1); i++) {
		column[i]->index = (INDEX **)calloc(1, sizeof(INDEX *));
		column[i]->index[0] = (INDEX *)calloc(1, sizeof(INDEX));
		column[i]->index[0]->key = column[i]->name;
		column[i]->nb_index = 1;
	}
}
