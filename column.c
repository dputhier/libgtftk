/*
 * common.c
 *
 *  Created on: Jan 9, 2017
 *      Author: fafa
 */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "libgtftk.h"

extern int split_ip(char ***tab, char *s, char *delim);
extern char *trim_ip(char *s);
extern int compare_row_list(const void *p1, const void *p2);
extern int nb_column;
extern COLUMN **column;

void *convert_string(char *token, void *def, GTF_ROW *row) {
	char *s;

	if (!strcmp(token, "."))
		s = strdup(def);
	else
		s = strdup(token);
	return s;
}

char *convert_from_string(void *f, void *def) {
	char *ret = (char *)calloc(100, sizeof(char));

	if (!strcmp((char *)f, (char *)def))
			ret = strdup(".");
		else {
			snprintf(ret, 100, "%s", (char *)f);
			ret = realloc(ret, strlen(ret) + 1);
		}
	return ret;
}

void *convert_int(char *token, void *def, GTF_ROW *row) {
	int *i = (int *)calloc(1, sizeof(int));

	if (!strcmp(token, "."))
		*i = *(int *)def;
	else
		*i = atoi(token);
	return i;
}

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

void *convert_float(char *token, void *def, GTF_ROW *row) {
	float *f = (float *)calloc(1, sizeof(float));

	if (!strcmp(token, "."))
		*f = *(float *)def;
	else
		*f = atof(token);
	return f;
}

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

void *convert_char(char *token, void *def, GTF_ROW *row) {
	char *c = (char *)calloc(1, sizeof(char));

	if (!strcmp(token, "."))
		*c = *(char *)def;
	else
		*c = *token;
	return c;
}

char *convert_from_char(void *f, void *def) {
	char *ret = (char *)calloc(100, sizeof(char));

	if (*(char *)f == *(char *)def)
			ret = strdup(".");
		else {
			snprintf(ret, 100, "%c", *(char *)f);
			ret = realloc(ret, strlen(ret) + 1);
		}
	return ret;
}

void *convert_attributes(char *token, void *def, GTF_ROW *row) {
	char **attribute, **key_value;
	int j, na;
	ATTRIBUTES *attributes;

	while (*(token + strlen(token) - 1) == ' ') *(token + strlen(token) - 1) = 0;
	na = split_ip(&attribute, token, ";");

	row->nb_attributes = na;
	attributes = (ATTRIBUTES *)calloc(1, sizeof(ATTRIBUTES));
	attributes->attr = (ATTRIBUTE **)calloc(na, sizeof(ATTRIBUTE *));
	attributes->nb = na;
	for (j = 0; j < na; j++) {
		split_ip(&key_value, attribute[j], "\"");
		attributes->attr[j] = (ATTRIBUTE *)calloc(1, sizeof(ATTRIBUTE));
		attributes->attr[j]->key = strdup(trim_ip(key_value[0]));
		attributes->attr[j]->value = strdup(trim_ip(key_value[1]));
		free(key_value);
	}
	free(attribute);
	return attributes;
}

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

void print_int(void *token, FILE *output, void *col, char delim) {
	if (((COLUMN *)col)->default_value != NULL)
		if (*(int *)(token) == *(int *)((COLUMN *)col)->default_value)
			delim != 0 ? fprintf(output, ".%c", delim) : fprintf(output, ".");
		else
			delim != 0 ? fprintf(output, "%d%c", *(int *)(token), delim) : fprintf(output, "%d", *(int *)(token));
	else
		delim != 0 ? fprintf(output, "%d%c", *(int *)(token), delim) : fprintf(output, "%d", *(int *)(token));
}

void print_float(void *token, FILE *output, void *col, char delim) {
	if (((COLUMN *)col)->default_value != NULL)
		if (*(float *)(token) == *(float *)((COLUMN *)col)->default_value)
			delim != 0 ? fprintf(output, ".%c", delim) : fprintf(output, ".");
		else
			delim != 0 ? fprintf(output, "%f%c", *(float *)(token), delim) : fprintf(output, "%f", *(float *)(token));
	else
		delim != 0 ? fprintf(output, "%f%c", *(float *)(token), delim) : fprintf(output, "%f", *(float *)(token));
}

void print_string(void *token, FILE *output, void *col, char delim) {
	if (((COLUMN *)col)->default_value != NULL)
		if (!strcmp((char *)(token), (char *)((COLUMN *)col)->default_value))
			delim != 0 ? fprintf(output, ".%c", delim) : fprintf(output, ".");
		else
			delim != 0 ? fprintf(output, "%s%c", (char *)(token), delim) : fprintf(output, "%s", (char *)(token));
	else
		delim != 0 ? fprintf(output, "%s%c", (char *)(token), delim) : fprintf(output, "%s", (char *)(token));
}

void print_char(void *token, FILE *output, void *col, char delim) {
	if (((COLUMN *)col)->default_value != NULL)
		if (*(char *)(token) == *(char *)((COLUMN *)col)->default_value)
			delim != 0 ? fprintf(output, ".%c", delim) : fprintf(output, ".");
		else
			delim != 0 ? fprintf(output, "%c%c", *(char *)(token), delim) : fprintf(output, "%c", *(char *)(token));
	else
		delim != 0 ? fprintf(output, "%c%c", *(char *)(token), delim) : fprintf(output, "%c", *(char *)(token));
}

void print_attribute(void *token, char *attr, FILE *output, void *r, char delim) {
	int k;
	if (((GTF_ROW *)r)->nb_attributes > 0) {
		for (k = 0; k < ((ATTRIBUTES *)token)->nb; k++) {
			if (!strcmp(attr, ((ATTRIBUTES *)token)->attr[k]->key)) {
				delim != 0 ? fprintf(output, "%s%c", ((ATTRIBUTES *)token)->attr[k]->value, delim) : fprintf(output, "%s", ((ATTRIBUTES *)token)->attr[k]->value);
				break;
			}
		}
		if (k == ((ATTRIBUTES *)token)->nb) delim != 0 ? fprintf(output, "NA%c", delim) : fprintf(output, "NA");
	}
}
void print_attributes(void *token, FILE *output, void *r, char delim) {
	int k;
	if (((GTF_ROW *)r)->nb_attributes > 0) {
		fprintf(output, "%s \"%s\";", ((ATTRIBUTES *)token)->attr[0]->key, ((ATTRIBUTES *)token)->attr[0]->value);
		for (k = 1; k < ((ATTRIBUTES *)token)->nb; k++)
			fprintf(output, " %s \"%s\";", ((ATTRIBUTES *)token)->attr[k]->key, ((ATTRIBUTES *)token)->attr[k]->value);
	}
}

void print_row(FILE *output, GTF_ROW *r, char delim) {
	int i;

	for (i = 0; i < nb_column - 1; i++)
		if (column[i]->type == 'A')
			column[i]->print(r->data[i], output, r, delim);
		else
			column[i]->print(r->data[i], output, column[i], delim);
	if (column[i]->type == 'A')
		column[i]->print(r->data[i], output, r, 0);
	else
		column[i]->print(r->data[i], output, column[i], 0);
	fprintf(output, "\n");
}

void make_index_string(int i, GTF_ROW *row, COLUMN *col) {
	ROW_LIST *test_row_list = calloc(1, sizeof(ROW_LIST)), *row_list, **find_row_list;
	test_row_list->token = row->data[col->num];
	find_row_list = tfind(test_row_list, &(col->index[0]->data), compare_row_list);
	if (find_row_list == NULL) {
		row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
		row_list->token = row->data[col->num];
		row_list->nb_row = 1;
		row_list->row = (int *)calloc(1, sizeof(int));
		row_list->row[row_list->nb_row - 1] = i;
		tsearch(row_list, &(col->index[0]->data), compare_row_list);
	}
	else {
		row_list = *((ROW_LIST **)find_row_list);
		row_list->nb_row++;
		row_list->row = (int *)realloc(row_list->row, row_list->nb_row * sizeof(int));
		row_list->row[row_list->nb_row - 1] = i;
	}
	free(test_row_list);
}

void make_index_attribute(int i, GTF_ROW *row, COLUMN *col) {
	ROW_LIST *test_row_list = calloc(1, sizeof(ROW_LIST)), *row_list, **find_row_list;
	int k, j;

	for (k = 0; k < col->nb_index; k++)
		for (j = 0; j < row->nb_attributes; j++)
			if (!strcmp(col->index[k]->key, ((ATTRIBUTES *)row->data[col->num])->attr[j]->key)) {
				test_row_list->token = ((ATTRIBUTES *)row->data[col->num])->attr[j]->value;
				find_row_list = tfind(test_row_list, &(col->index[k]->data), compare_row_list);
				if (find_row_list == NULL) {
					row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
					row_list->token = ((ATTRIBUTES *)row->data[col->num])->attr[j]->value;
					row_list->nb_row = 1;
					row_list->row = (int *)calloc(1, sizeof(int));
					row_list->row[row_list->nb_row - 1] = i;
					tsearch(row_list, &(col->index[k]->data), compare_row_list);
				}
				else {
					row_list = *((ROW_LIST **)find_row_list);
					row_list->nb_row++;
					row_list->row = (int *)realloc(row_list->row, row_list->nb_row * sizeof(int));
					row_list->row[row_list->nb_row - 1] = i;
				}
				break;
			}
}

void make_index_int(int i, GTF_ROW *row, COLUMN *col) {

}

void make_index_float(int i, GTF_ROW *row, COLUMN *col) {

}

void make_index_char(int i, GTF_ROW *row, COLUMN *col) {
	ROW_LIST *test_row_list = calloc(1, sizeof(ROW_LIST)), *row_list, **find_row_list;
	test_row_list->token = row->data[col->num];
	find_row_list = tfind(test_row_list, &(col->index[0]->data), compare_row_list);
	if (find_row_list == NULL) {
		row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
		row_list->token = row->data[col->num];
		row_list->nb_row = 1;
		row_list->row = (int *)calloc(1, sizeof(int));
		row_list->row[row_list->nb_row - 1] = i;
		tsearch(row_list, &(col->index[0]->data), compare_row_list);
	}
	else {
		row_list = *((ROW_LIST **)find_row_list);
		row_list->nb_row++;
		row_list->row = (int *)realloc(row_list->row, row_list->nb_row * sizeof(int));
		row_list->row[row_list->nb_row - 1] = i;
	}
	free(test_row_list);
}

COLUMN *make_column(char type, int i, void *dv, char *name) {
	COLUMN *column = (COLUMN *)calloc(1, sizeof(COLUMN));
	column->type = type;
	column->num = i;
	column->name = strdup(name);
	column->index = NULL;
	column->nb_index = 0;
	if (type == 'S') {
		if (dv != NULL) column->default_value = strdup((char *)dv);
		column->convert = convert_string;
		column->convert_to_string = convert_from_string;
		column->make_index = make_index_string;
		column->print = print_string;
	}
	else if (type == 'I') {
		column->default_value = (int *)calloc(1, sizeof(int));
		if (dv != NULL) *((int *)(column->default_value)) = *(int *)dv;
		column->convert = convert_int;
		column->convert_to_string = convert_from_int;
		column->make_index = make_index_int;
		column->print = print_int;
	}
	else if (type == 'F') {
		column->default_value = (float *)calloc(1, sizeof(float));
		if (dv != NULL) (*((float *)(column->default_value))) = *(float *)dv;
		column->convert = convert_float;
		column->convert_to_string = convert_from_float;
		column->make_index = make_index_float;
		column->print = print_float;
	}
	else if (type == 'C') {
		column->default_value = (char *)malloc(sizeof(char));
		if (dv != NULL) (*((char *)(column->default_value))) = *(char *)dv;
		column->convert = convert_char;
		column->convert_to_string = convert_from_char;
		column->make_index = make_index_char;
		column->print = print_char;
	}
	else if (type == 'A') {
		if (dv != NULL) column->default_value = strdup((char *)dv);
		column->convert = convert_attributes;
		column->convert_to_string = convert_from_attributes;
		column->make_index = make_index_attribute;
		column->print = print_attributes;
	}
	return column;
}

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

	/*column[8]->nb_index = 3;
	column[8]->index = calloc(3, sizeof(INDEX *));
	for (i = 0; i < 3; i++) column[8]->index[i] = (INDEX *)calloc(1, sizeof(INDEX));
	column[8]->index[0]->key = strdup("gene_id");
	column[8]->index[1]->key = strdup("transcript_id");
	column[8]->index[1]->key = strdup("ccds_id");*/

	for (i = 0; i < (nb_column - 1); i++) {
		column[i]->index = (INDEX **)calloc(1, sizeof(INDEX *));
		column[i]->index[0] = (INDEX *)calloc(1, sizeof(INDEX));
		column[i]->index[0]->key = column[i]->name;
		column[i]->nb_index = 1;
	}
}
