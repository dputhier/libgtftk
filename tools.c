/*
 * tools.c
 *
 *  Created on: Jan 10, 2017
 *      Author: fafa
 */

#include "libgtftk.h"

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

void split_key_value(char *s, char **key, char **value) {
	int k = 0;
	while (*s == ' ') s++;
	while (*(s + k) != ' ') k++;
	*(s + k) = 0;
	*key = strdup(s);
	s += k + 1;
	while (*s != '"') s++;
	s++;
	k = 0;
	while (*(s + k) != '"') k++;
	*(s + k) = 0;
	*value = strdup(s);
}

int compare_row_list(const void *p1, const void *p2) {
	ROW_LIST *rl1 = ((ROW_LIST *)p1);
	ROW_LIST *rl2 = ((ROW_LIST *)p2);
	return strcmp(rl1->token, rl2->token);
}

int comprow(const void *m1, const void *m2) {
	 int *r1 = (int *)m1;
	 int *r2 = (int *)m2;
	 return *r1 - *r2;
}

int add_row_list(ROW_LIST *src, ROW_LIST *dst) {
	int i;
	for (i = 0; i < src->nb_row; i++)
		if (dst == NULL) {
			dst = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
			dst->row = (int *)calloc(1, sizeof(int));
			dst->row[dst->nb_row] = src->row[i];
			dst->nb_row++;
		}
		else if (bsearch(&(src->row[i]), dst->row, dst->nb_row, sizeof(int), comprow) == NULL) {
			dst->row = (int *)realloc(dst->row, (dst->nb_row + 1) * sizeof(int));
			dst->row[dst->nb_row] = src->row[i];
			dst->nb_row++;
		}
	return dst->nb_row;
}
