/*
 * libgtftk.h
 *
 *  Created on: Jan 10, 2017
 *      Author: fafa
 */

#ifndef GTFTOOLKIT_GTFTK_SRC_LIB_LIBGTFTK_H_
#define GTFTOOLKIT_GTFTK_SRC_LIB_LIBGTFTK_H_

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <search.h>
#include <string.h>

#define SEQID "seqid"

const char *COLUMN_NAME[1];

typedef struct GTF_READER {
	char *filename;
	int gz;
	gzFile gzfile;
	FILE *plainfile;
} GTF_READER;

typedef struct ATTRIBUTE {
	char *key, *value;
} ATTRIBUTE;

typedef struct ATTRIBUTES {
	ATTRIBUTE **attr;
	int nb;
} ATTRIBUTES;

typedef struct GTF_ROW {
	void **data;
	int nb_attributes;
} GTF_ROW;

typedef struct GTF_DATA {
	int size;
	GTF_ROW **data;
} GTF_DATA;

typedef struct INDEX {
	char *key;
	void *data;
} INDEX;

typedef struct COLUMN {
	int num;
	char type;
	char *name;
	void *default_value;
	INDEX **index;
	int nb_index;
	void *(*convert)(char *, void *, GTF_ROW *);
	char *(*convert_to_string)(void *, void *);
	void (*make_index)(int , GTF_ROW *, struct COLUMN *);
	void (*print)(void *, FILE *, void *, char);
} COLUMN ;

typedef struct ROW_LIST {
	char *token;
	int nb_row;
	int *row;
} ROW_LIST;

GTF_DATA *load_GTF(char *input);
GTF_DATA *select_by_key(GTF_DATA *gtf_data, char *key, char *value, int not);
void print_gtf_data(GTF_DATA *gtf_data);

#endif /* GTFTOOLKIT_GTFTK_SRC_LIB_LIBGTFTK_H_ */
