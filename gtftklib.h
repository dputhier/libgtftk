/*
 * gtftklib.h
 *
 *  Created on: Apr 4, 2016
 *      Author: fafa
 */
#include <stdio.h>
#include <stdlib.h>

#ifndef GTFTOOLKIT_SRC_LIB_GTFTKLIB_H_
#define GTFTOOLKIT_SRC_LIB_GTFTKLIB_H_
#ifndef GTFTOOLKIT_GTFTK_SRC_LIB_LIBGTFTK_H_
#ifndef GTFTK_H_
typedef struct FTP_DATA {
	char *URL;
	char *base_dir;
	char *user;
	char *password;
	void *data;
} FTP_DATA;

typedef struct ROW {
	void **data;
	int nb_attributes;
} ROW ;

typedef struct INDEX {
	char *key;
	void *data;
} INDEX ;

typedef struct COLUMN {
	int num;
	char type;
	char *name;
	void *default_value;
	INDEX **index;
	int nb_index;
	void *(*convert)(char *, void *, ROW *);
	char *(*convert_to_string)(void *, void *);
	void (*make_index)(int , ROW *, struct COLUMN *);
	void (*print)(void *, FILE *, void *, char);
} COLUMN ;

typedef struct GTF_READER {
	char *filename;
	int gz;
	gzFile gzfile;
	FILE *plainfile;
} GTF_READER;
#endif /* GTFTK_H_ */
#endif /* GTFTOOLKIT_GTFTK_SRC_LIB_LIBGTFTK_H_ */
typedef struct GTF_ROW_ATTRIBUTES {
	int size;
	char **key;
	char **value;
} GTF_ROW_ATTRIBUTES;

typedef struct GTF_ROW {
	char **field;
	GTF_ROW_ATTRIBUTES *attribute;
} GTF_ROW;

typedef struct GTF_RESULT {
	int size;
	GTF_ROW **data;
} GTF_RESULT;

typedef struct GTF_DATA {
	int size;
	GTF_ROW **data;
} GTF_DATA;

typedef struct TAB_RESULT {
	int size;
	char **data;
} TAB_RESULT;

typedef struct FASTA {
	char *header;
	char *sequence;
} FASTA;

typedef struct FASTA_RESULT {
	int size;
	union {
		FASTA **data;
		char *fasta_file;
	};
} FASTA_RESULT;

extern FTP_DATA *ftp_data;

GTF_RESULT *select_by_key_lib(char *gtf_filename, char *key, char *value, int not);
GTF_RESULT *select_by_key_from_gtfresult_lib(GTF_RESULT *gtf_result, char *key, char *value, int not);
GTF_RESULT *select_by_genomic_location_lib(char *gtf_filename, char *chr, int start, int end);
TAB_RESULT *extract_data_lib(char *gtf_filename, char *key);
FASTA_RESULT *get_fasta_lib(char *gtf_filename, char *genome_file, int intron, int rc);
TAB_RESULT *get_feature_list_lib(char *gtf_filename);
TAB_RESULT *get_attribute_list_lib(char *gtf_filename);
TAB_RESULT *get_seq_list_lib(char *gtf_filename);
void free_ptr(void *p);
void free_gtfresult_ptr(GTF_RESULT *p);
GTF_RESULT *select_by_transcript_size_lib(char *gtf_filename, int min, int max);

#endif /* GTFTOOLKIT_SRC_LIB_GTFTKLIB_H_ */
