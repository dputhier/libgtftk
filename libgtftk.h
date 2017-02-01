/*
 * libgtftk.h
 *
 *  Created on: Jan 10, 2017
 *      Author: fafa
 *
 *  Header file for the gtftk library.
 *  Contains all the structure definitions and the prototype declarations.
 */

#ifndef GTFTOOLKIT_GTFTK_SRC_LIB_LIBGTFTK_H_
#define GTFTOOLKIT_GTFTK_SRC_LIB_LIBGTFTK_H_

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <search.h>
#include <string.h>

/*
 * This structure describes the input (i.e. a GTF file or a gzipped GTF file).
 * It is created by the get_gtf_reader function in get_reader.c source file.
 * gzFile or plainfile are set depending on the kind of input (gzip or plain).
 * Even if these two elements are exclusive, they are not in an union to be
 * sure to be compatible with the most of native interfaces.
 */
typedef struct GTF_READER {
	// The file name with its path (or "-" for standard input)
	char *filename;

	// A boolean to tell if it is a gzipped file
	int gz;

	// The gzip file descriptor
	gzFile gzfile;

	// The plain file descriptor
	FILE *plainfile;
} GTF_READER;

/*
 * The structure that represents an attribute (key/value) in the last column of
 * a GTF file.
 */
typedef struct ATTRIBUTE {
	char *key, *value;
} ATTRIBUTE;

/*
 * A set of ATTRIBUTE
 */
typedef struct ATTRIBUTES {
	ATTRIBUTE **attr;
	int nb;
} ATTRIBUTES;

/*
 * A structure to store a row from a GTF file. data is a table of 9 pointers
 * to store the 9 values in a row. Each pointer has a specific type, depending
 * on the column. For example, data[2] ("feature" column) is a char *, data[3]
 * ("start" column) is an int * and data[8] ("attributes" column) is an
 * ATTRIBUTES *.
 */
typedef struct GTF_ROW {
	// a pointer on the 9 values of a GTF file row
	void **data;

	// the rank number of the row in the GTF file
	int rank;
} GTF_ROW;

/*
 * The same structure as ROW_CHAR but with string values. This structure is
 * used by some client languages (python) that can't handle easily the tables
 * with mixed pointer types. A convenient function is provided in the library
 * to translate GTF_ROW into GTF_ROW_CHAR : gtf_row_to_char
 */
typedef struct GTF_ROW_CHAR {
	// a pointer on the 9 string-ed values of a GTF file row
	char **data;

	// the rank number of the row in the GTF file
	int rank;
} GTF_ROW_CHAR;

/*
 * This is the structure that holds data in GTF format. It is also the
 * structure used as input/output for most of the functions of the library. To
 * start using the library, one must call the loadGTF() function with a GTF
 * file name in parameter and gets a pointer on a GTF_DATA in return. Then, all
 * the other functions must be called with this GTF_DATA pointer as input.
 * Their result can be another GTF_DATA pointer that can be used as input for
 * another function of the library.
 */
typedef struct GTF_DATA {
	// the number of rows
	int size;

	// a table of pointers on the rows
	GTF_ROW **data;
} GTF_DATA;

/*
 * This structure represents an index on a column, or on an attribute of the
 * last column.
 */
typedef struct INDEX {
	// the name of a column (feature, seqid, ...) or of an attribute (gene_id,
	// transcript_id, ...)
	char *key;

	/* the pointer on a binary tree created with tsearch C function in
	 * the index_row function in column.c source file. This tree contains
	 * ROW_LIST elements described later in this file and that contains a
	 * token and the associated list of row numbers (the rows containing the
	 * token as the value of the key (a column name or an attribute name).
	 */
	void *data;
} INDEX;

/*
 * This is a class-like structure that modelize a column of a GTF file. It
 * contains obvious information like the name of the column, its type and the
 * indexes made on it. It contains also pointers on functions that are intended
 * to work with the corresponding values in the GTF data. There is a particular
 * function for each type of data. For example, the convert function takes a
 * string value and a default value as parameters and returns and integer
 * pointer for columns that contains integers (start, end ...), a float pointer
 * for "score" column and an ATTRIBUTES pointer for the last column. Pointers
 * on functions are initialized in the make_column function in the column.c
 * source file, depending on the type of each column.
 */
typedef struct COLUMN {
	// the rank number of the column
	int num;

	/* the type of the column, that is 'A' for attributes, 'I' for integer,
	 * 'C' for char, 'F' for float and 'S' for string
	 */
	char type;

	/* the column name : seqid, source, feature, start, end, score, strand,
	 * phase or attributes
	 */
	char *name;

	// the default value to print if no value is available (".", -1 ...)
	void *default_value;

	/* a table of indexes. It contains only one pointer for each column except
	 * the attributes column for which there is as needed indexes (one can
	 * index data on several attributes)
	 */
	INDEX **index;

	// the number of indexes in the previous table
	int nb_index;

	/* this function convert a string value (read in a GTF file) in the
	 * appropriate type (the column type : int, float, char, ATTRIBUTES ...).
	 * The corresponding functions are implemented in the column.c source file.
	 * Parameters : the value and the default value as strings
	 * Returns : the typed value
	 */
	void *(*convert)(char *token, void *def);

	/* do the reverse job of the previous function : convert a typed value in
	 * a column into a string value (for indexing in index_gtf function in
	 * the load_gtf.c source file).
	 * Parameters : the typed value and default value (int, char ...)
	 * Returns : the char representation of the value (or the default value)
	 */
	char *(*convert_to_string)(void *value, void *def);

	/*
	 * this function add a row in the ROW_LIST element corresponding to the
	 * value in the specified index. This function is called by the index_gtf
	 * function that is implemented in the load_gtf.c source file.
	 * Parameters:
	 * 		row_nb:		the row number
	 * 		value:		the concerned value
	 * 		index:		the index in which the value to be searched and modified
	 */
	void (*index_row)(int row_nb, char *value, INDEX *index);

	/*
	 * prints in output a typed token with the given delimiter according to
	 * the type of the column col. The value of col is useless for attributes.
	 * Parameters:
	 * 		token:		a typed value
	 * 		output:		the output to print in (a file or stdout)
	 * 		col:		a column or NULL for the attributes column
	 * 		delim:		a char as delimiter ('\t' for a GTF output)
	 */
	void (*print)(void *token, FILE *output, void *col, char delim);
} COLUMN ;

/*
 * A list of row numbers associated with a token (the values in the 8 first
 * columns or the values associated to an attribute in the last column). This
 * structure is used in the indexes as elements.
 */
typedef struct ROW_LIST {
	/* the token that is contained in the rows. For example, this can be "gene"
	 * or "transcript" for an index on the column feature, or "protein_coding"
	 * and "lincRNA" for an index on the attribute "gene_biotype".
	 */
	char *token;

	// the number of rows
	int nb_row;

	// the table of row numbers
	int *row;
} ROW_LIST;

/*
 * Prototypes for the visible functions (callable by external cient)
 */
GTF_DATA *load_GTF(char *input);
GTF_DATA *select_by_key(GTF_DATA *gtf_data, char *key, char *value, int not);
void print_gtf_data(GTF_DATA *gtf_data);
GTF_DATA *select_by_transcript_size(GTF_DATA *gtf_data, int min, int max);
GTF_ROW_CHAR *gtf_row_to_char(GTF_ROW *row);
GTF_DATA *select_by_number_of_exon(GTF_DATA *gtf_data, int min, int max);
GTF_DATA *select_by_genomic_location(GTF_DATA *gtf_data, char *chr, int begin_gl, int end_gl);

#endif /* GTFTOOLKIT_GTFTK_SRC_LIB_LIBGTFTK_H_ */
