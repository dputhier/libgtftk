/*
 * merge_attr.c
 *
 *  Created on: Nov 14, 2017
 *      Author: dputhier (inspired from fafa code...)
 *		Objective: Merge two keys into a destination key using sep as separator.
 *
 */

#include "libgtftk.h"

/*
 * external functions declaration
 */
extern GTF_DATA *clone_gtf_data(GTF_DATA *gtf_data);
extern int update_attribute_table(GTF_ROW * row);
extern  int split_ip(char ***tab, char *s, char *delim);
extern void add_attribute(GTF_ROW *row, char *key, char *value);

char *concatenate_str(char *a, char *b, char *c){
  int size = strlen(a) + strlen(b) + strlen(c) + 1;
  char *str = malloc(size);
  strcpy (str, a);
  strcat (str, b);
  strcat (str, c);

  return str;
}

/*
 * global variables declaration
 */
extern COLUMN **column;

__attribute__ ((visibility ("default")))
GTF_DATA *merge_attr(GTF_DATA *gtf_data, char *features, char *keys, char *dest_key, char *sep) {

	int i, j, nb, ok;
	int nb_keys;
	char **key_list;
	char *value_1;
	char *value_2;
	char *str_concat;


	key_list = (char **)malloc(2 * sizeof(char *));
	nb_keys = split_ip(&key_list, strdup(keys), ",");

	/*
	 * reserve memory for the GTF_DATA structure to return
	*/

	GTF_DATA *ret = clone_gtf_data(gtf_data);

	GTF_ROW *row;

	ATTRIBUTE *pattr;

	for (i = 0; i < ret->size; i++) {

		value_1 = "";
		value_2 = "";

		row = ret->data[i];
		ok = (*features == '*');
		if (!ok) ok = (strstr(features, row->field[2]) != NULL);
		if (ok) {
			nb = row->attributes.nb;

			for(j=0; j < nb; j++){

				pattr = row->attributes.attr[j];

				if (strcmp(key_list[0], pattr->key) == 0)  {
					if(strcmp(value_1, "") == 0){
						value_1 = strdup(pattr->value);
					}
				}

				if (strcmp(key_list[1], pattr->key) == 0) {
					if(strcmp(value_2, "") == 0){
						value_2 = strdup(pattr->value);
					}
				}

				if(strcmp(value_1, "") != 0 && strcmp(value_2, "") != 0){
					str_concat = concatenate_str(value_1, sep, value_2);
					add_attribute(row, dest_key, str_concat);
				}
				pattr = pattr->next;

			}

		}

		update_attribute_table(row);
	}

	return ret;

}

