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
extern void add_attribute(GTF_ROW *row, char *key, char *value);
extern int split_ip(char ***tab, char *s, char *delim);
extern void *bookmem(int nb, int size, char *file, const char *func, int line);
extern char *dupstring(const char *s, char *file, const char *func, int line);
extern void *rebookmem(void *ptr, int size, char *file, const char *func, int line);
extern void freemem(void *ptr, char *file, const char *func, int line);

/*
 * global variables declaration
 */
extern COLUMN **column;
extern int nb_column;

__attribute__ ((visibility ("default")))
GTF_DATA *merge_attr(GTF_DATA *gtf_data, char *features, char *keys, char *dest_key, char *sep) {
	int i, j, k, c, p, nb_attr, new_size, ok;
	int nb_requested_key = 0;
	char **key_list;
	char *none = dupstring(".", __FILE__, __func__, __LINE__);
	char *new_buffer;

	typedef struct hit {
		char *key;
		char *value;
	} hit;

	nb_requested_key = split_ip(&key_list, dupstring(keys, __FILE__, __func__, __LINE__), ",");

	hit *hits = bookmem(nb_requested_key, sizeof(hit *), __FILE__, __func__, __LINE__);

	/*
	 * reserve memory for the GTF_DATA structure to return
	*/

	GTF_DATA *ret = clone_gtf_data(gtf_data);

	GTF_ROW *row;

	ATTRIBUTE *pattr;

	for (i = 0; i < ret->size; i++)	{
		new_size = 1;
		new_buffer = (char *)bookmem(1, sizeof(char), __FILE__, __func__, __LINE__);
		*new_buffer = 0;
		for (p = 0; p < nb_requested_key; p++)	{
			hits[p].key = key_list[p];
			hits[p].value = none;
		}

		row = ret->data[i];
		ok = (*features == '*');
		if (!ok) ok = (strstr(features, row->field[2]) != NULL);
		if (ok)	{
			nb_attr = row->attributes.nb;
			for (k = 0; k < nb_requested_key; k++)	{
				for (j = 0; j < nb_attr; j++) {
					pattr = row->attributes.attr[j];
					if (strcmp(key_list[k], pattr->key) == 0) {
						for (p = 0; p < nb_requested_key; p++) {
							if (strcmp(pattr->key, hits[p].key) == 0) {
								//hits[p].value = strdup(pattr->value);
								hits[p].value = dupstring(pattr->value, __FILE__, __func__, __LINE__);
							}
						}
					}
				}
				for (p = 0; p < nb_requested_key; p++) {
					if (strcmp(hits[p].value, ".") == 0)	{
						for (c = 0; c < nb_column; c++)	{
						    if (!strcmp(column[c]->name, hits[p].key)) {
						    	//hits[p].value = strdup(row->field[c]);
						    	hits[p].value = dupstring(row->field[c], __FILE__, __func__, __LINE__);
						    	break;
						    }
						}

					}
					if (p < (nb_requested_key -1)) {
						new_size += strlen(hits[p].value) + strlen(sep);
						//new_buffer = (char *)malloc(new_size);
						new_buffer = (char *)rebookmem(new_buffer, new_size, __FILE__, __func__, __LINE__);
						strcat(new_buffer, hits[p].value);
						strcat(new_buffer, sep);
					}
					else {
						new_size += strlen(hits[p].value);
						//new_buffer = (char *)malloc(new_size);
						new_buffer = (char *)rebookmem(new_buffer, new_size, __FILE__, __func__, __LINE__);
						strcat(new_buffer, hits[p].value);
					}
				}
				pattr = pattr->next;
			}
			add_attribute(row, dest_key, new_buffer/*str_concat*/);
			freemem(new_buffer, __FILE__, __func__, __LINE__);
		}
		update_attribute_table(row);
	}

	return ret;
}
