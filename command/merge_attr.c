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

//<<<<<<< HEAD
//char *concatenate_str(char *a, char *b, char *c){
//  int size = strlen(a) + strlen(b) + strlen(c) + 1;
//  char *str = malloc(size);
//  strcpy (str, a);
//  strcat (str, b);
//  strcat (str, c);
//
//  return str;
//}
//=======
/*
 * global variables declaration
 */
extern COLUMN **column;
extern int nb_column;
//>>>>>>> c2d0f09934535ced750c90f997c24b68c0686fbc

__attribute__ ((visibility ("default")))
GTF_DATA *merge_attr(GTF_DATA *gtf_data, char *features, char *keys, char *dest_key, char *sep) {

//<<<<<<< HEAD
//	int i, j, nb, ok;
//=======
	int i, j, k, c, p, nb_attr, new_size, ok;
	int nb_requested_key = 0;
//>>>>>>> c2d0f09934535ced750c90f997c24b68c0686fbc
	char **key_list;
	char *str_concat = NULL;
	char none[] = ".";
	char void_str[] = "";
	char * new_buffer = NULL;

	typedef struct hit{
		char *key;
		char *value;
	} hit;

//<<<<<<< HEAD
//	key_list = (char **)malloc(2 * sizeof(char *));
//=======
	//key_list = (char **)malloc(2 * sizeof(char *));
	//nb_requested_key = split_ip(&key_list, strdup(keys), ",");
	nb_requested_key = split_ip(&key_list, dupstring(keys, __FILE__, __func__, __LINE__), ",");

	//hit *hits = malloc(nb_requested_key * sizeof *hits);
	hit *hits = bookmem(nb_requested_key, sizeof(hit *), __FILE__, __func__, __LINE__);
//>>>>>>> c2d0f09934535ced750c90f997c24b68c0686fbc

	/*
	 * reserve memory for the GTF_DATA structure to return
	*/

	GTF_DATA *ret = clone_gtf_data(gtf_data);

	GTF_ROW *row;

	ATTRIBUTE *pattr;

	for (i = 0; i < ret->size; i++)
	{

		for(p=0; p < nb_requested_key; p++)
		{
			hits[p].key = key_list[p];
			hits[p].value = none;
		}

		row = ret->data[i];
		ok = (*features == '*');
		if (!ok) ok = (strstr(features, row->field[2]) != NULL);
		if (ok)
		{
			nb_attr = row->attributes.nb;

			for(k=0; k < nb_requested_key; k++)
			{
				for(j=0; j < nb_attr; j++)
				{
					pattr = row->attributes.attr[j];
					if (strcmp(key_list[k], pattr->key) == 0)
					{
						for(p=0; p < nb_requested_key; p++)
						{
							if(strcmp(pattr->key, hits[p].key) ==0)
							{
								//hits[p].value = strdup(pattr->value);
								hits[p].value = dupstring(pattr->value, __FILE__, __func__, __LINE__);
							}
						}
					}
				}

				str_concat = void_str;

				for(p=0; p < nb_requested_key; p++){

					if(strcmp(hits[p].value, ".") ==0)
					{
						for (c = 0; c < nb_column; c++)
						{

						    if(!strcmp(column[c]->name, hits[p].key))
						    {
						    	//hits[p].value = strdup(row->field[c]);
						    	hits[p].value = dupstring(row->field[c], __FILE__, __func__, __LINE__);
						    	break;
						    }
						}

					}
					if(p < (nb_requested_key -1))
					{
						new_size = strlen(hits[p].value)  + strlen(str_concat) + strlen(sep) + 1;
						//new_buffer = (char *)malloc(new_size);
						new_buffer = (char *)bookmem(new_size, sizeof(char), __FILE__, __func__, __LINE__);
						strcpy(new_buffer, str_concat);
						strcat(new_buffer, hits[p].value);
						strcat(new_buffer, sep);
						str_concat = new_buffer;
					}else
					{
						new_size = strlen(hits[p].value)  + strlen(str_concat) + 1;
						//new_buffer = (char *)malloc(new_size);
						new_buffer = (char *)bookmem(new_size, sizeof(char), __FILE__, __func__, __LINE__);
						strcpy(new_buffer, str_concat);
						strcat(new_buffer, hits[p].value);
						str_concat = new_buffer;
					}

				}



				pattr = pattr->next;

			}

			add_attribute(row, dest_key, str_concat);
		}

		update_attribute_table(row);
	}

	return ret;

}

