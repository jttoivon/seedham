/*
	Copyright 2010-2011 Kamiar Kanani <kamiar.kanani@gmail.com>

    This file is part of multifast.

    multifast is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    multifast is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with multifast.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "node.h"

/* reallocation step for node::matched_strings array */
#define REALLOC_CHUNK_MATCHSTR 8

/* reallocation step for node::outgoing array */
#define REALLOC_CHUNK_OUTGOING 8
/* TODO: For different depth of node, number of outgoing edges differs
considerably, It is efficient to use different chunk size for 
different depths */

/* Private Functions */
void   node_init              (NODE * thiz);
int    node_edge_compare      (const void * l, const void * r);
int    node_has_matchstr      (NODE * thiz, STRING * newstr);



/******************************************************************************
FUNCTION: node_create

DESCRIPTION:
	Create the node
******************************************************************************/
struct node * node_create(void)
{
	NODE * thiz;

	thiz = (NODE *) malloc (sizeof(NODE));
	node_init(thiz);
	#ifdef DEBUG_DISPLAY_AC
	node_assign_id(thiz);
	#endif

	return thiz;
}


/******************************************************************************
FUNCTION: node_init

DESCRIPTION:
	Initialize node
******************************************************************************/
void node_init(NODE * thiz)
{
	memset(thiz, 0, sizeof(NODE));

	thiz->outgoing_max = REALLOC_CHUNK_OUTGOING;
	thiz->outgoing = (struct edge *) malloc (thiz->outgoing_max*sizeof(struct edge));

	thiz->matched_strings_max = REALLOC_CHUNK_MATCHSTR;
	thiz->matched_strings = (STRING *) malloc (thiz->matched_strings_max*sizeof(STRING));
}


/******************************************************************************
FUNCTION: node_release

DESCRIPTION:
	Release node
******************************************************************************/
void node_release(NODE * thiz)
{
	free(thiz->matched_strings);
	free(thiz->outgoing);
	free(thiz);
}


/******************************************************************************
FUNCTION: node_find_next

DESCRIPTION: 
	Find out the target node from here for a given Alpha.
	it uses linear search. this function is used in
	preprocessig stage in which edge array is not sorted.
******************************************************************************/
NODE * node_find_next(NODE * thiz, ALPHA alpha)
{
	int i;

	for (i=0; i < thiz->outgoing_degree; i++)
	{
		if(thiz->outgoing[i].alpha == alpha)
			return (thiz->outgoing[i].next);
	}

	return NULL;
}


/******************************************************************************
FUNCTION: node_findbs_next

DESCRIPTION: 
	Find out the target node from here for a given Alpha.
	it uses Binary Search. this function is used after 
	preprocessing stage in which we sort edges.
******************************************************************************/
NODE * node_findbs_next (NODE * thiz, ALPHA alpha)
{
	int min, max, mid;
	ALPHA amid;

	min = 0;
	max = thiz->outgoing_degree - 1;

	while (min <= max)
	{
		mid = (min+max) >> 1;
		amid = thiz->outgoing[mid].alpha;
		if (alpha > amid)
			min = mid + 1;
		else if (alpha < amid)
			max = mid - 1;
		else
			return (thiz->outgoing[mid].next);
	}

	return NULL;
}


/******************************************************************************
FUNCTION: node_has_matchstr

DESCRIPTION:
	Determine if a final node contains an string in its accepted string list
	Return values: 1 = has, 0 = hasn't
******************************************************************************/
int node_has_matchstr (NODE * thiz, STRING * newstr)
{
	int i;
	unsigned int j;
	STRING * str;

	for (i=0; i < thiz->matched_strings_num; i++)
	{
		str = &thiz->matched_strings[i];

		if (str->length != newstr->length)
			continue;

		for (j=0; j<str->length; j++)
			if(str->str[j] != newstr->str[j])
				continue;

		if (j == str->length)
			return 1;
	}

	return 0;
}


/******************************************************************************
FUNCTION: node_create_next

DESCRIPTION:
	Create next node for the given alpha
******************************************************************************/
NODE * node_create_next (NODE * thiz, ALPHA alpha)
{
	NODE * next;
	next = node_find_next (thiz, alpha);
	if (next)
	/* The edge already exists */
		return NULL;
	/* Otherwise register new edge */
	next = node_create ();
	node_register_outgoing(thiz, next, alpha);

	return next;
}


/******************************************************************************
FUNCTION: node_register_matchstr

DESCRIPTION:
	Add the string to the list of accepted strings
******************************************************************************/
void node_register_matchstr(NODE * thiz, STRING * str)
{
	/* Check if the new string already exists in the node list */
	if (node_has_matchstr(thiz, str))
		return;

	/* Manage memory */
	if (thiz->matched_strings_num >= thiz->matched_strings_max)
	{
		thiz->matched_strings_max += REALLOC_CHUNK_MATCHSTR;
		thiz->matched_strings = (STRING *) realloc 
			(thiz->matched_strings, thiz->matched_strings_max*sizeof(STRING));
	}

	thiz->matched_strings[thiz->matched_strings_num].str = str->str;
	thiz->matched_strings[thiz->matched_strings_num].length = str->length;
	thiz->matched_strings[thiz->matched_strings_num].id = str->id;
	thiz->matched_strings_num++;
}


/******************************************************************************
FUNCTION: node_register_outgoing

DESCRIPTION:
	Establish an edge between two nodes
******************************************************************************/
void node_register_outgoing (NODE * thiz, NODE * next, ALPHA alpha)
{
	if(thiz->outgoing_degree >= thiz->outgoing_max)
	{
		thiz->outgoing_max += REALLOC_CHUNK_OUTGOING;
		thiz->outgoing = (struct edge *) realloc 
			(thiz->outgoing, thiz->outgoing_max*sizeof(struct edge));
	}

	thiz->outgoing[thiz->outgoing_degree].alpha = alpha;
	thiz->outgoing[thiz->outgoing_degree++].next = next;
}


/******************************************************************************
FUNCTION: node_assign_id

DESCRIPTION: 
	every node has a unique ID only for debugging purpose. it is not
	necessary in normal execution.
******************************************************************************/
#ifdef DEBUG_DISPLAY_AC
void node_assign_id (NODE * thiz)
{
	static int unique_id = 100;
	thiz->id = unique_id ++;
}
#endif

/******************************************************************************
FUNCTION: node_edge_compare

DESCRIPTION: Comparison function for qsort. see man qsort.
******************************************************************************/
int node_edge_compare (const void * l, const void * r)
{
	/* 
	According to man page:
	The comparison function must return an integer less than, equal to, or 
	greater than zero if the first argument is considered to be respectively 
	less than,  equal  to,  or greater than the second.  If two members 
	compare as equal, their order in the sorted array is undefined.

	NOTE: Because edge alphas are unique in every node we ignore equivalence
	case for sake of performance.
	*/
	if ( ((struct edge *)l)->alpha >= ((struct edge *)r)->alpha )
		return 1;
	else
		return -1;
}


/******************************************************************************
FUNCTION: node_sort_edges

DESCRIPTION: sort edges alpha.
******************************************************************************/
void node_sort_edges (NODE * thiz)
{
	qsort ((void *)thiz->outgoing, thiz->outgoing_degree, sizeof(struct edge),
			node_edge_compare);
}

