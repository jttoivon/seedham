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
#include <stdlib.h>
#include <string.h>

#include "aho_corasick.h"

/* Allocation step for automata::all_nodes array */
#define REALLOC_CHUNK_ALLNODES 200


/* Private Functions */
void   ac_automata_register_nodeptr  (AC_AUTOMATA * thiz, NODE * node);
void   ac_automata_set_failure       (AC_AUTOMATA * thiz, NODE * node, ALPHA * alphas);
void   ac_automata_dfs_traverse      (AC_AUTOMATA * thiz, NODE * node, ALPHA * alphas);
void   ac_automata_union_matchstrs   (NODE * node);



/******************************************************************************
FUNCTION: ac_automata_init
PARAMS:
	AC_AUTOMATA * thiz: Pointer to the automata
	MATCH_CALBACK mc: callback function

DESCRIPTION:
	Initialize automata; allocate memories for automata.
******************************************************************************/
void ac_automata_init (AC_AUTOMATA * thiz, MATCH_CALBACK mc)
{
	memset (thiz, 0, sizeof(AC_AUTOMATA));
	thiz->root = node_create ();
	thiz->all_nodes_max = REALLOC_CHUNK_ALLNODES;
	thiz->all_nodes = (NODE **) malloc (thiz->all_nodes_max*sizeof(NODE *));
	thiz->match_callback = mc;
	ac_automata_register_nodeptr (thiz, thiz->root);
	ac_automata_reset (thiz);
	thiz->total_strings = 0;
	thiz->accept_strings = 1;
}


/******************************************************************************
FUNCTION: ac_automata_reset

DESCRIPTION: 
	reset automata state. when you finished with the input string,
	you must reset automata state for new input, otherwise it will not work.
******************************************************************************/
void ac_automata_reset (AC_AUTOMATA * thiz)
{
	thiz->current_node = thiz->root;
	thiz->base_position = 0;
}


/******************************************************************************
FUNCTION: ac_automata_register_nodeptr

DESCRIPTION:
	Add node pointer to the list (all_nodes) for further use.
******************************************************************************/
void ac_automata_register_nodeptr (AC_AUTOMATA * thiz, NODE * node)
{
	if(thiz->all_nodes_num >= thiz->all_nodes_max)
	{
		thiz->all_nodes_max += REALLOC_CHUNK_ALLNODES;
		thiz->all_nodes = (NODE**)realloc (thiz->all_nodes, thiz->all_nodes_max*sizeof(NODE *));
	}

	thiz->all_nodes[thiz->all_nodes_num++] = node;
}


/******************************************************************************
FUNCTION: ac_automata_add_string

RETUERNS:
	ACERR_NONE on success
	ACERR_LONG_STRING when string length is longer than AC_PATTRN_MAX_LENGTH
	ACERR_DUPLICATE_STRING on duplicte strings
	ACERR_ZERO_STRING on zero length string

DESCRIPTION:
	Add string to the automata.
	CAUTION: If the given string be larger than AC_PATTRN_MAX_LENGTH, it will
	be cropped without any warning.
******************************************************************************/
AC_ERROR ac_automata_add_string (AC_AUTOMATA * thiz, STRING * str)
{
	unsigned int i;
	NODE * n = thiz->root;
	NODE * next;
	ALPHA alpha;

	if(!thiz->accept_strings)
		return ACERR_STRING_CLOSED;

	if (!str->length)
		return ACERR_ZERO_STRING;

	if (str->length > AC_PATTRN_MAX_LENGTH)
		return ACERR_LONG_STRING;

	for (i=0; i<str->length; i++)
	{
		alpha = str->str[i];
		if ((next = node_find_next(n, alpha)))
		{
			n = next;
			continue;
		}
		else
		{
			next = node_create_next(n, alpha);
			next->depth = n->depth + 1;
			n = next;
			ac_automata_register_nodeptr(thiz, n);
		}
	}

	if(n->final)
		return ACERR_DUPLICATE_STRING;

	n->final = 1;
	node_register_matchstr(n, str);
	thiz->total_strings++;

	return ACERR_NONE;
}


/******************************************************************************
FUNCTION: ac_automata_release

DESCRIPTION:
	Release all alocated memories to the automata
******************************************************************************/
void ac_automata_release (AC_AUTOMATA * thiz)
{
	unsigned int i;
	NODE * n;

	for (i=0; i < thiz->all_nodes_num; i++)
	{
		n = thiz->all_nodes[i];
		node_release(n);
	}

	free(thiz->all_nodes);
}


/******************************************************************************
FUNCTION: ac_automata_union_matchstrs

DESCRIPTION:
	Accepted string in any node consists of its own strings plus strings of
	its failure node. this function collects all strings in failure node.
******************************************************************************/
void ac_automata_union_matchstrs (NODE * node)
{
	unsigned int i;
	NODE * m = node;

	while ((m = m->failure_node))
	{
		for (i=0; i < m->matched_strings_num; i++)
			node_register_matchstr(node, &(m->matched_strings[i]));

		if (m->final)
			node->final = 1;
	}
}


/******************************************************************************
FUNCTION: ac_automata_set_failure

DESCRIPTION:
	find failure node for the given node.
******************************************************************************/
void ac_automata_set_failure (AC_AUTOMATA * thiz, NODE * node, ALPHA * alphas)
{
	unsigned int i, j;
	NODE * m;

	for (i=1; i < node->depth; i++)
	{
		m = thiz->root;
		for (j=i; j < node->depth && m; j++)
			m = node_find_next (m, alphas[j]);
		if (m)
		{
			node->failure_node = m;
			break;
		}
	}
	if (!node->failure_node)
		node->failure_node = thiz->root;
}


/******************************************************************************
FUNCTION: ac_automata_dfs_traverse

DESCRIPTION:
	Traverse all automata nodes using DFS (Deapth First Search), meanwile
	it set the failure node for every node it passes through.
*****************************************************************************/
void ac_automata_dfs_traverse (AC_AUTOMATA * thiz, NODE * node, ALPHA * alphas)
{
	unsigned int i;
	NODE * next;

	for (i=0; i < node->outgoing_degree; i++)
	{
		alphas[node->depth] = node->outgoing[i].alpha;
		next = node->outgoing[i].next;

		/* At every node look for its failure node */
		ac_automata_set_failure (thiz, next, alphas);

		/* Recursively call itself to traverse all nodes */
		ac_automata_dfs_traverse (thiz, next, alphas);
	}
}


/******************************************************************************
FUNCTION: ac_automata_locate_failure

DESCRIPTION:
	Locate the failure node for all nodes.
******************************************************************************/
void ac_automata_locate_failure (AC_AUTOMATA * thiz)
{
	unsigned int i;
	ALPHA alphas[AC_PATTRN_MAX_LENGTH];
	NODE * node;

	ac_automata_dfs_traverse (thiz, thiz->root, alphas);

	for (i=0; i < thiz->all_nodes_num; i++)
	{
		node = thiz->all_nodes[i];
		ac_automata_union_matchstrs (node);
		node_sort_edges (node);
	}

	thiz->accept_strings = 0; /* do not accept strings any more */
}


/******************************************************************************
FUNCTION: ac_automata_search

DESCRIPTION:
	Search for patterns inside the given input, as it finds a match
	it will call the Callback functions to report it to caller.
******************************************************************************/
void ac_automata_search (AC_AUTOMATA * thiz, STRING * str, void * param)
{
	unsigned long position;
	NODE * current;
	NODE * next;

	if(thiz->accept_strings)
		/* you must call ac_automata_locate_failure() first */
		return;

	position = 0;
	/* reload status variable(s) */
	current = thiz->current_node;

	/* This is the main search loop. it must be keep 
	   as lightwaight as possible. */
	while (position < str->length)
	{
		if(!(next = node_findbs_next(current, str->str[position])))
		{
			if(current->failure_node /* we are not in root */)
				current = current->failure_node;
			else
				position++;
		}
		else
		{
			current = next;
			position++;
		}

		if(current->final && next)
		/* We check 'next' to find out if we came here after a alpha transition 
		   or due to a fail. in second case we should not report matching,
		   because it was reported in previous node */
		{
			thiz->match.position = position + thiz->base_position;
			thiz->match.match_num = current->matched_strings_num;
			thiz->match.matched_strings = current->matched_strings;
			/* do callback: we found a match */
			if (thiz->match_callback(&thiz->match, param))
				break;
		}
	}

	/* save status variables */
	thiz->current_node = current;
	thiz->base_position += position;
}


/******************************************************************************
FUNCTION: ac_automata_dbg_show

DESCRIPTION:
	Print automata to autput in human readable form.
	used for debugging purpose.
******************************************************************************/
#ifdef DEBUG_DISPLAY_AC
void ac_automata_dbg_show (AC_AUTOMATA * thiz)
{
	unsigned int i, j;
	NODE * n;
	struct edge * e;
	STRING sid;

	for (i=0; i<thiz->all_nodes_num; i++)
	{
		if(!i)
			printf("--------------------------------------\n");
		n = thiz->all_nodes[i];
		printf("NODE(%d)/----FAIL---> NODE(%d)\n", n->id, (n->failure_node)?n->failure_node->id:0);
		for (j=0; j<n->outgoing_degree; j++)
		{
			e = &n->outgoing[j];
			printf("         |----(%c)----> NODE(%d)\n", e->alpha, e->next->id);
		}
		printf("ACCEPTED STRING: {");
		for (j=0; j<n->matched_strings_num; j++)
		{
			sid = n->matched_strings[j];
			printf("%ld ", sid.id);
		}
		printf("}\n");
		printf("--------------------------------------\n");
	}
}
#endif

