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

#ifndef _NODE_H_
#define _NODE_H_

#include "config.h"
#include "ac_types.h"

/* Forward Declaration */
struct edge;

typedef struct node
/* The Node of the Automata */
{
	int id; /* Node ID : for debugging purpose */
	short int final; /* 0: no ; 1: yes, it is a final node */
	struct node * failure_node; /* The failure node of this node */
	unsigned short depth; /* depth: distance between this node to the root */

	/* Matched Strings */
	STRING * matched_strings; /* Array of matched strings */
	unsigned short matched_strings_num; /* Number of matched string at this node */
	unsigned short matched_strings_max; /* Max capacity of allocated memory for 'matched_strings' */

	/* Outgoing Edges */
	struct edge * outgoing; /* Array of outgoing edges */
	unsigned short outgoing_degree; /* Number of outgoing edges */
	unsigned short outgoing_max; /* Max capacity of allocated memory for 'outgoing' */
} NODE;

struct edge
/* The Edge of the Node */
{
	ALPHA alpha; /* Edge alpha */
	struct node * next; /* Target of the edge */
};

/* Public Functions */
NODE * node_create            (void);
NODE * node_create_next       (NODE * thiz, ALPHA alpha);
void   node_register_matchstr (NODE * thiz, STRING * str);
void   node_register_outgoing (NODE * thiz, NODE * next, ALPHA alpha);
NODE * node_find_next         (NODE * thiz, ALPHA alpha);
NODE * node_findbs_next       (NODE * thiz, ALPHA alpha);
void   node_release           (NODE * thiz);
void   node_assign_id         (NODE * thiz);
void   node_sort_edges        (NODE * thiz);

#endif
