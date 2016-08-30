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

#ifndef _AUTOMATA_H_
#define _AUTOMATA_H_

#include "config.h"
#include "node.h"

typedef struct
{
	/* The root of the Aho-Corasick trie */
	NODE * root;

	/* maintain all node pointers in this automata.
	   it will be used to traverse or release all nodes.
	*/
	NODE ** all_nodes;
	unsigned int all_nodes_num; /* Number of all nodes in the automata */
	unsigned int all_nodes_max; /* Max capacity of allocated memory for *all_nodes */

	MATCH match; /* Any match is writen in here */
	MATCH_CALBACK match_callback; /* Match callback function */

	/* this is a flag for string acceptance check 
	   after ac_automata_locate_failure() it is set to 0 to indicate end of strings
	   it is set to 1 in ac_automata_init()
	*/
	unsigned int accept_strings;

	/* 
	   following members keep automata state for sake of repeated call for
	   ac_automata_search(); it needed when we deal with larg input string.
	   in that case we break input string into smaller chunk and repeatedly
	   call ac_automata_search().
	*/
	NODE * current_node; /* Pointer to current node while searching */
	unsigned long base_position; /* Represents the position of current chunk related to whole input */

	/* Statistic Variables */
	unsigned long total_strings; /* Total Strings in the Automata */

} AC_AUTOMATA;


/* Public Functions */
void     ac_automata_init           (AC_AUTOMATA * thiz, MATCH_CALBACK mc);
AC_ERROR ac_automata_add_string     (AC_AUTOMATA * thiz, STRING * str);
void     ac_automata_locate_failure (AC_AUTOMATA * thiz);
void     ac_automata_search         (AC_AUTOMATA * thiz, STRING * str, void * param);
void     ac_automata_reset          (AC_AUTOMATA * thiz);
void     ac_automata_release        (AC_AUTOMATA * thiz);

#ifdef DEBUG_DISPLAY_AC
void     ac_automata_dbg_show       (AC_AUTOMATA * thiz);
#endif


#endif

