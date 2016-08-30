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

#ifndef _AC_TYPES_H_
#define _AC_TYPES_H_

/* ALPHA

   We use fix alpha size in this implementation
   if you are intended to use variable size codings
   you must convert it to a fix one (see iconv library).

   define your alpha type in below.
*/
typedef char ALPHA;

/* String Identifier 

   Every string may have an identifier to distinguish it
   in Aho-Corasik automata. It may be useful while retrieving 
   strings from a database.
   The defined type must be big enough to enumerate all strings.
*/
typedef unsigned long STRINGID;

/* STRING type */
typedef struct
{
	const ALPHA * str; /* Array of ALPHAs */
	unsigned int length; /* Length of string */
	STRINGID id; /* String identifier (Optional) */
} STRING;

/* Maximum Pattern Length */
#define AC_PATTRN_MAX_LENGTH 256

/* Any Match is Reported in the following structure */
typedef struct
{
	/* Array of Matched Strings */
	STRING * matched_strings;
	/* Matched Position:
	   determine the last element (alpha) of string(s) */
	long position;
	/* Number of matched string (length of sid array) */
	unsigned int match_num;
} MATCH;

/* MATCH_CALBACK: 
	Callback function type 
	We use a callback function to report a match accurance
	to the caller.
*/
typedef int (*MATCH_CALBACK)(MATCH *, void *);

/* Error Numbers */
typedef enum
{
	ACERR_NONE = 0,
	ACERR_DUPLICATE_STRING,
	ACERR_LONG_STRING,
	ACERR_ZERO_STRING,
	ACERR_STRING_CLOSED,
} AC_ERROR;

#endif

