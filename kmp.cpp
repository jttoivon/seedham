/*

    SeedHam is a program to learn DNA binding motifs from SELEX datasets.
    Copyright (C) 2016, 2017  Jarkko Toivonen

    SeedHam is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    SeedHam is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/
#include <string>
#include <cstdlib>

void preKmp(const char *x, int m, int kmpNext[]) 
{
   int i, j;

   i = 0;
   j = kmpNext[0] = -1;
   while (i < m) {
      while (j > -1 && x[i] != x[j])
         j = kmpNext[j];
      i++;
      j++;
      if (x[i] == x[j])
         kmpNext[i] = kmpNext[j];
      else
         kmpNext[i] = j;
   }
}



// search function
int
KMP_func(const std::string& t, const std::string& p, bool find_first) 
{
   int i, j;
   int* kmpNext;
   const char* pattern = p.c_str();
   int m = p.length();
   const char* y = t.c_str();
   int n = t.length();

   kmpNext = (int*)malloc((m+1)*sizeof(int));


   int occurence_counter=0;
   /* Preprocessing */
   preKmp(pattern, m, kmpNext);

   /* Searching */
   i = j = 0;
   while (j < n) {
      while (i > -1 && pattern[i] != y[j])
         i = kmpNext[i];
      i++;
      j++;
      if (i >= m) {
	if (find_first) {
	  free(kmpNext);
	  return j - i;
	}
        //printf("%i\n",j - i);
	++occurence_counter;
        i = kmpNext[i];
      }
   }
   free(kmpNext);
   if (find_first)
     return -1;   // no items found
   else
     return occurence_counter;
}


// returns number of occurences
int
KMP(const std::string& t, const std::string& p) 
{
  return KMP_func(t, p, false);
}


// returns position of the first occurence or -1
int
KMP_first_occurence(const std::string& t, const std::string& p) 
{
  return KMP_func(t, p, true);
}
