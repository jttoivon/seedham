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
#ifndef AHO_CORASICK_WRAPPER_HPP
#define AHO_CORASICK_WRAPPER_HPP

#include "ahocorasick/aho_corasick.h"

#include <string>
#include <vector>

class aho_corasick
{
public:

  aho_corasick(int match_handler(MATCH* m, void* param)) {
    ac_automata_init(&aca, match_handler);
  }

  void
  add_string(const std::string& s, unsigned long id)
  {
    strings.push_back(s);

    STRING tmp_shell;
    tmp_shell.str = strings.back().c_str();
    tmp_shell.id = id;
    tmp_shell.length = s.length();
    
    shells.push_back(tmp_shell);

    ac_automata_add_string (&aca, &tmp_shell);
  }

  /*
  void
  search(const std::string& s) 
  {
    ac_automata_locate_failure(&aca);
    STRING tmp_str;

    tmp_str.str = s.c_str();
    tmp_str.length = s.length();
    ac_automata_search(&aca, &tmp_str, (void *)&id_to_count);
  }
  */

  void
  search(const std::string& s, void* param) 
  {
    ac_automata_locate_failure(&aca);   /// HUOM!! ONKO OK KUTSUA TÄTÄ KAHDESTI
    STRING tmp_str;

    tmp_str.str = s.c_str();
    tmp_str.length = s.length();
    ac_automata_search(&aca, &tmp_str, param);
  }

  ~aho_corasick()
  {
    ac_automata_release(&aca);
  }

private:
  AC_AUTOMATA aca;
  std::vector<std::string> strings;
  std::vector<STRING> shells;   // char*, id, length triples
};

#endif // AHO_CORASICK_WRAPPER_HPP
