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
