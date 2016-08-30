#define TIMING

//extern "C" {
//#include "suffix_array/interface.h"
//}


//#include "ahocorasick/aho_corasick.h"
#include "aho_corasick_wrapper.hpp"
#include "suffix_array_wrapper.hpp"

#include "timing.hpp"
#include "my_assert.hpp"
//#include "fvector.hpp"
#include "combinatorics.hpp"
#include "lambda.hpp"

#include "matrix.hpp"
#include "iupac.hpp"
#include "matrix_tools.hpp"
#include "common.hpp"
#include "probabilities.hpp"
#include "data.hpp"
#include "kmp.hpp"
#include "bndm.hpp"
#include "dependence.hpp"
#include "packed_string.hpp"
#include "parameters.hpp"
#include "dinucleotide.hpp"
#include "kmer_tools.hpp"
#include "multinomial_helper.hpp"

#include <cfloat>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/multi_array.hpp>
#include <boost/unordered_map.hpp>

//int character_count;
//int digram_count=0;

bool statistics=false; // print statistics
bool find_dependence_matrix=false;
bool print_neighbourhood=false;


bool automatic_limit_for_palindromic_index = false;
bool use_multimer=false;   // allow multiple occurrences per sequence
bool use_background_correction=false;
bool use_bernoulli_background=false;
bool use_dinucleotide_model=false;
bool use_cell_correction=false;
bool contains_N = false;


prior<double> pseudo_counts;


std::vector<double> background_frequencies(4);
std::vector<double> background_probabilities(4);
matrix<double> background_frequency_matrix(4,4);   // for background noise
matrix<double> background_probability_matrix(4,4); // for background noise


// this accepts multiple occurences per sequence
void
distribution_of_startpositions(const std::vector<std::string>& sequences, const std::string& str)
{
  int L=sequences[0].length();
  std::vector<int> pos_frequencies_forward(L,0);
  std::vector<int> pos_frequencies_backward(L,0);
  std::vector<int> pos_frequencies_palindrome(L,0);
  std::string rev_str = reverse_complement(str);
  int k=str.length();

  int sum=0;
  for (int i=0; i < sequences.size(); ++i) {
    const std::string& line = sequences[i];
    for (int j=0; j < L-k+1; ++j) {
      if (line.substr(j, k) == str) {
	if (line.substr(j, k) == rev_str)         // palindrome occurrence
	  ++pos_frequencies_palindrome[j];
	else                                      // forward occurrence
	  ++pos_frequencies_forward[j];
      }
      else if (line.substr(j, k) == rev_str)      // backward occurrence
	++pos_frequencies_backward[j];
    }
  }
  
  printf("MCP occurs on %i sequences\n", sum);
  printf("Distribution of start positions of string %s:\n", str.c_str()); 
  printf("D");
  for (int i=0; i < L; ++i)
    printf("\t%i", i);
  printf("\n");
  printf("Forw");
  for (int i=0; i < L; ++i)
    printf("\t%i", pos_frequencies_forward[i]);
  printf("\n");
  printf("Back");
  for (int i=0; i < L; ++i)
    printf("\t%i", pos_frequencies_backward[i]);
  printf("\n");
  printf("Palind");
  for (int i=0; i < L; ++i)
    printf("\t%i", pos_frequencies_palindrome[i]);
  printf("\n");

  //  for (int i=0; i < L; ++i)
  //    printf("%.2f ", (float)pos_frequencies[i]/sum);
}



// print occurence statistics for all possible k-mers
void
print_hamming_statistics2(const std::vector<std::string>& sequences, const std::string& str, int delta)
{
  printf("\nStatistics 2\n");
  printf("============\n");
  int k=str.length();
  int lines = sequences.size();

  std::vector<std::set<int> > seqs((int)pow(4,k)); // in which sequences are the k-mer contained
  std::vector<int> all;       // how many times each k-mer occurs

  for (int l=0; l <= delta; ++l) {
    all.assign((int)pow(4,k), 0);
    for  (int i=0; i < lines; ++i) {
      const std::string& line = sequences[i];
      for (int j=0; j < line.length()-k+1; ++j) {
	int d = hamming_distance(str, line.substr(j,k));
	if (d==l) {          // found occurence
	  ++all[dna_to_number(line.substr(j,k))];
	  seqs[dna_to_number(line.substr(j,k))].insert(i);
	}

      }
    }
    printf("\nHamming distance %i\n", l);
    int count=0;    // the number of unique k-mers contained in data
    int count2=0;   // the number of k-mers that exactly only once per sequence
    double ratio=0;
    for (int i=0; i < (int)pow(4,k); ++i) {
      if (all[i]==0)
	continue;
      std::string s = number_to_dna(i,k);
      assert(all[i]>=0);
      printf("String %s occurs %4i times on %4zu sequences\n", s.c_str(),
	     all[i], seqs[i].size());
      ++count;

      if (seqs[i].size() == all[i])  // occurs only once per sequence
	++count2;

      ratio+=(double)all[i]/seqs[i].size();
      seqs[i].clear();
    }
    printf("On average, %lg occurences per sequence\n", ratio/count);
    printf("%.3lg%% of the above strings occur only once per sequence\n",
	   (double)count2/count*100);
  }
}

void
find_count_of_hamming_neighbourhood(const std::vector<std::string>& sequences, const std::string& seed, int delta)
{
  int k = seed.length();
  std::map<std::string, int> counts;
  int lines = sequences.size();
  for  (int i=0; i < lines; ++i) {
    const std::string& line = sequences[i];
    for (int j=0; j < line.length()-k+1; ++j) {
      std::string s = line.substr(j,k);
      int d = hamming_distance(seed, s);
      if (d <= delta) {
	++counts[s];
      }
    }
  }
  std::string s;
  int count;
  BOOST_FOREACH(boost::tie(s, count), counts) {
    printf("%s\t%i\n", s.c_str(), count);
  }
}

// print occurence statistics for all Hamming neighbourhoods H(str, d) for 0<=d<=delta
void
print_hamming_statistics(const std::vector<std::string>& sequences, const std::string& str, int delta)
{
  printf("\nStatistics 1\n");
  printf("============\n");
  int k=str.length();
  int lines = sequences.size();
  std::vector<int> occs(lines, 0);


  for (int l=0; l <= delta; ++l) {
    int all=0;
    int seqs=0;
    occs.assign(lines, 0);
    for  (int i=0; i < lines; ++i) {
      const std::string& line = sequences[i];
      for (int j=0; j < line.length()-k+1; ++j) {
	int d = hamming_distance(str, line.substr(j,k));
	if (d==l) {
	  ++all;
	  if (occs[i] == 0)
	    ++seqs;
	  ++occs[i];
	}

      }
    }
    printf("\n%i occurences of the Hamming neighbourhood H(str,%i) on %i sequences\n", 
	   all, l, seqs);
    int count2 = std::count(occs.begin(), occs.end(), 1);
//     for (int i=0; i < lines; ++i)
//       if (occs[i]==1)
// 	++count2;
    printf("On average, %lg occurences per sequence\n", (double)all/seqs);
    printf("%.3lg%% of sequences have exactly one occurence\n",
	   (double)count2/seqs*100);
  }
}






double
l2_norm(const std::vector<double>& v)
{
  double s = 0;
  
  for (int i=0; i < v.size(); ++i)
    s += v[i]*v[i];

  return sqrt(s);
}

struct quad
{
  quad(int p, int c, int d) : pos(p), column(c), dir(d) {}
  quad() : pos(-1), column(0), dir(0) {} 
  int pos;
  int column;  // the related PWM column, -1 for the consensus/seed sequence
  int dir;
};

std::string
substring(const std::string& s, int pos, int len, int dir)
{
  std::string t = s.substr(pos,len);
  return dir == 1 ? t : reverse_complement(t);
}


void
find_snips_monomer_helper(std::vector<quad>& positions, const std::string& query, int j, 
			  const std::vector<std::string>& sequences)
{
  int k = query.length();
  int lines = sequences.size();
  for (int i=0; i < lines; ++i) {
    const std::string& line = sequences[i];
    int c1 = BNDM_with_joker(line, query);
    int c2 = use_two_strands ? BNDM_with_joker(line, reverse_complement(query)) : 0;
    if (c1 + c2 == 1 or (c1 == 1 && c2 == 1 && is_palindromic(query))) {
      std::string s = query;
      int dir = (c1 == 0 ? -1 : 1);
      if (dir == -1)
	s=reverse_complement(s);
      size_t first_pos = std::search(line.begin(), line.end(), s.begin(), s.end(), iupac_match) - line.begin();
      assert(iupac_string_match(line.substr(first_pos,k), s));

      switch (positions[i].pos) {
      case -2:                                  // rejected
	continue;
      case -1: positions[i]=quad(first_pos, j, dir);   // was untouched
	continue;
      default:                                           // another occurrence on the same line
	if (positions[i].pos != first_pos)
	  positions[i].pos = -2;  // reject
      }

    }
    else if (c1+c2 > 1)   // definitely reject
      positions[i].pos = -2;
  } // for lines
}

// Number of subsequences used to build the pwm
int
get_total_multinomial_count(const dmatrix& count_pwm, const std::string& seed)
{
  int total_count = 0;
  bool seed_count_added = false;
  const int k = count_pwm.get_columns();
  std::string nucs = "ACGT";
  for (int j=0; j < k; ++j) {      // iterate through all string positions
    for (int a=0; a < 4; ++a) {    // iterate through all characters
      if (not seed_count_added and nucs[a] == seed[j]) {
	total_count += count_pwm(a, j);
	seed_count_added = true;
      }
      else if (not iupac_match(nucs[a], seed[j]))
	total_count += count_pwm(a, j);   // one point mutation counts
    }
  }
  return total_count;
}

// monomeric version
// Reject sequences with multiple query occurences
// Compute the counts for the multinomial1 matrix
matrix<double>
find_snips_monomer(const std::string& seed, const std::vector<std::string>& sequences, int hamming_distance)
{
  assert(hamming_distance == 1);
  TIME_START(t);
  int lines = sequences.size();
  std::vector<quad> positions(lines);

  int k = seed.length();
  char nucs[] = "ACGT";
  printf("Multinomial-1 is using seed %s\n", seed.c_str());

  matrix<double> result(4,k);
  find_snips_monomer_helper(positions, seed, -1, sequences);

  for (int j=0; j < k; ++j) {        // iterate through all string positions
    std::string query = seed;
    for (int a=0; a < 4; ++a) {      // iterate through all characters
      if (nucs[a] == seed[j])  // the consensus count was already computed
	continue;
      query[j]=nucs[a];
      find_snips_monomer_helper(positions, query, j, sequences);
    }
  }
  
  // extract subsequences from reads that have only one hit
  std::vector<std::string> sites;
  for (int i=0; i < lines; ++i) {
    int pos = positions[i].pos;
    if (pos >= 0) {
      std::string s = sequences[i].substr(pos, k);
      sites.push_back(s);
      //      printf("##%s\n", s.c_str());
    }
    // if (pos == -1)
    //   printf("#%s\n", sequences[i].c_str());
    // if (pos == -2)
    //   printf("*%s\n", sequences[i].c_str());
  }

  printf("Number of sites is %zu\n", sites.size());
  int seed_count=0;
  for (int i=0; i < sites.size(); ++i) {
    bool is_palindrome = is_palindromic(sites[i]);
    if (iupac_string_match(sites[i], seed))
      ++seed_count;
    if (use_two_strands and iupac_string_match(reverse_complement(sites[i]), seed) and (count_palindromes_twice or not is_palindrome))
      ++seed_count;
  }

  for (int j=0; j < k; ++j) {        // iterate through all string positions
    std::string query = seed;
    for (int a=0; a < 4; ++a) {      // iterate through all characters
      if (nucs[a] == seed[j])  // the consensus count was already computed
	continue;
      query[j]=nucs[a];
      for (int i=0; i < sites.size(); ++i) {
	bool is_palindrome = is_palindromic(sites[i]);
	if (iupac_string_match(sites[i], query))
	  ++result(a, j);
	if (use_two_strands and iupac_string_match(reverse_complement(sites[i]), query) and (count_palindromes_twice or not is_palindrome))
	  ++result(a, j);
      }
    }
  }

  printf("Seed %s count = %i\n", seed.c_str(), seed_count);

  // add seed counts to the pwm
  for (int j=0; j < k; ++j) {      // iterate through all string positions
    for (int a=0; a < 4; ++a) {    // iterate through all characters
      if (nucs[a] == seed[j]) 
	result(a, j) = seed_count;
    }
  }

  int total_count = get_total_multinomial_count(result, seed);
  printf("Total multinomial-1 count is %d\n", total_count);

  int match_count=0;
  int nomatch_count=0;
  int too_many_matches_count=0;
  for (int i=0; i < lines; ++i)
    switch (positions[i].pos) {
    case -1:
      ++nomatch_count;
      break;
    case -2:
      ++too_many_matches_count;
      break;
    default:
      ++match_count;
    }

  printf("%i sequences matched with hamming distance 1, rate = %g\n",
	 match_count, (double)match_count/lines);

  printf("%i sequences didn't match with hamming distance 1, error rate = %g\n",
	 nomatch_count, (double)nomatch_count/lines);

  printf("%i sequences had too many matches with hamming distance 1, rate = %g\n",
	 too_many_matches_count, (double)too_many_matches_count/lines);

  // print the alignment to file descriptor 3, if it is open
//   if (print_alignment) {
//     FILE* fp = fdopen(3, "a");
//     if (fp != NULL) {
//       for (int t = 0; t < alignment.size(); ++t)
// 	fprintf(fp, "%s\n", alignment[t].c_str());
//       fclose(fp);
//     }
//   }

  TIME_PRINT("Multinomial-1 algorithm took %.2f seconds.\n", t);

  return result;
}  // find_snips_monomer




// helper function for find_multinomial2
int 
find_base_counts(const std::string& base, int fixed_pos, 
		 const std::string& str1, const std::string& str2)
{

  int k = base.length();
  int count=0;
  std::string nucs = "ACGT";

  // this is either consensus or its one point mutation
  count += BNDM_with_joker(str1, base);  
  if (use_two_strands)
    count += BNDM_with_joker(str2, base);


  // these are one or two point mutations of the consensus
  for (int i=0; i < k; ++i) {        // iterate through all string positions, except fixed_pos
    if (i == fixed_pos)
      continue;

    for (int a=0; a<4; ++a) {
      std::string temp = base;
      //      if (to_int(base[i]) == a)      // these were already counted
      if (iupac_match(nucs[a], base[i]))      // these were already counted
	continue;
      temp[i]=nucs[a];
      //temp[i]='.';                     // dot matches any nucleotide
      count += BNDM_with_joker(str1, temp);
      if (use_two_strands)
	count += BNDM_with_joker(str2, temp);
    }
  }

  return count;
}

// this computes the multinomial2 matrix counts
matrix<double>
find_multinomial2(const std::string& seed, const std::vector<std::string>& sequences, int hamming_distance)
{
  assert(hamming_distance == 2);
  TIME_START(t);
  int k = seed.length();
  char nucs[] = "ACGT";
  std::string str1;
  std::string str2;

  str1=join(sequences, '#');
  str2=join_rev(sequences, '#');

  matrix<double> result(4,k);
  for (int i=0; i < k; ++i) {        // iterate through all string positions
    std::string temp = seed;
    for (int j=0; j < 4; ++j) {      // iterate through all characters
      temp[i]=nucs[j];
      result(j, i) = find_base_counts(temp, i, str1, str2);
    }

  }
  printf("Multinomial-2 algorithm took\n");
  TIME_CHECK(t);
  printf("seconds\n");
  return result;
}

typedef boost::unordered_map<int, int> id_to_count_type;
//typedef std::vector<int> id_to_count_type;
typedef boost::unordered_map<big_int, std::vector<boost::tuple<int, int> > > code_to_tuple_type;
typedef boost::unordered_map<std::string, std::vector<boost::tuple<int, int> > > string_to_tuple_type;

// callback function for Aho-Corasick automaton
int match_handler_counts(MATCH* m, void* param)
{

  //printf ("ending at @ position %ld string(s) ", m->position-4);

  assert(m->match_num == 1);
  big_int id = m->matched_strings[0].id;  // This is ordinal id, not coded DNA sequence

  id_to_count_type* result = (id_to_count_type*)param;
  (*result)[id] += 1;

  /* to find all matches always return 0 */
  return 0;

}

// Callback function for Aho-Corasick automaton.
// Add all the occurrences positions to the list of the corresponding pattern
int 
match_handler_positions(MATCH* m, void* param)
{

  assert(m->match_num == 1);   // Could be larger than zero, if one pattern is a suffix of another
  big_int id = m->matched_strings[0].id;  // This is ordinal id, not coded DNA sequence

  std::vector<std::vector<long int> >* result = (std::vector<std::vector<long int> >*)param;
  (*result)[id].push_back(m->position);

  /* to find all matches always return 0 */
  return 0;

}

// this computes the multinomial2 matrix counts using AC-automaton
matrix<double>
find_multinomial2_v2(const std::vector<std::string>& sequences, const std::string& consensus)
{
  TIME_START(t);
  int k = consensus.length();
  char nucs[] = "ACGT";
  std::string str1;
  std::string str2;

  str1=join(sequences, '#');
  str2=join_rev(sequences, '#');
  matrix<double> result(4,k);
  big_int myid = dna_to_number(consensus);

  code_to_tuple_type code_to_tuple;

  for (int j=0; j < k; ++j) {        // iterate through all string positions
    std::string temp = consensus;
    int shift1 = 2*(k-j-1);
    big_int myid2 = myid - ((big_int)to_int(consensus[j]) << shift1);
    for (int a=0; a < 4; ++a) {      // iterate through all characters
      temp[j]=nucs[a];
      big_int myid3 =  myid2 + ((big_int)a << shift1);
      //my_assert(myid3, dna_to_number(temp));
      code_to_tuple[myid3].push_back(boost::make_tuple(j, a));
      

      for (int h=0; h < k; ++h) {        // iterate through all string positions, except fixed_pos
	if (h == j)
	  continue;
	std::string temp2 = temp;
	int shift2 = 2*(k-h-1);
	big_int myid4 = myid3 - ((big_int)to_int(consensus[h]) << shift2);
	for (int b=0; b<4; ++b) {
	  if (to_int(temp[h]) == b)
	    continue;
	  big_int myid5 = myid4 + ((big_int)b << shift2);
	  temp2[h]=nucs[b];
	  //my_assert(myid5, dna_to_number(temp2));

	  code_to_tuple[myid5].push_back(boost::make_tuple(j, a));

	}
      }
    }
  }

  std::vector<std::string> strings;
  std::vector<big_int> id_to_code;
  big_int code;
  BOOST_FOREACH(boost::tie(code, boost::tuples::ignore), code_to_tuple) {
    strings.push_back(number_to_dna(code, k));
    id_to_code.push_back(code);
  }

  aho_corasick my_aca(match_handler_counts);
  for (int i=0; i < strings.size(); ++i) {
    my_aca.add_string(strings[i], i);
  }
  big_int seed_id = 0;

  id_to_count_type id_to_count;

  my_aca.search(str1, &id_to_count);
  if (use_two_strands) {
    my_aca.search(str2, &id_to_count);
  }
    
  printf("Seed %s count = %i\n", consensus.c_str(), id_to_count[seed_id]);
  int id;
  int count;
  BOOST_FOREACH(boost::tie(id, count), id_to_count) {
    std::vector<boost::tuple<int,int> >& temp = code_to_tuple[id_to_code[id]];
    for (int i=0; i < temp.size(); ++i) {
      int j, a;
      boost::tie(j,a) = temp[i];
      result(a, j) += count;
    }
  }

  printf("Multinomial-2 algorithm took\n");
  TIME_CHECK(t);
  printf("seconds\n");
  return result;
} // find_multinomial2_v2


int 
myskip(int i, int skip)
{
  assert(i >= 0);
  assert(i <= 2);
  assert(skip >= 0);
  assert(skip < 4);

  return i >= skip ? i+1 : i;
}


std::string
bit_representation_helper(void* f, int bytes)
{
  
  int number_of_bits = 8*bytes;
  std::string bit_repr(number_of_bits, '-');
  int pos=0;
  unsigned char* start = (unsigned char*)f;
  //  unsigned char* start = (unsigned char*)&f;
  for (int i = bytes-1; i >= 0; --i) {
    unsigned char c = start[i];
    for (int j = 7; j >= 0; --j) {
      unsigned char mask = 1 << j;
      bit_repr[pos] = c & mask ? '1' : '0';
      ++pos;
    }
  }
  return std::string(bit_repr);
}


// this computes the multinomial-n matrix counts by scanning all possible windows
dmatrix
find_multinomial_n_scan(const std::string& seed, const std::vector<std::string>& sequences, int n)
{
  TIME_START(t);
  const int k = seed.length();
  const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  //  char nucs[] = "ACGT";
  dmatrix result(4, k);
  for (int i=0; i < sequences.size(); ++i) {
    int max_dir = use_two_strands ? 2 : 1;
    for (int dir=0; dir < max_dir; ++dir) {
      const std::string& line = dir == 0 ? sequences[i] : reverse_complement(sequences[i]);
      for (int j=0; j < L-k+1; ++j) {
	std::string query = line.substr(j, k);
	bool is_palindrome = is_palindromic(query);
	if (dir==1 and is_palindrome and not count_palindromes_twice)
	  continue;
	int hd = iupac_hamming_dist(query, seed, n);
	if (hd > n)
	  continue;
	for (int pos=0; pos < k; ++pos) {
	  char c = line[j+pos];
	  int cc = iupac_match(c, seed[pos]) ? 0 : 1;   // 1 if mismatch
	  if (hd-cc <= n-1)
	    ++result(to_int(c), pos);
	}
      }
    }
  }

  TIME_PRINT("Multinomial-n scanning algorithm took %.2f seconds\n", t);
  return result;
}



// this computes the multinomial-n matrix counts using suffix array (was AC-automaton)
dmatrix
find_multinomial_n(const std::string& seed, const std::vector<std::string>& sequences, int n)
{
  TIME_START(t);
  const int k = seed.length();
  //const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  //char nucs[] = "ACGT";
  std::string str1;
  std::string str2;

  str1=join(sequences, '#');
  if (use_two_strands) {
    str1.append("#");
    str1 += join_rev(sequences, '#');
  }

  dmatrix result(4, k);

  code_to_tuple_type code_to_tuple;
  string_to_tuple_type string_to_tuple;

  bool use_suffix_array=true;

  if (use_suffix_array) {
     suffix_array sa(str1);
     result = find_multinomial_n_suffix_array(seed, sequences, sa, n, use_multimer).get<0>();
  }
  else {  // use Aho-Corasick automaton
    std::vector<std::pair<std::string, std::vector<boost::tuple<int, int> > > >  patterns;
    patterns = get_n_neighbourhood_in_vector(seed, n);
    printf("There are %zu distinct sequences in the Hamming-%i neighbourhood\n", patterns.size(), n);

    //    aho_corasick my_aca(match_handler_positions);
    aho_corasick my_aca(match_handler_counts);
    for (int i=0; i < patterns.size(); ++i) {
      my_aca.add_string(patterns[i].first, i);
    }
    
    id_to_count_type ordinal_id_to_count(patterns.size());
    std::vector<std::vector<long int> > positions;
    //    my_aca.search(str1, &positions);
    my_aca.search(str1, &ordinal_id_to_count);
    
    printf("Seed %s count = %i\n", seed.c_str(), ordinal_id_to_count[0]);
    int count;
    int ordinal_id;
    BOOST_FOREACH(boost::tie(ordinal_id, count), ordinal_id_to_count) {
      std::vector<boost::tuple<int,int> >& temp = patterns[ordinal_id].second;
      std::string query = patterns[ordinal_id].first;
      bool is_palindrome = is_palindromic(query);
      if (is_palindrome and not count_palindromes_twice)
	count /= 2;
      int j, a;
      //      if (use_palindromic_correction)
      //	count *= palindromic_correction(pattern, seed, seed_rev);

      BOOST_FOREACH(boost::tie(j,a), temp) 
	result(a, j) += count;
    }
  }

  printf("Multinomial-n algorithm took\n");
  TIME_CHECK(t);
  printf("seconds\n");
  return result;
} // find_multinomial_n


// this computes the init_motif matrix counts
// Chooses the first occurence with lowest Hamming distance for each sequence
matrix<double>
find_hamming_neighbourhood(const std::vector<std::string>& sequences,
			   const std::string& str, int delta)
{
  int k=str.length();
  int lines = sequences.size();
  std::vector<int> alignments(lines, -1);
  std::vector<int> direction(lines, 0);


  // brute force scan of all k-grams
  for (int l=0; l <= delta; ++l) {
    for  (int i=0; i < lines; ++i) {
      if (alignments[i] != -1)
	continue;
      const std::string& line = sequences[i];
      for (int j=0; j < line.length()-k+1; ++j) {
	std::string substr = line.substr(j,k);
	if (hamming_distance(str, substr)==l) {
	  alignments[i] = j;
	  direction[i] = 1;
	  break;
	}
	else if (use_two_strands && 
		 hamming_distance(str, reverse_complement(substr))==l) 
	  {
	  alignments[i] = j;
	  direction[i] = -1;
	  break;	  
	  }
      }
    }
  }



  int count = std::count(alignments.begin(), alignments.end(), -1);
  //assert(count == count2);
  printf("%i sequences didn't match with hamming distance %i, error rate = %g\n",
	 count, delta, (double)count/lines);

  std::vector<int> indexes(lines-count,-1);
  
  // create subset of matching sequences as set of indexes
  int n=0;
  for (int i=0; i < lines; ++i) 
    if (alignments[i]!=-1) {
      indexes[n]=i;
      ++n;
    }
  assert(n==indexes.size());

  // print the alignment to file descriptor 3, if it is open
  if (print_neighbourhood) {
    FILE* fp = fdopen(3, "a");
    if (fp != NULL) {
      for (int i=0; i < indexes.size(); ++i) {
	int j = alignments[indexes[i]];
	int d = direction[indexes[i]];
	std::string str = sequences[indexes[i]].substr(j, k);
	if (d == -1)
	  str=reverse_complement(str);
	fprintf(fp, "%s\n", str.c_str());
      }
      fclose(fp);
    }
  }

  //if (print_neighbourhood) {
  //  printf("%s, Hamming neighbourhood %d\n", header.c_str(), delta);
  //  for (int i=0; i < indexes.size(); ++i)
  //    printf("%s\n", sequences[indexes[i]].c_str());
  //}

  matrix<double> motif = 
    calculate_frequencies_in_subset(sequences, indexes.begin(), indexes.end(), 
				    alignments, direction, k); 


  return motif;

}





dmatrix
make_consensus_matrix(const std::string& consensus, int count)
{
  int k = consensus.length();
  dmatrix result(4, k);
  for (int i=0;i<k;++i)
    result.set_column(i, iupac_probability(consensus[i]));
  //    result(to_int(consensus[i]),i)=count;
  return result;
}


// return an interval of length k with the highest sum of elements
boost::tuple<double, int>
best_interval(const std::vector<double>& v, int k)
{
  double max = -DBL_MAX;
  int argmax = 0;
  for (int i=0; i < v.size()-k; ++i) {
    double sum=0;
    for (int j=0; j<k; ++j)
      sum += v[i+j];
    if (sum > max) {
      max = sum;
      argmax = i;
    }
  }
  return boost::make_tuple(max, argmax);
}

// print subinterval of lengths 3...k with maximum information content
void
print_submotifs(const dmatrix& used_motif)
{
  const int min_len=3;
  int k = used_motif.get_columns();
  std::vector<double> ic(k);

  // starting position of maximal subinterval for each length
  std::vector<int> start(k-min_len+1);
  std::vector<double> bg(4, 0.25);   // even background distribution

  for (int i=0; i < k; ++i) {
    ic[i] = information_content(used_motif.column(i), bg);
    //printf("ic[%i]=%f\n", i, ic[i]);  // for debug
  }

  int current_pos;
  double current_sum;
  // find the position of the shortest interval
  boost::tie(current_sum, current_pos) = best_interval(ic, min_len);
  start[0]=current_pos;

  for (int len=min_len+1; len <= k; ++len) { // iterate over the rest of subinterval lengths
    //printf("Current pos: %i Current sum: %f\n", current_pos, current_sum);  // for debug
    double sum1=-DBL_MAX;
    double sum2=-DBL_MAX;
    if (current_pos > 0) {            // extend interval to left
      sum1=current_sum+ic[current_pos-1];
      // printf("sum1: %f\n", sum1);  // debug
    }
    if (current_pos+len-1 < k) {      // extend interval to right
      sum2=current_sum+ic[current_pos+len-1];
      // printf("sum2: %f\n", sum2); // debug
    }
    if (sum1 > sum2) {
      current_sum = sum1;
      --current_pos;
    } else {
      current_sum = sum2;
    }
    start[len-min_len]=current_pos;
  }

  printf("Maximum IC subintervals:\n");
  for (int i=0; i < start.size(); ++i)
    printf("%i,%i ", start[i], i+min_len);
  printf("\n");
  return;
}

void
reverse_complement_sequences(std::vector<std::string>& sequences)
{
  int lines = sequences.size();
  for (int i=0; i < lines; ++i)
    sequences[i] = reverse_complement(sequences[i]);
}

void
print_ics(const dmatrix& m)
{
  std::vector<double> bg(4, 0.25);   // even background distribution
  int k = m.get_columns();

  std::vector<double> ic(k);
  for (int i=0; i < k; ++i) 
    ic[i] = information_content(m.column(i), bg);
  printf("Information content by columns\n");
  printf("%s\n", print_vector(ic, "\t", 2).c_str());
  printf("Average information content is %.2f\n", sum(ic) / k);
}

typedef dmatrix (*func_ptr_t)(const std::string&, const std::vector<std::string>&, int);


dmatrix
find_model(func_ptr_t func_ptr, const std::string& name, const std::string& seed, int hamming_distance,
	   std::vector<std::string>& sequences, std::vector<std::string>& sequences_bg, const std::vector<double>& bg,
	   double lambda)
{
  //int k = seed.length();
  int lines = sequences.size();
  int lines_bg = sequences_bg.size();
  //int L = sequences[0].length();
  dmatrix motif;
  motif = func_ptr(seed, sequences, hamming_distance);
  
  if (use_background_correction) {
    printf("\n");
    write_matrix(stdout, motif, to_string("Signal %s counts before correction:\n", name.c_str()), "%.0f");

    dmatrix motif_bg;
    motif_bg = func_ptr(seed, sequences_bg, hamming_distance);
    write_matrix(stdout, motif_bg, to_string("Background %s counts before correction:\n", name.c_str()), "%.0f");

    dmatrix m_new;
    if (use_cell_correction)
      m_new = (double)lines_bg/lines*motif - lambda*motif_bg;
    else
      m_new = motif - ((double)lines/lines_bg*lambda)*motif_bg;

    m_new.apply(cut);
    motif = m_new;
    printf("\n");
  } else if (use_bernoulli_background) {
    assert(use_multimer);
    write_matrix(stdout, motif, to_string("Signal %s counts before correction:\n", name.c_str()), "%.0f");
    dmatrix mean_matrix;
    int total_count;
    boost::tie(mean_matrix, total_count) =
      find_multinomial_n_background(seed, sequences, bg,
				    hamming_distance, use_multimer);
    write_matrix(stdout, mean_matrix, "Background matrix:\n", "%.0f");
    dmatrix m_new;
    m_new = motif - mean_matrix;
 
    m_new.apply(cut);
    motif = m_new;
    printf("\n");
  }

  write_matrix(stdout, motif, to_string("%s motif matrix counts:\n", name.c_str()), "%.0f");
  dmatrix norm=motif;
  normalize_matrix_columns(norm);
  write_matrix(stdout, norm, to_string("%s motif matrix:\n", name.c_str()), "%.6f"); 
  print_ics(norm);
  printf("\n");

  return motif;
}


dinuc_model
find_dinucleotide_model(func_ptr_t func_ptr, const std::string& name, const std::string& seed, int hamming_distance,
			std::vector<std::string>& sequences, std::vector<std::string>& sequences_bg,
			double lambda)
{
  int k = seed.length();
  int lines = sequences.size();
  int lines_bg = sequences_bg.size();
  int L = sequences[0].length();
  dmatrix motif = motif = func_ptr(seed, sequences, hamming_distance);

  if (use_background_correction) {
    printf("\n");
    write_matrix(stdout, motif, to_string("Signal %s counts before correction:\n", name.c_str()), "%.0f");

    dmatrix motif_bg;
    motif_bg = func_ptr(seed, sequences_bg, hamming_distance);
    write_matrix(stdout, motif_bg, to_string("Background %s counts before correction:\n", name.c_str()), "%.0f");

    dmatrix m_new;
    if (use_cell_correction)
      m_new = (double)lines_bg/lines*motif - lambda*motif_bg;
    else
      m_new = motif - ((double)lines/lines_bg*lambda)*motif_bg;

    m_new.apply(cut);
    motif = m_new;
    printf("\n");
  } else if (use_bernoulli_background) {
    // BERNOULLI CORRECTION ISN'T CURRENTLY CORRECT SINCE IT DOESN'T CONSIDER N's OR OTHER IUPACs, OR HD>2
    assert(use_multimer);
    write_matrix(stdout, motif, to_string("Signal %s counts before correction:\n", name.c_str()), "%.0f");
    double mean = lines * (use_two_strands ? 2 : 1) * (L-k+1) * pow(4, -k);
    
    dmatrix mean_matrix(4, k);
    mean_matrix.fill_with(mean);
    dmatrix m_new;
    m_new = motif - mean_matrix;
 
    m_new.apply(cut);
    motif = m_new;
    printf("\n");
  }

  write_matrix(stdout, motif, to_string("%s dinucleotide motif matrix counts:\n", name.c_str()), "%.0f");
  dinuc_model dm;
  dm.init(motif);

  return dm;
}

int main(int argc, char* argv[])
{

  print_command_line(argc, argv);

  //setvbuf(stdout, NULL, _IOLBF, 1024);
  using namespace boost;
  using std::string;
  using std::cout;

  std::vector<std::string> sequences;
  bool use_reverse_strand=false;
  use_two_strands = true;
  int hamming_radius=1;
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  // 
  // Do the initialization
  //
  ///////////////////////////////////////////////////////////////////////////////////////////

  // Declare the supported options.
  po::options_description nonhidden("Allowed options");
  nonhidden.add_options()
    ("help", "produce help message")
    ("hamming-radius", po::value<int>(), m("Maximum Hamming radius", 
					     hamming_radius).c_str())
    ("statistics", "Print matching statistics")
    ("single-strand", m("Assume sequences can come from either strand", use_two_strands).c_str())
    ("reverse-strand", m("Assume sequences come from reverse strand", false).c_str())
    ("count-palindromes-twice", m("Count palindromes twice", count_palindromes_twice).c_str())
    ("dinucleotide-model", m("Compute dinucleotide model", use_dinucleotide_model).c_str())
    ("use-cell-correction", m("Use different formula to do background correction", use_cell_correction).c_str())
    ("use-bernoulli-background", m("Use bernoulli model to subtract background", use_bernoulli_background).c_str())
    ("pseudo-counts", "Use pseudo counts, default: no")
    // ("prior", po::value<std::string>(), 
    //  "Choose either addone or dirichlet prior")
    ("background", po::value<std::string>(),      
     "Filename of the background sequences, enables also background correction")
    ("output-matrix", po::value<std::string>(), "Name of the output matrixfile, default: none")
    ("multimer", m("Use multimer version for multinomial-1/2", use_multimer).c_str())
    //    ("palindromic-correction", m("Correct mixing of direction by Esko's method", use_palindromic_correction).c_str())
    ("palindromic-index",  po::value<std::string>(), m("Require that seed has high-enough palindromic index: automatic (depending on Hamming radius), 0, 1, ...", palindromic_index_limit).c_str())
    ("dependence-matrix", m("Find dependence matrix of result", find_dependence_matrix).c_str())
    ("print-neighbourhood", m("Print the Hamming neighbourhood of consensus", 
			      print_neighbourhood).c_str())
    ("print-alignment", "Print the alignment, default: no")
    ;

  po::options_description hidden("Positional parameters");
  hidden.add_options()
    ("seed", po::value<std::string>(), "Motif length or seed")
    ("seqs", "Name of the fasta file containing sequences")
    ;

  po::options_description desc;
  desc.add(hidden).add(nonhidden);
  
  po::positional_options_description p;
  p.add("seed", 1);
  p.add("seqs", 1);

  string seqsfile;
  string background_file;
  string seed_param;
  string prior_parameter;
  string matrixfile = "-";
  int k;

  try {
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
	      options(desc).positional(p).run(), vm);
    po::notify(vm);

    // are all positional parameters present
    bool all_positional = vm.count("seqs")&&
      vm.count("seed");

    if (vm.count("help") || not all_positional) {
      cout << "usage: " << argv[0] 
	   << " [options] [ motif_width | seed ] fastafile \n";
      cout << nonhidden << "\n";
      return vm.count("help") ? 0 : 1;
    }

    if (vm.count("statistics")) {
      statistics=true;
    }   

    if (vm.count("single-strand"))
      use_two_strands = false;

    if (vm.count("count-palindromes-twice"))
      count_palindromes_twice = true;

    if (vm.count("reverse-strand")) {
      use_two_strands = false;
      use_reverse_strand = true;
    }

    if (vm.count("use-cell-correction"))
      use_cell_correction = true;

    if (vm.count("dinucleotide-model"))
      use_dinucleotide_model = true;

    if (vm.count("pseudo-counts"))
      use_pseudo_counts = true;


    if (vm.count("hamming-radius"))
      hamming_radius = vm["hamming-radius"].as< int >();

    if (vm.count("palindromic-index")) {
      std::string param = vm["palindromic-index"].as<std::string>();
      if (param == "automatic") {
	automatic_limit_for_palindromic_index = true;
	palindromic_index_limit = conflict_free_palindromic_index(hamming_radius);
      }
      else {
	int i = atoi(param);
	assert(i >= 0);
	palindromic_index_limit = i;
      }
    }
    
    //    if (vm.count("palindromic-correction"))
    //      use_palindromic_correction = true;

    if (vm.count("multimer"))
      use_multimer = true;

    if (vm.count("use-bernoulli-background")) {
      assert(vm.count("background") == 0);
      use_bernoulli_background = true;
    }

    if (vm.count("dependence-matrix"))
      find_dependence_matrix = true;

    if (vm.count("print-neighbourhood"))
      print_neighbourhood = true;

    if (vm.count("print-alignment"))
      print_alignment = true;


    
    
    // if (vm.count("prior")) {
    //   prior_parameter = vm["prior"].as< string >();
    //   if (prior_parameter == "addone" || prior_parameter == "dirichlet")
    // 	use_pseudo_counts=true;
    //   else {
    // 	fprintf(stderr, "Invalid prior: %s\n", prior_parameter.c_str());
    // 	exit(1);
    //   }
    // }   
    
    if (vm.count("background")) {
      background_file = vm["background"].as< string >();
      use_background_correction=true;
    }   
    

    if (vm.count("output-matrix")) 
      matrixfile = vm["output-matrix"].as< string >();
    
    seqsfile = vm["seqs"].as< string >();

    k = atoi(vm["seed"].as< std::string >().c_str());
    if (k == 0) {                    // a seed was given instead of motif width
      seed_param = vm["seed"].as< string >();
      boost::to_upper(seed_param);
      assert(is_iupac_string(seed_param));
      k = seed_param.length();
    }

  }
  catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    std::cerr << desc << "\n";
    exit(1);
  }

  //  initialize_to_int();


  int lines, bad_lines;


  boost::tie(lines, bad_lines) = read_sequences(seqsfile, sequences, true);
  



  check_data(sequences, "ACGTN");
  assert(sequences.size() == lines);
  if (use_reverse_strand)
    reverse_complement_sequences(sequences);



  int L = sequences[0].length();
  printf("Read %zu good lines from file %s\n", 
	  sequences.size(), seqsfile.c_str());
  printf("Input contained %d bad lines\n", bad_lines);
  printf("Sequence length is %i\n", L);
  printf("Motif length is %i\n", k);
  printf("Max hamming distance %i\n", hamming_radius);
  printf("Use two dna strands: %s\n", use_two_strands ? "yes" : "no");
  printf("Count palindromes twice: %s\n", count_palindromes_twice ? "yes" : "no");
  printf("Use reverse strand: %s\n", use_reverse_strand ? "yes" : "no");
  printf("Use bernoulli background: %s\n", use_bernoulli_background ? "yes" : "no");
  //  printf("Use palindromic correction: %s\n", use_palindromic_correction ? "yes" : "no");
  printf("Palindromic index limit is %i\n", palindromic_index_limit);
  //  printf("Require directional seed: %s\n", require_directional_seed ? "yes" : "no");
#ifdef USE_HASH
  printf("Use hashing: yes\n");
#else
  printf("Use hashing: no\n");
#endif
  //printf("Sequences contain %i characters\n", character_count);

  std::vector<double> background_frequencies;
  std::vector<int> character_frequencies;
  boost::tie(background_frequencies, background_frequency_matrix, character_frequencies) = count_background(sequences);
  if (character_frequencies['N'] > 0) {
    contains_N = true;
  }
  printf("Data contains non-specific nucleotides N: %s\n", yesno(contains_N));
  
  background_probabilities = normalize_vector_copy(background_frequencies);

  double CG=background_probabilities[1]+background_probabilities[2];
  printf("Background distribution: %s\n", print_vector(background_probabilities, " ", 2).c_str());
  printf("CG content: %lg\n", CG);

  // if (prior_parameter == "addone")
  //   pseudo_counts.use_add_one();
  // else if (prior_parameter == "dirichlet")
  //   pseudo_counts.use_dirichlet(0.01, background_probabilities);
  
  
  printf("Use positional model for background: %s\n", 
	 use_positional_background ? "yes" : "no");
  printf("Use multimer version for multinomial-1: %s\n", use_multimer ? "yes" : "no");

  // is this useless???
  printf("Use pseudo counts: %s\n", use_pseudo_counts ? "yes" : "no");

  std::string seed;
  int seed_count=0;
  if (seed_param == "") {
    if (use_multimer)
      boost::tie(seed,seed_count) = most_common_pattern_multimer(sequences, k, seed_param, contains_N, hamming_radius);
    else
      boost::tie(seed,seed_count) = most_common_pattern_monomer(sequences, k, seed_param, hamming_radius);
    printf("Using seed %s with count %i\n", seed.c_str(), seed_count);
  }
  else
    seed = seed_param;

  printf("Palindromic index of seed %s is %i\n", seed.c_str(), ::hamming_distance(seed, reverse_complement(seed)));
  
  /*
  dmatrix seed_matrix = make_consensus_matrix(seed, seed_count);
  write_matrix(stdout, seed_matrix, "seed matrix counts:\n", "%.0f");
  normalize_matrix_columns(seed_matrix);
  write_matrix(stdout, seed_matrix, "seed matrix:\n", "%.6f");
  */

  if (false) {
    dmatrix result = align_all(sequences);
    write_matrix(stdout, result, to_string("All aligned counts:\n"), "%.0f");
    dmatrix norm=result;
    normalize_matrix_columns(norm);
    write_matrix(stdout, norm, to_string("All aligned matrix:\n"), "%.6f");
  }

  ///////////////////////////////////////////////////////////////////////////////////////////
  //
  // Start the computations:
  //
  ///////////////////////////////////////////////////////////////////////////////////////////

  if (not use_multimer and L < 60)
    distribution_of_startpositions(sequences, seed);

  int lines_bg=0, bad_lines_bg;
  std::vector<std::string> background_seqs;
  double lambda=0.0;
  if (use_background_correction) {
    printf("\n");
    boost::tie(lines_bg, bad_lines_bg) = read_sequences(background_file, background_seqs);
    check_data(background_seqs);
    printf("Read %i good lines from background file %s\n", 
	   lines_bg, background_file.c_str());
    printf("Input contained %d bad lines\n", bad_lines_bg);

    printf("Computing background lambda:\n");
    lambda = jussi_lambda(background_seqs, sequences);
    if (lambda > 1)
      lambda = 1.0;
    printf("Background lambda is %f\n", lambda);

  }

  /////////////////////////////////////
  //
  // compute the init_motif matrix
  //
  /////////////////////////////////////

  if (false) {
    dmatrix init_motif(4,k);

    if (is_nucleotide_string(seed)) {
      init_motif = find_hamming_neighbourhood(sequences, seed, hamming_radius);
      write_matrix(stdout, init_motif, "initmotif matrix counts:\n", "%.0f");
      normalize_matrix_columns(init_motif);
      write_matrix(stdout, init_motif, "initmotif matrix:\n", "%.6f");
    }
  }
 
  /////////////////////////////////////
  //
  // compute the multinomial-1 matrix
  //
  /////////////////////////////////////
  
  printf("\n");

  func_ptr_t func_ptr = use_multimer ? find_snips_multimer : find_snips_monomer;
  dmatrix multinomial1_motif;

  if (k <= 64) 
    multinomial1_motif = find_model(func_ptr, "multinomial-1", seed, 1, sequences, background_seqs, background_probabilities,
					    lambda);



  /////////////////////////////////////
  //
  // compute the multinomial-2 matrix
  //
  /////////////////////////////////////

  //  dmatrix multinomial2_motif = find_model(find_multinomial2, "multinomial2", seed, 2, sequences, background_seqs,
  //				  lambda);



  /////////////////////////////////////
  //
  // compute the multinomial-n matrix
  //
  /////////////////////////////////////
 
  if (k >= 20 or hamming_radius >= 4) 
    func_ptr = find_multinomial_n_scan;
  else
    func_ptr = find_multinomial_n;

  dmatrix multinomial_n_motif= find_model(func_ptr, "multinomial-n", seed, hamming_radius, sequences, background_seqs, background_probabilities,
				  lambda);

  /////////////////////////////////////
  //
  // compute the dinucleotide model
  //
  /////////////////////////////////////

  if (use_dinucleotide_model) {
    typedef dmatrix (*func_ptr_t)(const std::string&, const std::vector<std::string>&, int);
    int k = seed.length();
    func_ptr_t func_ptr;
    if (k >= 20 or hamming_radius >= 4) 
      func_ptr = dinucleotide_counts_scan;
    else
      func_ptr = dinucleotide_counts_suffix_array;

    dinuc_model dm = find_dinucleotide_model(func_ptr, "dinucleotide-n", seed, hamming_radius, sequences, background_seqs,
					     lambda);

    //dinuc_model dm(sequences, seed, hamming_distance);
    //    dm.print("","");
    dm.print();
  }


  /////////////////////////////////////
  //
  // compute the positional background
  //
  /////////////////////////////////////

  if (not use_multimer and L < 60) {
    positional_background = count_positional_background(sequences);
    
    normalize_matrix_columns(positional_background);
    assert(is_column_stochastic_matrix(positional_background));
    write_matrix(stdout, positional_background, "Positional background probility matrix:\n", "%.6f");
  }


  /////////////////////////////////////
  //
  // compute matrix dependencies
  //
  /////////////////////////////////////

  if (find_dependence_matrix) {
    //    prior<double> pseudo_counts;
    pseudo_counts.use_dirichlet(0.01, background_probabilities);
    
    dmatrix temp = multinomial1_motif;
    pseudo_counts.add(temp);
    normalize_matrix_columns(temp);

    dmatrix dependence_matrix = get_dependence_matrix(sequences, seed, temp);
    write_matrix(stdout, dependence_matrix, "Chi square values:\n", "%10lf");
    dependence_matrix.apply(chi_square(15));  // turn chi^2 values to p-values
    write_matrix(stdout, dependence_matrix, "P-values 1:\n", "%10.2le");
  }




  // print_submotifs(used_motif);

  // store result matrix for EM-algorithm and for dimer computations
  if (matrixfile != "-") {
    write_matrix_file(matrixfile, multinomial_n_motif, "%.0f");
    /*
    if (hamming_distance == 1)
      write_matrix_file(matrixfile, multinomial1_motif);
    else
      write_matrix_file(matrixfile, multinomial_n_motif);
    */
  }

  if (statistics) {
    print_hamming_statistics(sequences, seed, hamming_radius);
    print_hamming_statistics2(sequences, seed, hamming_radius);
  }
  
  return 0;
}


