#include "dinucleotide.hpp"
#include "iupac.hpp"
#include "common.hpp"
#include "matrix.hpp"
#include "parameters.hpp"
#include "matrix_tools.hpp"
#define TIMING 1
#include "timing.hpp"
#include "data.hpp"
#include "my_assert.hpp"
#include "aho_corasick_wrapper.hpp"
#include "suffix_array_wrapper.hpp"

#include <cstring>
#include <set>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>

//typedef boost::unordered_map<int, int> id_to_count_type;
//id_to_count_type id_to_count;

//typedef boost::unordered_map<big_int, std::vector<boost::tuple<int, int, int, int> > > id_to_tuple_type;
//id_to_tuple_type id_to_tuple;

extern bool use_multimer;

dinuc_model
reverse_complement(const dinuc_model& dm)
{
  // complement the two-nucleotide string
  int transform[] = {15, 11, 7, 3, 14, 10, 6, 2,
		     13, 9, 5, 1, 12, 8, 4, 0};


  dmatrix dm_new = dm.dm;
  int len = dm.length() - 1;

  for (int j=0; j < len; ++j) {
    for (int c=0; c < 16; ++c) {
      dm_new(transform[c],len-j-1) = dm.dm(c, j);
    }
  }
  
  dinuc_model result;
  result.init(dm_new);

  return result;
}

dinuc_model
pwm_to_dinucleotide(const dmatrix& pwm)
{
  int rows, cols;
  boost::tie(rows, cols) = pwm.dim();
  assert(rows == 4);

  dmatrix result(16, cols-1);

  for (int j=0; j < cols-1; ++j) {
    for (int c = 0; c < 16; ++c) {
      result(c, j) = pwm(c>>2,j) * pwm(c&3,j+1);
    }
  }

  dinuc_model dinuc;
  dinuc.init(result);
  
  return dinuc;
}

namespace {

  // The set of distinct strings in the 2-Hamming neighbourhood.
  // The resulting set should be of size 1 + 3*k + binom(k,2)*3^2
  /*
  std::set<std::string>
  two_hamming_neighbourhood(const std::string& s)
  {
    std::set<std::string> result;
    const char nucs[] = "ACGT";
    int k = s.length();

    for (int i=0; i < k-1; ++i) {
      for (int j=i+1; j < k; ++j) {
	std::string temp = s;
	for (int c=0; c < 16; ++c) {  // all dinucleotides
	  temp[i] = nucs[c >> 2];
	  temp[j] = nucs[c & 3];
	  result.insert(temp);	
	}
      }
    }

    assert(result.count(s));
    return result;
  }
  
  // callback function for Aho-Corasick automaton
  int match_handler2(MATCH* m, void* param)
  {
    assert(m->match_num == 1);
    big_int id = m->matched_strings[0].id;
    std::map<big_int, int>* result = (std::map<big_int, int>*)param;
    (*result)[id] += 1;

    // to find all matches always return 0
    return 0;

  }

  std::map<big_int, int>
  string_counts(const std::vector<std::string>& sequences, 
		const std::set<std::string>& strings)
  {
    typedef std::set<std::string>::iterator iterator;
    aho_corasick ac(match_handler2);
    for (iterator it=strings.begin(); it != strings.end(); ++it) {
      const std::string& temp = *it;
      ac.add_string(temp, dna_to_number(temp));
    }

    std::string str1=join(sequences, '#');
    std::string str2=join_rev(sequences, '#');
    std::map<big_int, int> result;
    ac.search(str1, &result);
    ac.search(str2, &result);

    return result;
  }
  */
  
  /*
  // The dinucleotide counts in all possible two adjacent positions
  dmatrix
  dinucleotide_counts(const std::string& seed, const std::vector<std::string>& sequences)
  {
    int k = seed.length();
    dmatrix result(16, k-1);
    const char nucs[] ="ACGT";

    std::set<std::string> strings = two_hamming_neighbourhood(seed);
    // printf("Size of 2-Hamming neighbourhood is %zu distinct strings\n", strings.size());
    std::map<big_int, int> counts = string_counts(sequences, strings);

    for (int i=0; i < k-1; ++i) {
      std::string temp = seed;
      for (int c=0; c < 16; ++c) {
	temp[i] = nucs[c >> 2];
	temp[i+1] = nucs[c & 3];
	result(c, i) = counts[dna_to_number(temp)];
      }
    }

    return result;
  }
  */
} // end nameless namespace

  // For this to be unbiased n has to be at least 2
dmatrix
dinucleotide_counts_suffix_array(const std::string& seed, const std::vector<std::string>& sequences, int n)
{
  TIME_START(t);
  const int k = seed.length();
  const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  char nucs[] = "ACGT";
  std::string str1;
  std::string str2;

  str1=join(sequences, '#');
  if (use_two_strands) {
    str1.append("#");
    str1 += join_rev(sequences, '#');
  }

  dmatrix result(16, k-1);

  //    code_to_tuple_type code_to_tuple;

  typedef boost::unordered_map<std::string, std::vector<boost::tuple<int, int> > > string_to_tuple_type;

  string_to_tuple_type string_to_tuple;

  std::vector<std::string> complements(k);  // These are set complements, not nucleotide complements
  std::vector<int> bases(k);
  unsigned long long N_mask=0;                // bitmask for positions that contain 'N'. Those positions cannot contain an error
  for (int i=0; i < k; ++i) {
    complements[i] = complement_set(seed[i]);
    bases[i]=complements[i].length() - 1;
    if (seed[i]=='N')
      N_mask |=  (1 << (k-1-i));
  }
  for (int j=0; j < k-1; ++j) {        // iterate through all dinucleotide start positions
    std::string temp = seed;         
    for (int a=0; a < 16; ++a) {      // this handles hamming distances 0, 1 and 2
      if (n == 0 and not (iupac_match(nucs[a>>2], seed[j]) and iupac_match(nucs[a&3], seed[j+1])))
	continue;
      if (n == 1 and not (iupac_match(nucs[a>>2], seed[j]) or iupac_match(nucs[a&3], seed[j+1])))
	continue;
      temp[j]=nucs[a>>2];
      temp[j+1]=nucs[a&3];
      string_to_tuple[temp].push_back(boost::make_tuple(j, a));
    }

    for (int error=1; error < n-1; ++error) { // errors outside position j, handles hamming distances 3 <= hd <= n
      // bitvector c has 1-bit for each member of the subset, rightmost bit is bit number k-1
      unsigned long long c = (1ull<<error)-1;  // initially rightmost 'error' bits are 1
      int mycount = 0;
      // iterate through all subsets c of {0, ..., k-1} that have size 'error'
      while (c < (1ull<<k)) {   // Superset has only k elements
	assert(__builtin_popcount(c) == error);
	if (((c & (1ull << (k-1-j)))) == 0 and ((c & (1ull << (k-2-j)))) == 0 and ((c & N_mask) == 0))  { // j and j+1 don't belong to the subset, and subset positions don't contain 'N'
	  ++mycount;
	  std::vector<int> P;  // positions that contain error
	  int number_of_combinations = 1;
	  std::string temp = seed;
	  for (int pos=0; pos < k; ++pos) {
	    if ((c & (1ull << (k-1-pos))) != 0) {   // pos belongs to set c
	      P.push_back(pos);
	      temp[pos] = complements[pos][0];  // initialize to first string that has mismatches at these positions
	      number_of_combinations *= complements[pos].length();
	    }
	  }
	  
	  std::vector<int> y(error, 0);
	  y[error-1]=-1;  // Initialize
	  for (int r=0; r < number_of_combinations; ++r) {
      
	    int i;
	    for (i=error-1; y[i] == bases[P[i]]; --i) {
	      y[i]=0;
	      temp[P[i]] = complements[P[i]][y[i]];
	    }
	    y[i]++;
	    temp[P[i]] = complements[P[i]][y[i]];

	    for (int a=0; a < 16; ++a){
	      temp[j] = nucs[a>>2];
	      temp[j+1] = nucs[a&3];
	      string_to_tuple[temp].push_back(boost::make_tuple(j, a));
	    }
	  }  // end for r

	}    // end if j not in c

	unsigned long long a = c&-c;
	unsigned long long b = c+a;   // update bitvector c. This is "Gosper's hack"
	c = (c^b)/4/a|b;

      } // end foreach subset c


    }  // end for error
  }

  suffix_array sa(str1);
  printf("Number of patterns %lu\n", string_to_tuple.size());
  unsigned long seed_count=0;
  unsigned long total_count=0;
  //unsigned long number_of_sites=0;

  if (use_multimer) {
    unsigned long count = sa.count_iupac(seed);
    total_count += count;
    seed_count = count;

    std::vector<boost::tuple<int, int> > pairs;
    std::string pattern;
    BOOST_FOREACH(boost::tie(pattern, pairs), string_to_tuple) {
      unsigned long count = sa.count_iupac(pattern);
      if (not iupac_string_match(pattern, seed))
	total_count += count;
      for (int i=0; i < pairs.size(); ++i) {
	int j, a;
	boost::tie(j,a) = pairs[i];
	result(a, j) += count;
      }
    }
  } else {  // allow only single occurrence per read
    typedef string_to_tuple_type::iterator iterator;
    int lines = sequences.size();
    std::vector<std::set<int> > hit_positions(lines);
    std::vector<std::vector<iterator> > hit_patterns(lines);
    int divisor = L + 1;
    std::vector<boost::tuple<int, int> > pairs;
    std::string pattern;
    for (iterator it=string_to_tuple.begin(); it != string_to_tuple.end(); ++it) {
      std::vector<long int> positions;
      pattern = it->first;
      sa.locate_iupac(pattern, positions);
      BOOST_FOREACH(int pos, positions) {
	int i = pos / divisor;  // index of the read containing the pos
	int j = pos % divisor;
	if (i >= lines) {       // handle reverse complement
	  i -= lines;
	  j = L - j - 1;
	}
	hit_positions[i].insert(j);
	//	  if (hits[i] == 1)
	hit_patterns[i].push_back(it);
      }
    }
    for (int i=0; i < lines; ++i) {
      if (hit_positions[i].size() != 1)  // only one site per sequence allowed, BUT DOES THIS HANDLE PALINDROMES CORRECTLY!!!!!
	continue;
      for (int t=0; t < hit_patterns[i].size(); ++t) { // but several different patterns can map to this site
	boost::tie(pattern, pairs) = *(hit_patterns[i][t]);
	if (not iupac_string_match(pattern, seed))
	  total_count += 1;
	else
	  seed_count += 1;
	for (int s=0; s < pairs.size(); ++s) {
	  int j, a;
	  boost::tie(j,a) = pairs[s];
	  result(a, j) += 1;
	}
      }
    }
    // for (int a=0; a < 4; ++a)             // get the seed count, THIS IS NOT CORRECT!!!!!!!!!
    // 	if (iupac_match(nucs[a], seed[0]))
    // 	  seed_count += result(a,0);
    total_count += seed_count;
  }
  printf("Seed %s count = %lu\n", seed.c_str(), seed_count);
  printf("Total dinucleotide-n count is %lu\n", total_count);
    

  TIME_PRINT("Dinucleotide-n algorithm took %.2f seconds\n", t);

  return result;
} // dinucleotide_counts_suffix_array


// this computes the dinucleotide-n matrix counts by scanning all possible windows
dmatrix
dinucleotide_counts_scan(const std::string& seed, const std::vector<std::string>& sequences, int n)
{
  TIME_START(t);
  const int k = seed.length();
  const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  //  char nucs[] = "ACGT";
  dmatrix result(16, k-1);
  for (int i=0; i < sequences.size(); ++i) {
    int max_dir = use_two_strands ? 2 : 1;
    for (int dir=0; dir < max_dir; ++dir) {
      const std::string& line = dir == 0 ? sequences[i] : reverse_complement(sequences[i]);
      for (int j=0; j < L-k+1; ++j) {
	int hd = iupac_hamming_dist(line.substr(j, k), seed, n);
	if (hd > n)
	  continue;
	for (int pos=0; pos < k-1; ++pos) {
	  char c1 = line[j+pos];
	  char c2 = line[j+pos+1];
	  int cc = iupac_match(c1, seed[pos]) ? 0 : 1;
	  cc += iupac_match(c2, seed[pos+1]) ? 0 : 1;
	  if (hd-cc <= n-2)
	    ++result((to_int(c1)<<2) + to_int(c2), pos);
	}
      }
    }
  }

  TIME_PRINT("Dinucleotide-n scanning algorithm took %.2f seconds\n", t);
  return result;
}

namespace {

  /*
  // NOT USED ANYWHERE!!!!!!!!!!!
  // tells the 2-bit alphabet Hamming distance between x and y
  int
  number_of_unequal_nucleotides(unsigned x, unsigned y)
  {
    const int l = sizeof(unsigned)*4;
    unsigned even = 0;
    for (int i=0; i < l; ++i) {
      even <<= 2;
      even += 1;
    }
    unsigned odd = even << 1;
  
    unsigned tmp = x^y;
    unsigned result = (tmp & even) | ((tmp & odd) >> 1);
    return __builtin_popcount(result);
  }
  */

  // transforms dinucleotide counts in the dmatrix(16,k-1) to 
  // initial_probabilities dmatrix(4,k)
  dmatrix
  initial_probabilities(const dmatrix& dm)
  {
    assert(dm.get_rows() == 16);
    int k = dm.get_columns() + 1;
    dmatrix result(4, k);

    for (int i=0; i < k-1; ++i) {     
      for (int c=0; c < 16; ++c) {
	result(c >> 2, i) += dm(c, i);   // first character of dinucleotide 
      }
    }
    // handle the last position
    for (int c=0; c < 16; ++c) {
      result(c & 3, k-1) += dm(c, k-2);   // second character of dinucleotide 
    }

    normalize_matrix_columns(result);

    return result;
  }


  // transforms dinucleotide counts in the dmatrix(16,k-1) to 
  // conditional_probabilities dmatrix(16,k-1)
  dmatrix
  conditional_probabilities(const dmatrix& dm)
  {
    assert(dm.get_rows() == 16);
    int k = dm.get_columns()+1;
    dmatrix result(16, k-1);

    // do normalization
    for (int i=0; i < k-1; ++i) {     
      for (int first=0; first < 4; ++first) {
	double sum = 0;
	for (int second=0; second < 4; ++second)
	  sum += dm((first<<2) + second, i);
	if (sum > 0)
	  for (int second=0; second < 4; ++second)
	    result((first<<2) + second, i) = dm((first<<2) + second, i) / sum;
      }
    }

    return result;
  }

} // end empty namespace



void
dinuc_model::init(const dmatrix& dm_)
{
  dm=dm_;
  ip = initial_probabilities(dm);
  cp = conditional_probabilities(dm);

}

dinuc_model::dinuc_model(const std::string& filename) 
{
  init(read_matrix_file(filename));
}

dinuc_model::dinuc_model(const std::vector<std::string>& sequences, const std::string& seed, int n) 
{
  typedef dmatrix (*func_ptr_t)(const std::string&, const std::vector<std::string>&, int);
  int k = seed.length();
  func_ptr_t func_ptr;
  if (k >= 20 or n >= 4) 
    func_ptr = dinucleotide_counts_scan;
  else
    func_ptr = dinucleotide_counts_suffix_array;


  init(func_ptr(seed, sequences, n));
}

double
dinuc_model::cond(int i, int c1, int c2) const
{
  int dinuc = (c1 << 2) + c2;
    
  return cp(dinuc, i);
}

int
dinuc_model::length() const
{ return ip.get_columns(); }

void
dinuc_model::print() const
{
  const char* headers[] = {"AA","AC","AG","AT",
			   "CA","CC","CG","CT",
			   "GA","GC","GG","GT",
			   "TA","TC","TG","TT"};
  int k = length();
  for (int i=0; i < k-1; ++i) {
    printf("\t%i", i);
  }
  printf("\n");
   
  for (int c=0; c < 16; ++c) {
    printf("%s", headers[c]);
    for (int i=0; i < k-1; ++i) {
      printf("\t%i", (int)dm(c, i));
    }
    printf("\n");
  }

  printf("Initial probabilities:\n");
  ip.print();
  printf("Conditional probabilities:\n");
  cp.print();
}

double
dinuc_model::score(const std::string& kmer, int start_pos) const
{
  assert(start_pos < length());
  std::string temp;
  double score = 0.0;
  if (start_pos < 0) {
    start_pos = -start_pos;
    temp = kmer.substr(start_pos);
  }
  else
    temp = kmer;
  int len = temp.length();
  score += log2(ip(to_int(temp[0]), start_pos) / 0.25);
  for (int i=0; i < std::min(len, length() - start_pos)-1; ++i)
    score += log2(cond(start_pos+i, to_int(temp[i]), to_int(temp[i+1])) / 0.25);

  return score;
} 


