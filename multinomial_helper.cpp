#define TIMING
#include "timing.hpp"

#include "multinomial_helper.hpp"
#include "bndm.hpp"
#include "common.hpp"
#include "parameters.hpp"
#include "iupac.hpp"
#include "matrix_tools.hpp"
#include "kmer_tools.hpp"


#include <boost/foreach.hpp>

//extern bool use_palindromic_correction;

int palindromic_index_limit=0;
//bool use_palindromic_correction = false;
int low_count_limit = 20;

int
conflict_free_palindromic_index(int hamming_radius)
{
  int limit;
  switch (hamming_radius) {
  case 0: limit = 0; break;
  case 1: limit = 2; break;
  default: limit = 2*hamming_radius + 1;
  }
  return limit;
}

std::vector<std::string>
remove_masked_areas(const std::vector<std::string>& sequences, int k)
{
  std::vector<std::string> result;
  int lines = sequences.size();
  for (int i=0; i < lines; ++i) {
    const std::string& line = sequences[i];
    std::string current;
    for (int j=0; j < line.length(); ++j) {
      if (line[j] == 'N') {
	if (current.length() >= k)
	  result.push_back(current);
	if (current.length() > 0)
	  current.clear();
      }
      else
	current.push_back(line[j]);
    }
    if (current.length() >= k)
      result.push_back(current);
      
  }
  return result;
}

// The triples are of (code, count, palindromic index) type 
// Primary sort key is the count, the secondary sort key is the palindromic index
bool
triple_comp(boost::tuple<big_int, int, int> a, boost::tuple<big_int, int, int> b)
{
  return a.get<1>() > b.get<1>() or (a.get<1>() == b.get<1>() and a.get<2>() > b.get<2>());
}

boost::tuple<std::string,int>
most_common_pattern_multimer(const std::vector<std::string>& sequences, int k, std::string seed,
			     bool contains_N, int hamming_radius)
{
  assert(k > 1 && k <= max_matrix_len); 
  //int strings = (int)pow(4,k); // if k==12 then this is about 16M

  //  std::vector<unsigned> number_of_occurrences(strings);
  count_container_t number_of_occurrences;
  bool count_palindromes_twice = use_two_strands;  // There is also global version of this
  if (contains_N)
    get_kmer_counts(remove_masked_areas(sequences, k), k, number_of_occurrences, use_two_strands, count_palindromes_twice);
  else
    get_kmer_counts(sequences, k, number_of_occurrences, use_two_strands, count_palindromes_twice);

  
  printf("Size of number_of_occurrences container %zu\n", number_of_occurrences.size());

  //count_t count;
  //code_t code=0;
  /*
  BOOST_FOREACH(boost::tie(code, count), number_of_occurrences)
    printf("*%llu\t%s\t%i\n", (unsigned long long) code, number_to_dna(code,k).c_str(), count);
  
  printf("Size of number_of_occurrences container %zu\n", number_of_occurrences.size());
  */
  
  code_t argmax;
  std::vector<boost::tuple<big_int, int, int> > v;
  code_t code;
  int count;
  BOOST_FOREACH(boost::tie(code, count), number_of_occurrences)
    v.push_back(boost::make_tuple(code, count, palindromic_index(number_to_dna(code, k))));
  std::sort(v.begin(), v.end(), triple_comp);   // compares according the second member of the pair: the count

  if (palindromic_index_limit > 0) {
    // the palindromic index of a seed needs to be at least 'limit' for the n-Hamming-neighbourhood to be conflict-free
    code_t code=0;
    code_t max_code = 0;
    int max_pi = -1;
    for (int i=0; i < v.size(); ++i) {
      code = v[i].get<0>();
      int pi = v[i].get<2>();
      int count = v[i].get<1>();
      if (count < low_count_limit and max_pi >= 0)   // If we have already one candidate and the current kmer is too small, then quit.
	break;
      if (pi > max_pi) {
	max_pi = pi;
	max_code = code;
      }
      if (max_pi >= palindromic_index_limit)
	break;                   // found a good seed
    }
    argmax = max_code;
  }
  else {
    argmax = v[0].get<0>();
  }
  
  //  int max_count = number_of_occurrences[argmax];
  std::string result;

  // between string and its reverse complement, choose lexicographically smaller
  if (seed.length() == k) {
    result = seed;
  }
  else {
    std::string result1 = number_to_dna(argmax,k);
    std::string result2 = reverse_complement(result1);
    if (use_two_strands)
      result = (result1 < result2) ? result1 : result2;
    else
      result = result1;
  }

  printf("Seed %s has %i occurences\n",
	 result.c_str(),number_of_occurrences[dna_to_number(result)]);

  return boost::make_tuple(result, number_of_occurrences[dna_to_number(result)]);
  //return result;
}


boost::tuple<std::string,int>
most_common_pattern_monomer(const std::vector<std::string>& sequences, int k, std::string seed,
			    int hamming_radius)
{
  assert(k > 1 && k <= max_matrix_len); 
  //  int strings = (int)pow(4,k); // if k==12 then this is about 16M
  //std::vector<int> number_of_occurences(strings);
  //typedef std::map<big_int, int> my_container;
  typedef boost::unordered_map<big_int, int> my_container; // items are (code,count) pairs
  my_container number_of_occurrences;
  int lines = sequences.size();
  int max_count=-1;
  big_int argmax=-1;
  big_int id;
  big_int id2;
  bool count_palindromes_twice = use_two_strands;  // There is also global version of this
  for  (int i=0; i < lines; ++i) {
    const std::string& line = sequences[i];
    //std::map<big_int, std::vector<int> >occurences;

    boost::unordered_map<big_int, std::vector<int> > occurrences;  // occurrences on this line
    std::set<big_int> ids;

    // find all occurrences in sequence
    for (int j=0; j < line.length()-k+1; ++j) {
      id = dna_to_number(line.substr(j,k));
      occurrences[id].push_back(j);
      ids.insert(id);
      if (use_two_strands) {
	id2 = dna_to_number(reverse_complement(line.substr(j,k)));
	occurrences[id2].push_back(j);
      }
    }
    // accept only subsequences that appear only once per sequence, or is a palindromic occurrence
    for (std::set<big_int>::iterator it=ids.begin(); it != ids.end(); ++it) {
      big_int id = *it;
      std::vector<int>& r = occurrences[id];
      if (r.size() == 1 || (count_palindromes_twice && r.size() == 2 && r[0] == r[1])) {
	big_int id2 = reverse_complement_2bitstring(id, k);
	//	std::string s = number_to_dna(id, k);
	//	big_int id2 = dna_to_number(reverse_complement(s));
	++number_of_occurrences[id];
	if (use_two_strands)
	  ++number_of_occurrences[id2];
	if (number_of_occurrences[id] > max_count) {
	  max_count = number_of_occurrences[id];
	  argmax = id;
	}
	if (use_two_strands && number_of_occurrences[id2] > max_count) {
	  max_count = number_of_occurrences[id2];
	  argmax = id2;
	}	  
      }

    }
  }

  // print top10 of strings
  // triples are (code, count, palindromic index)
  std::vector<boost::tuple<big_int, int, int> > v;
  my_container::iterator it;
  for (it=number_of_occurrences.begin(); it != number_of_occurrences.end(); ++it) {
    v.push_back(boost::make_tuple(it->first, it->second, palindromic_index(number_to_dna(it->first, k))));
  }
  std::sort(v.begin(), v.end(), triple_comp);   // compares according the second member of the pair: the count
  if (palindromic_index_limit > 0) {
    code_t code=0;
    code_t max_code = 0;
    int max_pi = -1;
    for (int i=0; i < v.size(); ++i) {
      code = v[i].get<0>();
      int pi = v[i].get<2>();
      int count = v[i].get<1>();
      if (count < low_count_limit and max_pi >= 0)   // If we have already one candidate and the current kmer is too small, then quit.
	break;
      if (pi > max_pi) {
	max_pi = pi;
	max_code = code;
      }
      if (max_pi >= palindromic_index_limit)
	break;                   // found a good seed
    }
    argmax = max_code;
  } else
    argmax = v[0].get<0>();

  // between string and its reverse complement, choose lexicographically smaller
  std::string result;
  std::string result1;
  if (seed.length() == k)
    result = seed;
  else {
    result1 = number_to_dna(argmax,k);
    std::string result2 = reverse_complement(result1);
    if (use_two_strands)
      result = (result1 < result2) ? result1 : result2;
    else
      result = result1;
  }
  printf("Seed %s has %i occurences\n",
	 result.c_str(),number_of_occurrences[dna_to_number(result)]);

  return boost::make_tuple(result, number_of_occurrences[dna_to_number(result)]);
}



// Do not reject sequences with multiple occurrences of query strings.
// Compute the counts for the multinomial1 matrix
boost::tuple<dmatrix,int,int>
find_snips_multimer_helper(const std::string& consensus, const std::vector<std::string>& sequences)
{
  
  std::string str1;
  std::string str2;

  //typedef boost::tuple<int,int,int> triple;  // (seqno, position, direction)
  //std::vector<triple> alignment;
  std::vector<std::string> alignment;

  //  int lines = sequences.size();

  str1=join(sequences, '#');
  str2=join_rev(sequences, '#');

  int k = consensus.length();
  char nucs[] = "ACGT";
  bool is_palindrome = is_palindromic(consensus);
  int consensus_count = BNDM_with_joker(str1, consensus);
  if (use_two_strands && (count_palindromes_twice || not is_palindrome))
    consensus_count += BNDM_with_joker(str2, consensus);

  if (print_alignment) {
    for (int t=0; t < consensus_count; ++t)
      alignment.push_back(consensus);
  }

  matrix<double> result(4, k);
  for (int j=0; j < k; ++j) {        // iterate through all string positions
    std::string temp = consensus;
    for (int a=0; a < 4; ++a) {      // iterate through all characters
      if (nucs[a] == consensus[j])   // the consensus count was already computed. NOTE: no iupac_match here, on purpose
	continue;
      temp[j]=nucs[a];
      is_palindrome = is_palindromic(temp);
      
      result(a,j) = BNDM_with_joker(str1,temp);
      if (use_two_strands && (count_palindromes_twice || not is_palindrome))
	result(a,j) += BNDM_with_joker(str2,temp);

      if (print_alignment) {
	for (int t=0; t < result(a,j); ++t)
	  alignment.push_back(temp);
      }
    }
  }


  int total_count = consensus_count;
  for (int j=0; j < k; ++j) {       // iterate through all string positions
    for (int a=0; a < 4; ++a) {     // iterate through all characters
      if (nucs[a] == consensus[j]) {
	result(to_int(consensus[j]), j) = consensus_count;   // add consensus count to the matrix
	if (print_alignment) {
	  for (int t=0; t < consensus_count; ++t)              // add consensus to alignment also
	    alignment.push_back(consensus);                    // for other columns
	}
      }
      else if (not iupac_match(nucs[a], consensus[j])) // the consensus count was already added to the total_count
	total_count += result(a, j);   // number of sequences used for the matrix
    }
  }


  // print the alignment to file descriptor 3, if it is open
  if (print_alignment) {
    FILE* fp = fdopen(3, "a");
    if (fp != NULL) {
      for (int t = 0; t < alignment.size(); ++t)
	fprintf(fp, "%s\n", alignment[t].c_str());
      fclose(fp);
    }
  }

  return boost::make_tuple(result, consensus_count, total_count);
} // find_snips_multimer


dmatrix
find_snips_multimer(const std::string& consensus, const std::vector<std::string>& sequences, int hamming_distance)
{
  TIME_START(t);
  assert(hamming_distance == 1);
  dmatrix result;
  int seed_count;
  int multinomial_count;
  boost::tie(result, seed_count, multinomial_count) = find_snips_multimer_helper(consensus, sequences);
  TIME_PRINT("Multinomial-1 algorithm took %.2f seconds.\n", t);
  printf("Seed count = %i\n", seed_count);
  printf("Total multinomial1 count is %d\n", multinomial_count);
  return result;
}

string_to_tuple_type
get_n_neighbourhood(const std::string&seed, int n)
{
  const int k = seed.length();
  //const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  char nucs[] = "ACGT";

  //code_to_tuple_type code_to_tuple;
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
  for (int j=0; j < k; ++j) {        // iterate through all string positions
    std::string temp = seed;         // Variable temp is not really needed. It's just for sanity check.
    //packed_string myid2(seed);
    for (int a=0; a < 4; ++a) {      // this handles hamming distances 0 and 1
      if (n == 0 and not iupac_match(nucs[a], seed[j]))
	continue;
      temp[j]=nucs[a];
      // myid2[j] = a;
      // my_assert(myid2.get_bits(), dna_to_number(temp));
      // code_to_tuple[myid2.get_bits()].push_back(boost::make_tuple(j, a));
      string_to_tuple[temp].push_back(boost::make_tuple(j, a));
    }

    for (int error=1; error < n; ++error) { // errors outside position j, handles hamming distances 2 <= hd <= n
      // bitvector c has 1-bit for each member of the subset, rightmost bit is bit number k-1
      unsigned long long c = (1ull<<error)-1;  // initially rightmost 'error' bits are 1
      int mycount = 0;
      // iterate through all subsets c of {0, ..., k-1} that have size 'error'
      while (c < (1ull<<k)) {   // Superset has only k elements
	assert(__builtin_popcountll(c) == error);
	if (((c & (1ull << (k-1-j)))) == 0 and ((c & N_mask) == 0))  { // j doesn't belong to the subset, and subset positions don't contain 'N'
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
	  //packed_string myid4(seed);
	  
	  std::vector<int> y(error, 0);
	  y[error-1]=-1;  // Initialize
	  for (int r=0; r < number_of_combinations; ++r) {
      
	    int i;
	    for (i=error-1; y[i] == bases[P[i]]; --i) {
	      y[i]=0;
	      // temp[P[i]] = nucs[myskip(y[i], skip[i])];
	      // myid4[P[i]] = myskip(y[i], skip[i]);
	      temp[P[i]] = complements[P[i]][y[i]];
	    }
	    y[i]++;
	    temp[P[i]] = complements[P[i]][y[i]];
	    // temp[P[i]] = nucs[myskip(y[i], skip[i])];
	    // myid4[P[i]] = myskip(y[i], skip[i]);

	    for (int a=0; a < 4; ++a){
	      temp[j] = nucs[a];
	      //myid4[j] = a;
	      //my_assert(myid4.get_bits(), dna_to_number(temp));
	      //code_to_tuple[myid4.get_bits()].push_back(boost::make_tuple(j, a));
	      string_to_tuple[temp].push_back(boost::make_tuple(j, a));
	    }
	  }  // end for r

	}    // end if j not in c

	unsigned long long a = c&-c;
	unsigned long long b = c+a;   // update bitvector c. This is "Gosper's hack"
	c = (c^b)/4/a|b;

      } // end foreach subset c


    }  // end for error
  } // end for j

  return string_to_tuple;
}

// Returns a vector with the seed pattern at the first index
std::vector<std::pair<std::string, std::vector<boost::tuple<int, int> > > >
get_n_neighbourhood_in_vector(const std::string&seed, int n)
{
  std::vector<std::pair<std::string, std::vector<boost::tuple<int, int> > > > result;
  string_to_tuple_type neigh = get_n_neighbourhood(seed, n);
  std::pair<std::string, std::vector<boost::tuple<int, int> > > t;
  string_to_tuple_type::iterator it = neigh.find(seed);   // Make sure that the pair corresponding to the seed is first on the vector
  result.push_back(*it);
  neigh.erase(it);
  BOOST_FOREACH(t, neigh) {
     result.push_back(t);
  }
  return result;
}

class iupac_probability_in_background
{
public:
  iupac_probability_in_background(const std::vector<double>& bg)
    : iupac_probabilities(256)
  {
    BOOST_FOREACH(char iupac_char, iupac_chars) {
      BOOST_FOREACH(char c, iupac_class(iupac_char)) {
	iupac_probabilities[(unsigned char)iupac_char] += bg[to_int(c)];
      }
    }
  }

  double
  operator()(const std::string& s) {
    assert(is_iupac_string(s));
    double result = 1.0;
    for (int i=0; i < s.length(); ++i)
      result *= iupac_probabilities[s[i]];
    return result;
  }
  
private:
  std::vector<double> iupac_probabilities;
};

double
palindromic_correction(const std::string& pattern, const std::string& seed, const std::string& seed_rev)
{
  double t=1.5;
  int hd1=hamming_distance(pattern,seed);
  int hd2=hamming_distance(pattern, seed_rev);
  double correction = pow(t, -hd1) /
    (pow(t, -hd1) + pow(t, -hd2));

  return correction;
}

boost::tuple<dmatrix,int>
find_multinomial_n_background(const std::string& seed, const std::vector<std::string>& sequences, const std::vector<double>& bg,
			      int n, bool use_multimer)
{
  const int k = seed.length();
  const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  dmatrix result(4, k);

  string_to_tuple_type string_to_tuple;
  string_to_tuple = get_n_neighbourhood(seed, n);

  //printf("Number of patterns %lu\n", string_to_tuple.size());
  unsigned long seed_count=0;
  unsigned long total_count=0;
  int lines = sequences.size();
  int sites = lines * (L-k+1)*2;
  iupac_probability_in_background iupac_prob(bg);
  
  if (use_multimer) {
    seed_count = sites * iupac_prob(seed);
    total_count += seed_count;

    std::vector<boost::tuple<int, int> > pairs;
    std::string pattern;
    std::string seed_rev = reverse_complement(seed);
    BOOST_FOREACH(boost::tie(pattern, pairs), string_to_tuple) {
      unsigned long count = sites * iupac_prob(pattern);
      //      if (use_palindromic_correction)
      //	count *= palindromic_correction(pattern, seed, seed_rev);
      
      if (not iupac_string_match(pattern, seed))
	total_count += count;
      int j, a;
      BOOST_FOREACH(boost::tie(j,a), pairs) {
	result(a, j) += count;
      }
    }
  } else {  // allow only single occurrence per read
    error(true, "Not implemented");
  }
  return boost::make_tuple(result, total_count);
} // find_multinomial_n_background



boost::tuple<dmatrix,int>
find_multinomial_n_suffix_array(const std::string& seed, const std::vector<std::string>& sequences, const suffix_array& sa, int n, bool use_multimer)
{
  const int k = seed.length();
  const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  //char nucs[] = "ACGT";
  dmatrix result(4, k);
  //code_to_tuple_type code_to_tuple;

  string_to_tuple_type string_to_tuple;
  string_to_tuple = get_n_neighbourhood(seed, n);

  printf("Number of patterns %lu\n", string_to_tuple.size());
  unsigned long seed_count=0;
  unsigned long total_count=0;
  //unsigned long number_of_sites=0;

  std::string seed_rev = reverse_complement(seed);
  if (use_multimer) {
    seed_count = sa.count_iupac(seed);
    total_count += seed_count;

    std::vector<boost::tuple<int, int> > pairs;
    std::string pattern;
    printf("#String\tColumn\tHamming distance\tPalindrome\tCount\tMatches at col\n");
    BOOST_FOREACH(boost::tie(pattern, pairs), string_to_tuple) {
      unsigned long count = sa.count_iupac(pattern);
      bool is_palindrome = is_palindromic(pattern);
      int hd = hamming_distance(seed, pattern);
      if (is_palindrome and use_two_strands and not count_palindromes_twice)
	count /= 2;
      //      if (use_palindromic_correction)
      //	count *= palindromic_correction(pattern, seed, seed_rev);

      if (not iupac_string_match(pattern, seed))
	total_count += count;
      int j, a;
      BOOST_FOREACH(boost::tie(j,a), pairs) {
	printf("#%s\t%i\t%i\t%s\t%zu\t%s\n", pattern.c_str(), j, hd,
	       yesno(is_palindrome), count, yesno(seed[j]==pattern[j]));
	result(a, j) += count;
      }
    }
  } else {  // allow only single occurrence per read
    typedef string_to_tuple_type::iterator iterator;
    int lines = sequences.size();
    std::vector<std::set<int> > hit_positions(lines);
    std::vector<std::vector<iterator> > hit_patterns(lines);
    int divisor = L + 1;    // Includes the separator '#'
    std::vector<boost::tuple<int, int> > pairs;
    std::string pattern;
    for (iterator it=string_to_tuple.begin(); it != string_to_tuple.end(); ++it) {
      std::vector<long int> positions;
      pattern = it->first;
      bool is_palindrome = is_palindromic(pattern);
      sa.locate_iupac(pattern, positions);
      BOOST_FOREACH(int pos, positions) {
	int i = pos / divisor;  // index of the read containing the pos
	int j = pos % divisor;
	if (i >= lines) {       // handle reverse complement
	  i = 2*lines - i - 1;
	  j = L - (j + k - 1) - 1;
	}
	if (hit_positions[i].count(j) < 1 or not is_palindrome or count_palindromes_twice) {
	  hit_positions[i].insert(j);
	  //	  if (hits[i] == 1)
	  hit_patterns[i].push_back(it);
	}
      }
    }
    for (int i=0; i < lines; ++i) {
      if (hit_positions[i].size() != 1)                    // Because hit_positions[i] is a set, this doesn't exclude, for instance, palindromes
	continue;
      for (int t=0; t < hit_patterns[i].size(); ++t) {
	boost::tie(pattern, pairs) = *(hit_patterns[i][t]);
	if (not iupac_string_match(pattern, seed))
	  total_count += 1;
	else
	  ++seed_count;
	int j, a;
	BOOST_FOREACH(boost::tie(j,a), pairs) {
	  result(a, j) += 1;
	}
      }
    }
    // for (int a=0; a < 4; ++a)             // get the seed count                      // THIS PROBABLY ISN'T CORRECT, CONTAINS ALSO EXTRA COUNTS, WORKS ONLY WHEN n=1
    //   if (iupac_match(nucs[a], seed[0]))
    // 	seed_count += result(a,0);
    total_count += seed_count;
  }

  printf("Seed %s count = %lu\n", seed.c_str(), seed_count);
  printf("Total multinomial-n count is %lu\n", total_count);

  return boost::make_tuple(result, total_count);
}


dmatrix
align_all(const std::vector<std::string>& sequences)
{
  int L=sequences[0].length();
  int lines = sequences.size();
  int k = L;
  dmatrix result(4, k);
  for (int i=0; i < lines; ++i) {
    const std::string& line = sequences[i];
    for (int j=0; j < k; ++j) 
      result(to_int(line[j]), j) += 1;
    // if (use_two_strands) {
    //   const std::string& line_rev = reverse_complement(sequences[i]);
    //   for (int j=0; j < k; ++j) 
    // 	result(to_int(line_rev[j]), j) += 1;
    // }
  }

  return result;
}
