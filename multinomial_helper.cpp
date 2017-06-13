#define TIMING
#include "timing.hpp"

#include "multinomial_helper.hpp"
#include "bndm.hpp"
#include "common.hpp"
#include "parameters.hpp"
#include "probabilities.hpp"
#include "iupac.hpp"
#include "huddinge.hpp"
#include "matrix_tools.hpp"
#include "kmer_tools.hpp"

#include <boost/foreach.hpp>

bool use_multimer=false;
background_correction_type background_correction = no_correction;
counting_type data_counting = all_occurrences;
counting_type background_counting = data_counting;


//extern bool use_palindromic_correction;
//bool use_palindromic_correction = false;

int cluster_threshold=4; // The shift between neighbouring occurrences in a cluster can be at most 'cluster_threshold'
//bool use_one_per_cluster=true;  // not used
bool use_weighted_hamming = false;

// Weighted Hamming distance

class esko_distance_type {
public:

  //  esko_distance() {}

  esko_distance_type(int k_) : k(k_)
  {
    weights.resize(k);
    double b = (k-1)/2.0;  // centre point
    for (int i=0; i < k; ++i) {
      if (use_weighted_hamming)
	weights[i] = (100 - floor(fabs(i-b)));
      else
	weights[i] = 1;  // uniform weights
    }

  }

  int
  distance(const std::string& s, const std::string& iupac) const
  {

    assert(s.length() == iupac.length());
    assert(k == s.length());
    double distance = 0;
  
    for (int i=0; i < k; ++i) {
      //      if (s[i] != t[i])
      if (not iupac_match(s[i], iupac[i]))
	distance += weights[i];
    }

    return distance;
  }

  template <typename T>
  int
  distance_with_bits(T x, T y) const
  {
    static_assert(std::is_unsigned<T>::value,
		  "hamming_distance_with_bits requires unsigned parameter type");
    // Note! Here we assume that the leftmost bits are cleared and contain no thrash.
    // We could also get length of string as parameter and mask out the extra bits, just to be sure, but this would slow things down.
    static const T evenmask = detail::get_even_mask<T>();
    static const T oddmask = evenmask << 1;
    T w = x ^ y;
    T result = ((w & oddmask) >> 1) | (w & evenmask);
    int distance = 0;
    T mask = 1;
    for (int i=0; i < k; ++i, mask<<=2) {
      if (result & mask)
	distance += weights[i];
    }
    return distance;
  }

  std::vector<int>
  get_weights() const
  {
    return weights;
  }
  
private:
  int k;
  std::vector<int> weights;
};

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



// this computes the multinomial-n matrix counts by scanning all possible windows
boost::tuple<dmatrix,int>
find_multinomial_n_scan(const std::string& seed, const std::vector<std::string>& sequences, int n)
{
  //  TIME_START(t);
  const int k = seed.length();
  const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  //  char nucs[] = "ACGT";
  dmatrix result(4, k);
  int total_count = 0;
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
	++total_count;
	for (int pos=0; pos < k; ++pos) {
	  char c = line[j+pos];
	  int cc = iupac_match(c, seed[pos]) ? 0 : 1;   // 1 if mismatch
	  if (hd-cc <= n-1)
	    ++result(to_int(c), pos);
	}
      }
    }
  }

  //  TIME_PRINT("Multinomial-n scanning algorithm took %.2f seconds\n", t);
  return boost::make_tuple(result, total_count);
}



// this computes the multinomial-n matrix counts using suffix array (was AC-automaton)
boost::tuple<dmatrix,int>
find_multinomial_n(const std::string& seed, const std::vector<std::string>& sequences, int n)
{
  //  TIME_START(t);
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
  int total_count;
  if (use_suffix_array) {
     suffix_array sa(str1);
     boost::tie(result, total_count) = find_multinomial_n_suffix_array(seed, sequences, sa, n, use_multimer);
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

  //  printf("Multinomial-n algorithm took\n");
  //  TIME_CHECK(t);
  //  printf("seconds\n");
  return boost::make_tuple(result, total_count);
} // find_multinomial_n





string_to_tuple_type
get_n_neighbourhood_mononucleotide_contributions(const std::string&seed, int n)
{
  const int k = seed.length();
  //const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  char nucs[] = "ACGT";

  //code_to_tuple_type code_to_tuple;
  string_to_tuple_type string_to_contributions;

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
    for (int a=0; a < 4; ++a) {      // number of errors outside position j is zero, handles Hamming distances 0 and 1
      if (n == 0 and not iupac_match(nucs[a], seed[j]))
	continue;
      temp[j]=nucs[a];
      // myid2[j] = a;
      // my_assert(myid2.get_bits(), dna_to_number(temp));
      // code_to_tuple[myid2.get_bits()].push_back(boost::make_tuple(j, a));
      string_to_contributions[temp].push_back(boost::make_tuple(j, a));
    }

    for (int error=1; error < n; ++error) { // errors outside position j, handles hamming distances 1 <= hd <= n
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
	      string_to_contributions[temp].push_back(boost::make_tuple(j, a));
	    }
	  }  // end for r

	}    // end if j not in c

	unsigned long long a = c&-c;
	unsigned long long b = c+a;   // update bitvector c. This is "Gosper's hack"
	c = (c^b)/4/a|b;

      } // end foreach subset c


    }  // end for error
  } // end for j

  
  // Just a sanity check for the case when the seed does not contain iupac characters
  if (is_nucleotide_string(seed)) {
    std::string pattern;
    std::vector<boost::tuple<int, int> > pairs;
    
    BOOST_FOREACH(boost::tie(pattern, pairs), string_to_contributions) {
      assert(pairs.size() == k or pairs.size() == n);
    }
  }
  return string_to_contributions;
}

// Returns a vector with the seed pattern at the first index
std::vector<std::pair<std::string, std::vector<boost::tuple<int, int> > > >
get_n_neighbourhood_in_vector(const std::string&seed, int n)
{
  std::vector<std::pair<std::string, std::vector<boost::tuple<int, int> > > > result;
  string_to_tuple_type neigh = get_n_neighbourhood_mononucleotide_contributions(seed, n);
  std::pair<std::string, std::vector<boost::tuple<int, int> > > t;
  string_to_tuple_type::iterator it = neigh.find(seed);   // Make sure that the pair corresponding to the seed is first on the vector
  result.push_back(*it);
  neigh.erase(it);
  BOOST_FOREACH(t, neigh) {
     result.push_back(t);
  }
  return result;
}



/*
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
*/

class min_hamming_distance_class
{
public:

  min_hamming_distance_class () {}
  
  min_hamming_distance_class(const std::string& seed)
  {
    init(seed);
  }

  void
  init(const std::string& seed)
  {
    k=seed.length();
    esko_distance_type esko_distance(k);
    if (is_nucleotide_string(seed)) {
      code_t s = dna_to_number<code_t>(seed);
      code_t s_rev_comp = reverse_complement_2bitstring(s, k);
      number_of_sequences = pow(4, k);
      v.resize(number_of_sequences);
      if (use_two_strands) {
	for (code_t code=0; code < number_of_sequences; ++code) {
	  v[code] = std::min(esko_distance.distance_with_bits(code, s), esko_distance.distance_with_bits(code, s_rev_comp));
	}
      }
      else {
	for (code_t code=0; code < number_of_sequences; ++code) {
	  v[code] = esko_distance.distance_with_bits(code, s);
	}
      }
    }
    else {
      std::string seed_rev = reverse_complement(seed);
      number_of_sequences = pow(4, k);
      v.resize(number_of_sequences);
      if (use_two_strands) {
	for (code_t code=0; code < number_of_sequences; ++code) {
	  v[code] = std::min(esko_distance.distance(number_to_dna(code, k), seed),
			     esko_distance.distance(number_to_dna(code, k), seed_rev));
	}
      }
      else {
	for (code_t code=0; code < number_of_sequences; ++code) {
	  v[code] = esko_distance.distance(number_to_dna(code, k), seed);
	}
      }
    }
  }
  
  unsigned short
  get(code_t code) const {
    assert(code < number_of_sequences);
    return v[code];
  }
  
private:
  int k;
  unsigned int number_of_sequences;
  std::vector<unsigned short> v;
};

// f[i][h] is the probability that sequence of length i has exactly h mismatches to the seed
class prefix_hamming_distance
{
public:
  prefix_hamming_distance(const std::string& seed, const std::vector<double>& bg = std::vector<double>(4, 0.25))
    : k(seed.length())
  {
    f.resize(boost::extents[k+1][k+1]);
    f[0][0] = 1.0;
    for (int i=1; i <= k; ++i) {
      f[i][0] = f[i-1][0] * bg[to_int(seed[i-1])];
      for (int h=1; h <= k; ++h)
	f[i][h] = f[i-1][h] * bg[to_int(seed[i-1])] + f[i-1][h-1] * (1.0 - bg[to_int(seed[i-1])]);
    }
  }

  double
  get(int i, int h) const
  {
    return f[i][h];
  }

  void
  print() const
  {
    print_array_with_default_headers(stdout, f);
  }
  
private:
  int k;
  boost::multi_array<double, 2> f;
};

// e[i][h] is the probability that sequence of length i has at most h mismatches to the seed
class prefix_hamming_distance_cumulative
{
public:
  prefix_hamming_distance_cumulative(const std::string& seed, const std::vector<double>& bg = std::vector<double>(4, 0.25))
    : k(seed.length())
  {
    prefix_hamming_distance f(seed, bg);
    printf("Prefix hamming distance:\n");
    f.print();
    e.resize(boost::extents[k+1][k+1]);
    for (int i=0; i <= k; ++i) {
      e[i][0] = f.get(i,0);
      for (int h=1; h <= k; ++h)
	e[i][h] = e[i][h-1] + f.get(i, h);
    }
  }

  double
  get(int i, int h) const
  {
    return e[i][h];
  }

  void
  print() const
  {
    print_array_with_default_headers(stdout, e);
  }
  
  
private:
  int k;
  boost::multi_array<double, 2> e;
};

class cluster_probability_array {
public:

  typedef boost::multi_array<double, 2>::extent_range range;

  cluster_probability_array() {}
  
  // d is the hamming radius, e is the smallest Hamming distance to the seed in the cluster.
  // e must be atteined only in the end of the cluster.
  cluster_probability_array(const std::string& seed, int d_, int e_, const std::vector<double>& q_, int lmax, const min_hamming_distance_class& f_)
  //    : d(d_), e(e_), s(dna_to_number(seed)), q(q_), f(f_)
  {
    init(seed, d_, e_, q_, lmax, f_);
  }
  
  void
  init(const std::string& seed, int d_, int e_, const std::vector<double>& q_, int lmax, const min_hamming_distance_class& f_)
  {
    d = d_;
    e = e_;
    s = dna_to_number<code_t>(seed);   // not needed
    q = q_;
    f = f_;
    k = seed.length();
    epsilon = k - cluster_threshold;
    imax = lmax - (k-epsilon);
    number_of_sequences = pow(4, k);
    a.resize(boost::extents[number_of_sequences][range(k, imax+1)][k-epsilon+1]);
    mask = 1;
    mask = (mask << (2*k)) - 1;

    initialize();
    compute();
  }

  double
  get(code_t code, int i, int hamdist) const
  {
    return a[code][i][hamdist];
  }
  
private:

  void
  initialize()
  {
    int number_of_sequences = pow(4, k);
    for (code_t code=0; code < number_of_sequences; ++code) {
      if (f.get(code) > d)
	a[code][k][1] = compute_bernoulli_probability(code, k, q);
    }
  }

  void
  compute()
  {
    for (int i=k; i < imax; ++i) {
      if (i < (2*k-epsilon)) {
	int hampos = i-k+1;
	for (code_t code=0; code < number_of_sequences; ++code) {
	  if (f.get(code) > d) {
	    double p = a[code][i][hampos];
	    if (p > 0.0) {
	      for (int x=0; x < 4; ++x) {
		code_t newcode = ((code << 2) + x) & mask;
		int j = f.get(newcode) <= d ? 0 : hampos + 1;
		if (j <= k-epsilon)
		  a[newcode][i+1][j] += p * q[x];
	      }
	    }
	  }
	}
      }
      else {
	for (int hampos=0; hampos <= k-epsilon; ++hampos) {
	  for (code_t code=0; code < number_of_sequences; ++code) {
	    double p = a[code][i][hampos];
	    if (p > 0.0 and f.get(code) > e) {
	      for (int x=0; x < 4; ++x) {
		code_t newcode = ((code << 2) + x) & mask;
		int j = f.get(newcode) <= d ? 0 : hampos + 1;
		if (j <= k-epsilon)
		  a[newcode][i+1][j] += p * q[x];
	      }
	    }
	  }
	}
      }
    }
  }



  int d;
  int e;
  code_t s;
  std::vector<double> q;
  min_hamming_distance_class f;
  int k;
  int imax;
  int epsilon;
  int number_of_sequences;
  code_t mask;
  boost::multi_array<double, 3> a;
};

class cluster_probability_type
{
public:

  cluster_probability_type() {}
  
  cluster_probability_type(const std::string& seed, int d, int e, const std::vector<double>& q_, int epsilon_, int max_cluster_len,
			   const min_hamming_distance_class& f, const min_hamming_distance_class& f_rev)
  //    : k(seed.length()), q(q_), epsilon(epsilon_), lmax(max_cluster_len),
  //      prefix(seed, d, e, q, lmax, f), suffix(reverse(seed), d, e, q, lmax, f_rev)
  {
    init(seed, d, e, q_, epsilon_, max_cluster_len, f, f_rev);
  }

  void
  init(const std::string& seed, int d, int e, const std::vector<double>& q_, int epsilon_, int max_cluster_len,
			   const min_hamming_distance_class& f, const min_hamming_distance_class& f_rev)
  {
    k = seed.length();
    q = q_;
    epsilon = epsilon_;
    lmax = max_cluster_len;
    prefix.init(seed, d, e, q, lmax, f);
    suffix.init(reverse(seed), d, e, q, lmax, f_rev);

  }
  
  double
  operator()(const std::string& u, int cluster_len)
  {
    code_t code = dna_to_number<code_t>(u);
    code_t code_rev = reverse_2bitstring(code, k);
    double p = 0.0;
    double div = compute_bernoulli_probability(code, k, q);
    for (int l1=2*k-epsilon; l1 <= cluster_len - (k-epsilon); ++l1) {
      p += prefix.get(code, l1, 0) * suffix.get(code_rev, cluster_len - l1 + k, 0) / div;
    }
    
    return p;
  }
  
private:
  int k;
  std::vector<double> q;
  int epsilon;
  int lmax;
  cluster_probability_array prefix;
  cluster_probability_array suffix;
};



// In the neighbourhood only the centermost hit has Hamming distance d to the seed, other hits
// must have Hamming distance higher that d.
class neighbourhood_probability_type
{
public:
  
  void
  //  init(const std::string& seed, int d, const std::vector<double>& q_, int epsilon_)
  init(const std::string& seed, int d, const dmatrix& m, int g)
  {
    k = seed.length();
    l = k + g;
    prefix.resize(pow(4, l));
    suffix.resize(pow(4, l));
    result.resize(pow(4, k));
    f.init(seed);  // Hamming distances between all k-mers and the seed
    dmatrix m_prefix = matrix_prefix(m, l);
    dmatrix m_suffix = matrix_suffix(m, g);
    //write_matrix(stdout, m_prefix, to_string("m_prefix:\n"), "%f");
    //write_matrix(stdout, m_suffix, to_string("m_suffix:\n"), "%f");
    assert(m_prefix.get_columns() == l);
    assert(m_suffix.get_columns() == g);
    
    code_t number_of_partial_sequences = pow(4, l);
    code_t maskk = ((code_t)1 << (2*k)) - 1;
    // Compute the prefix array
    for (code_t code=0; code < number_of_partial_sequences; ++code) {
      //      double p = compute_bernoulli_probability(code, l, q);
      double p = compute_normal_probability(code, l, m_prefix);
      int hd = f.get(code & maskk);
      for (int shift=1; shift <= g; ++shift) {
	if (f.get((code >> (2*shift)) & maskk) <= hd) {
	  p = 0.0;
	  break;
	}
      }
      prefix[code] = p;
    }

    // Compute the suffix array
    for (code_t code=0; code < number_of_partial_sequences; ++code) {
      //      double p = compute_bernoulli_probability(code, k-epsilon, q);  // Note! The center part is not included in the probability
      double p = compute_normal_probability(code, g, m_suffix);  // Note! The center part is not included in the probability
      int hd = f.get((code >> (2*g)) & maskk);
      for (int shift=0; shift < g; ++shift) {
	if (f.get((code >> (2*shift)) & maskk) <= hd) {
	  p = 0.0;
	  break;
	}
      }
      suffix[code] = p;
    }

    // Combine prefix and suffix
    code_t maskl = ((code_t)1 << (2*l)) - 1;   // bits corresponding to right-most l character positions are set to 1
    code_t number_of_full_sequences = pow(4, k+2*g);  
    for (code_t code=0; code < number_of_full_sequences; ++code) {
      code_t first = code >> (2*g);
      code_t second = code & maskl;
      code_t newcode = second >> (2*g);
      result[newcode] += prefix[first] * suffix[second];
    }
  }

  double
  operator()(const std::string& u) const
  {
    code_t code = dna_to_number<code_t>(u);
    return result[code];
  }
  
private:
  int k;
  int l;
  std::vector<double> prefix;
  std::vector<double> suffix;
  std::vector<double> result;
  min_hamming_distance_class f;
};

bool
test_esko_distance(const esko_distance_type& esko_distance, int k)
{
  std::string s(k, 'A');
  std::string r=s;
  code_t code_s = dna_to_number<code_t>(s);
  for (int i=0; i <= k; ++i) {
    if (i > 0)
      r[i-1] = 'C';
    code_t code_r = dna_to_number<code_t>(r);
    assert(esko_distance.distance(s,r) == i);
    assert(esko_distance.distance_with_bits(code_s, code_r) == i);
  }

  return true;
}

boost::multi_array<double, 1>
compute_conditional_neighbour_probabilities(int hamming_radius,
					    const prefix_hamming_distance_cumulative& e,
					    const prefix_hamming_distance_cumulative& e_reverse,
					    const std::string& seed,
					    const dmatrix& m,
					    const string_to_tuple_type& string_to_contributions)
{
  int k = seed.length();
  std::string pattern;
  boost::multi_array<double, 1> c;
  boost::multi_array<double, 1> p;
  typedef boost::multi_array<double, 1>::extent_range range;
  p.resize(boost::extents[range(-(k-1), k)]);
  c.resize(boost::extents[range(-(k-1), k)]);

  // p[0] is the probability of getting a sequence belonging to the hamming neighbourhood of the seed.
  // p[j] is the probability of getting a sequence of length k+|j| whose prefix and suffix of length k
  // are in the hamming neighbourhood of the seed.
  BOOST_FOREACH(boost::tie(pattern, boost::tuples::ignore), string_to_contributions) {
    double pm = compute_normal_probability(pattern, m);
    p[0] += pm;
    for (int j=-(k-1); j <= -1; ++j) {
      int h = hamming_distance(pattern.substr(0, k+j), seed.substr(-j, k+j)); // pattern prefix and seed suffix
      if (hamming_radius >= h)
	p[j] += pm * e.get(-j, hamming_radius - h);
    }
    for (int j=1; j <= k-1; ++j) {
      int h = hamming_distance(pattern.substr(j, k-j), seed.substr(0, k-j));  // pattern suffix and seed prefix
      if (hamming_radius >= h)
	p[j] += pm * e_reverse.get(j, hamming_radius - h);
    }
    
  }
  // c[j] is the probability of getting a hamming occurrence in position j on the condition that
  // there is a hamming occurrence in the position 0
  for (int j=-(k-1); j <= k-1; ++j) {
    c[j] = p[j] / p[0]; 
    //    c_sum += c[j];
  }
  return c;
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


// Multinomial-n pfm when data (single random sequence) is distributed according to Bernoulli distribution 'bg' 
dmatrix
get_expected_model_in_bg(const string_to_tuple_type& string_to_contributions, int k, const std::vector<double>& bg)
{
  dmatrix result(4, k);
  std::string pattern;
  std::vector<boost::tuple<int, int> > pairs;
  iupac_probability_in_background iupac_prob(bg);
  BOOST_FOREACH(boost::tie(pattern, pairs), string_to_contributions) {
    double p = iupac_prob(pattern);
	
    int j, a;
    BOOST_FOREACH(boost::tie(j,a), pairs) {
      result(a, j) += p;
    }
  }
  return result;
}

// Multinomial-n pfm when data (single random sequence) is distributed according to another pfm
dmatrix
get_expected_model_in_pwm(const string_to_tuple_type& string_to_contributions, int k, const dmatrix& m)
{
  dmatrix result(4, k);
  std::string pattern;
  std::vector<boost::tuple<int, int> > pairs;
  //  iupac_probability_in_background iupac_prob(bg);
  BOOST_FOREACH(boost::tie(pattern, pairs), string_to_contributions) {
    double p = compute_normal_probability(pattern, m);
	
    int j, a;
    BOOST_FOREACH(boost::tie(j,a), pairs) {
      result(a, j) += p;
    }
  }
  return result;
}



// If the pwm m is positioned in positions [0,k) and other positions are distributed according to bg distribution q,
// then return the matrix from positions [j,j+k).
dmatrix
get_shifted_window_pwm(const dmatrix& m, const dvector& q, int j)
{
  int k = m.get_columns();
  assert(-(k-1) <= j);
  assert(j <= k-1);
  dmatrix result(4, k);
  
  if (j < 0) {
    for (int i=0; i < -j; ++i)
      result.set_column(i, q);
    for (int i=-j; i < k; ++i)
      result.set_column(i, m.column(i+j));
  }
  else {
    for (int i=0; i < k-j; ++i)
      result.set_column(i, m.column(i+j));
    for (int i=k-j; i < k; ++i)
      result.set_column(i, q);
  }
  return result;
}

boost::tuple<double, double>
estimate_neighbour_count(const boost::multi_array<dmatrix,1>& corrected_shift_matrices, double lambda, int N,
			 const std::string& u, const std::vector<double>& q)
{
  int k = corrected_shift_matrices[0].get_columns();
  //  double bg_count = (1.0 - lambda*(2*k-1))*N*pow(4,-k);
  double bg_count = (1.0 - lambda*(2*k-1))*N * compute_bernoulli_probability(u, q);
  double p = 0.0;
  for (int j=-(k-1); j <= k-1; ++j) {
    p += compute_normal_probability(u, corrected_shift_matrices[j]);
  }
  double pfm_count = p*lambda*N;
  return boost::make_tuple(pfm_count, bg_count);
}

int
sign(double x)
{
  if (x < 0.0)
    return -1;
  else
    return 1;
}

std::vector<boost::tuple<double, double> >
find_roots(const std::vector<double>& counts, const std::vector<double>& r, double value)
{
  std::vector<boost::tuple<double, double> > solutions;
  int len = counts.size();
  assert(len == r.size());
  for (int i=0; i < len-1; ++i) {
    if (sign(counts[i]-value) != sign(counts[i+1]-value))
      solutions.push_back(boost::make_tuple(r[i], r[i+1]));
  }
  
  return solutions;
}

dmatrix
correct_pfm(const dmatrix& pfm, const dmatrix& bg_pfm, const boost::multi_array<dmatrix,1>& shifted_bg_pfms,
	    int N, double lambda)
{
  int k = pfm.get_columns();
  int g=cluster_threshold; // this is the maximum shift to be considered as neighbour
  //  dmatrix corrected = pfm - (1.0 - (2*g+1)*lambda)*N*bg_pfm;
  dmatrix full_bg_pfm = (1.0 - lambda*(2*(k+g)-1))*N*bg_pfm;
  for (int j=-(k+g-1); j <= k+g-1; ++j) {
    if (j == 0)
	  continue;
    full_bg_pfm += lambda * N * shifted_bg_pfms[j];
  }
  dmatrix corrected = pfm - full_bg_pfm;
  corrected.apply(cut);  // cut negative values to 0
  return normalize_matrix_columns_copy(corrected);
}

std::vector<double>
get_counts(const std::string& u, int observed_count, int N,
	   const std::vector<boost::multi_array<dmatrix,1> >& corrected_pfms,
	   const std::vector<double>& r, const std::vector<double>& q)
{
  std::vector<double> counts_total;
  int len = corrected_pfms.size();
  assert(len == r.size());
  for (int i=0; i < len; ++i) {
    double lambda = r[i];
    double count_signal, count_bg;
    boost::tie(count_signal, count_bg) = estimate_neighbour_count(corrected_pfms[i], lambda, N, u, q);
    counts_total.push_back(count_signal + count_bg);
  }
  return counts_total;
}


double 
estimate_signal_fraction(const dmatrix& pfm, dmatrix bg_pfm, const boost::multi_array<dmatrix,1>& shifted_bg_pfms,
			 const std::vector<std::string>& sequences,
			 const std::vector<double>& q)
{
  int n = sequences.size();
  int L = sequences[0].length();
  int k = pfm.get_columns();
  int N = n * (L-k+1);   // Number of sites of type 'all'
  int g = cluster_threshold;
                                                          
  std::vector<boost::multi_array<dmatrix,1> > corrected_pfms; // a pfm for each lambda and j (shift)
  //  dvector q(4, 0.25);  // uniform background
  dvector r;
  //  double max_lambda = 1.0;
  //  double max_lambda = 1.0 / (2*k-1);
  double max_lambda = 1.0 / (2*(k+g)-1);
  int number_of_intervals = 100;
  double step = max_lambda / (double)number_of_intervals;
  for (int i = 0; i <= number_of_intervals ; ++i) {
    r.push_back(i*step);
  }
  int len = r.size();
  for (int i=0; i < len; ++i) {     // Iterate through different lambda values
    double lambda = r[i];
    dmatrix corrected_pfm = correct_pfm(pfm, bg_pfm, shifted_bg_pfms, N, lambda);
    //typedef boost::multi_array<dmatrix, 1>::extent_range range;
    typedef boost::multi_array_types::extent_range range;

    boost::multi_array<dmatrix,1> shifted_pfms(boost::extents[range(-(k-1), k)]);
    for (int j=-(k-1); j <= k-1; ++j) {
      shifted_pfms[j] = get_shifted_window_pwm(corrected_pfm, q, j);
    }
    corrected_pfms.push_back(shifted_pfms);
  }

  typedef boost::tuple<code_t, int> elem_t;
  std::vector<elem_t> v; // kmers sorted by counts
  count_container_t counts;
  get_kmer_counts(sequences, k, counts, use_two_strands);
  code_t code;
  int count;
  BOOST_FOREACH(boost::tie(code, count), counts)
    v.push_back(boost::make_tuple(code, count));

  std::sort(v.begin(), v.end(), [](elem_t x, elem_t y) -> bool { return x.get<1>() > y.get<1>(); });   // sort key is the count

  std::string u;
  int seq_count;
  int solution_type = 0;  // 0 = no solution found yet, 1 = found unique solution, 2 = found two solutions
  BOOST_FOREACH(boost::tie(code, count), v) {
    std::string seq = number_to_dna(code, k);
    std::vector<double> counts_total = get_counts(seq, count, N, corrected_pfms, r, q);  // one count for each lambda
    std::vector<boost::tuple<double, double> > solutions = find_roots(counts_total, r, count);
    if (solutions.size() == 1) {   // we will prefer this
      solution_type = 1;
      u = seq;
      seq_count = count;
      break;
    }
    else if (solutions.size() == 2 and solution_type == 0) {  // but will settle for this if no better solution exists
      solution_type = 2;
      u = seq;
      seq_count = count;
    }
      
  }

  if (solution_type == 0)    // Fail-safe
    return 0.0;              // Assume everything is background
  
  std::vector<double> counts_total = get_counts(u, seq_count, N, corrected_pfms, r, q);
  std::vector<boost::tuple<double, double> > solutions = find_roots(counts_total, r, seq_count);

  double a,b;
  boost::tie(a, b) = solutions[0];  // choose the lowest solution if multiple exist
  return (a+b)/2.0;                 // return the average of end points as solution 
}



// This uses some initial value for 'm' to obtain count of occurrences.
// And iteratively uses that count to obtain better estimate for 'm'.
// In practise this doesn't work at all.
// Occurrences are defined to be all sequences within Hamming-radius from the seed.
double
estimate_real_count_version_2(int observed_count, const std::string& seed, const dmatrix& m, int hamming_radius, int L,
			      const string_to_tuple_type& string_to_contributions, const std::vector<double>& bg, double gamma)
{
  int k = seed.length();
  boost::multi_array<double, 1> c;
  //typedef boost::multi_array<double, 1>::extent_range range;
  typedef boost::multi_array_types::extent_range range;
  c.resize(boost::extents[range(-(k-1), k)]);
  std::string pattern;
  int max_hd = 5;   // At least this should cover most instances of a PWM, when the seed is the same as the consensus of the PWM
  std::vector<string_to_tuple_type> hamming_neighbourhood(max_hd+1);
  for (int H=0; H <= max_hd; ++H)
    hamming_neighbourhood[H] = get_n_neighbourhood_mononucleotide_contributions(seed, H);
  prefix_hamming_distance_cumulative e(seed, bg);
  prefix_hamming_distance_cumulative e_reverse(reverse(seed), bg);
  dvector d(max_hd+1, 0.0);
  BOOST_FOREACH(boost::tie(pattern, boost::tuples::ignore), string_to_contributions) {
    int hd = hamming_distance(pattern, seed);
    double p = compute_bernoulli_probability(pattern, bg);
    for (; hd <= max_hd; ++hd)
      d[hd] += p;
  }

  dmatrix m_bg = get_expected_model_in_bg(string_to_contributions, k, bg);
  double x = 0.0;
  dmatrix m_new=normalize_matrix_columns_copy(m);
  std::map<int, dmatrix> m_shifted;
  for (int H=1; H <= max_hd; ++H) {
    for (int round=0; round < 20; ++round) {
      for (int j=-(k-1); j <= k-1; ++j) {
	m_shifted[j] = get_expected_model_in_pwm(string_to_contributions, k, get_shifted_window_pwm(m_new, bg, j));
	//write_matrix(stdout, x*m_shifted[j], to_string("m_shifted %i\n", j));
      }
      m_new = m - ((L - x*(2*k-1)) * m_bg);
      for (int j=-(k-1); j <= k-1; ++j) {
	if (j==0)
	  continue;
	m_new -= x*m_shifted[j];
      }
      write_matrix(stdout, m, "m\n");
      write_matrix(stdout, (L - x*(2*k-1)) * m_bg, "times bg\n");
      write_matrix(stdout, m_new, "m_new\n");
      m_new.apply(cut);
      normalize_matrix_columns(m_new);
      write_matrix(stdout, m_new, "m_new normalized\n");
      print_ics(m_new);
      c = compute_conditional_neighbour_probabilities(H, e, e_reverse,
						      seed, m_new, hamming_neighbourhood[H]);
      double c_sum = sum(c);
      double a = observed_count - L*d[H];
      double b = c_sum - (2*k-1)*d[H];
      if (a < gamma or b < gamma)   // So that the division doesn't go unstable nor any component goes negative
	goto myexit;
      x = a / b;
      //x = std::min(x, (double)L/(2*k-1));
      printf("round=%i H=%i a=%f b=%f x=%f\n", round, H, a, b, x);
    }
  }
 myexit:
  
  return x;
}


// This uses some initial value for 'm' to obtain count of occurrences.
// Occurrences are defined to be all sequences within Hamming-radius from the seed.
double
estimate_real_count(int observed_count, const std::string& seed, const dmatrix& m, int hamming_radius, int L,
		    const string_to_tuple_type& string_to_contributions, const std::vector<double>& bg)
{
  int k = seed.length();
  //  boost::multi_array<double, 1> p;
  boost::multi_array<double, 1> c;
  // typedef boost::multi_array<double, 1>::extent_range range;
  // p.resize(boost::extents[range(-(k-1), k)]);
  // c.resize(boost::extents[range(-(k-1), k)]);
  std::string pattern;
  prefix_hamming_distance_cumulative e(seed, bg);
  prefix_hamming_distance_cumulative e_reverse(reverse(seed), bg);
  printf("prefix hamming distance cumulative\n");
  e.print();
  printf("prefix hamming distance cumulative reverse\n");
  e_reverse.print();
  
  double d = 0.0;
  BOOST_FOREACH(boost::tie(pattern, boost::tuples::ignore), string_to_contributions) {
    d += compute_bernoulli_probability(pattern, bg);
  }

  c = compute_conditional_neighbour_probabilities(hamming_radius, e, e_reverse,
						  seed, m, string_to_contributions);
  double c_sum = sum(c);

  printf("d is %f\n", d);
  printf("c_sum is %f\n", c_sum);
  // printf("p array is:\n");
  // print_array_with_default_headers(stdout, p, "%.5f");
  printf("c array is:\n");
  print_array_with_default_headers(stdout, c, "%.5f");
  double upstairs = observed_count - L*d;
  double downstairs = c_sum -(2*k-1)*d;
  printf("Upstairs is %f\n", upstairs);
  printf("Downstairs is %f\n", downstairs);
	
  double x = upstairs / downstairs;

  return x;
}

dmatrix
constant_column_matrix(const std::vector<double>& q, int length)
{
  dmatrix result(4, length);
  for (int i=0; i < length; ++i)
    result.set_column(i, q);

  return result;
}

boost::tuple<dmatrix,double>
neighbour_expected_pfm_in_pfm_distribution(const std::string& seed,
					   const dmatrix& m,
					   int hamming_radius)
{
  int k = seed.length();
  int g = cluster_threshold;
  assert(m.get_columns() == k+2*g);
  string_to_tuple_type string_to_contributions;
  string_to_contributions = get_n_neighbourhood_mononucleotide_contributions(seed, hamming_radius);

  neighbourhood_probability_type neighbourhood_probability;
  neighbourhood_probability.init(seed, hamming_radius, m, g);
  dmatrix result(4, k);

  std::vector<boost::tuple<int, int> > pairs;
  std::string pattern;
  double prob_sum = 0.0;
  BOOST_FOREACH(boost::tie(pattern, pairs), string_to_contributions) {

    double p=0.0;
 
    p = neighbourhood_probability(pattern);
    if (use_two_strands)
      p += neighbourhood_probability(reverse_complement(pattern));
    prob_sum += p;
 
    int j, a;
    BOOST_FOREACH(boost::tie(j,a), pairs) {
      result(a, j) += p;
    }
  }
  return boost::make_tuple(result, prob_sum);
 //return result;
}


boost::tuple<dmatrix,double>
find_multinomial_n_background(const std::string& seed, const std::vector<std::string>& sequences, const std::vector<double>& bg,
			      int hamming_radius, bool use_multimer)
{
  //TIME_START(t);
  const int k = seed.length();
  const int L = sequences[0].length();
  int g = cluster_threshold;
  assert(hamming_radius >= 0);
  assert(hamming_radius <= k);
  dmatrix result(4, k);
  esko_distance_type esko_distance(k);
  test_esko_distance(esko_distance, k);
  printf("Esko-distance weights are %s\n", print_vector(esko_distance.get_weights()).c_str());
  printf("Using background distribution %s\n", print_vector(bg).c_str());
  string_to_tuple_type string_to_contributions;
  string_to_contributions = get_n_neighbourhood_mononucleotide_contributions(seed, hamming_radius);

  std::vector<boost::tuple<int, int> > pairs;
  std::string pattern;

  //printf("Number of patterns %lu\n", string_to_contributions.size());
  //double seed_count=0;
  int epsilon = k - cluster_threshold;
  double total_count=0;
  int lines = sequences.size();
  int number_of_orientations = use_two_strands ? 2 : 1;

  iupac_probability_in_background iupac_prob(bg);
  //  int max_cluster_len = k + 4*(k-epsilon);  // This is crude approximation.
  int max_cluster_len = 200;  // This is crude approximation.
  int min_cluster_len = k + 2*(k-epsilon);  // By definition of cluster, the cluster length cannot be shorter than this

  std::string seed_rev = reverse_complement(seed);
    
  double prob_sum = 0.0;
  min_hamming_distance_class f;
  min_hamming_distance_class f_rev;
  neighbourhood_probability_type neighbourhood_probability;
  if (background_counting == choose_one_per_cluster) {
    f.init(seed);
    f_rev.init(reverse(seed));
  }
  else if (background_counting == neighbourhood_contains_one) {
    dmatrix m = constant_column_matrix(bg, k+2*g);
    neighbourhood_probability.init(seed, hamming_radius, m, g);
  }

  if (background_counting == choose_one_per_cluster) {
    std::vector<string_to_tuple_type> hamming_neighbours(100*hamming_radius+1);   // Bin the patterns according to Hamming distance to the seed
    BOOST_FOREACH(boost::tie(pattern, pairs), string_to_contributions) {
      hamming_neighbours[esko_distance.distance(seed, pattern)].insert(std::make_pair(pattern, pairs));
    }
    for (int e=0; e <= 100*hamming_radius; ++e) {
      if (hamming_neighbours[e].size() == 0)
	continue;
      cluster_probability_type cluster_probability;
      cluster_probability.init(seed, 100*hamming_radius, e, bg, epsilon, max_cluster_len, f, f_rev);
      BOOST_FOREACH(boost::tie(pattern, pairs), hamming_neighbours[e]) {

	double p;
	double count=0;
	for (int cluster_len=min_cluster_len; cluster_len <= max_cluster_len; ++cluster_len) {
	  int cluster_sites = lines * (L-cluster_len+1) * number_of_orientations;
	  p = cluster_probability(pattern, cluster_len);
	  prob_sum += p;
	  count += cluster_sites * p;
	}
	total_count += count;

	int j, a;
	BOOST_FOREACH(boost::tie(j,a), pairs) {
	  result(a, j) += count;
	}
      }
    }
  }
  else {  
    BOOST_FOREACH(boost::tie(pattern, pairs), string_to_contributions) {

      double p=0.0;
      switch (background_counting) {
      case choose_one_per_cluster:
	break;
      case neighbourhood_contains_one:
	p = neighbourhood_probability(pattern);
	if (use_two_strands)
	  p += neighbourhood_probability(reverse_complement(pattern));
	prob_sum += p;
	break;
      case sequence_contains_at_least_one:
	error(true, "Not implemented.");
	break;
      case all_occurrences:
	p = iupac_prob(pattern);
	if (use_two_strands)
	  p += iupac_prob(reverse_complement(pattern));
	prob_sum += p;
	break;
      case sequence_contains_one:
	error(true, "Not implemented.");
	break;
      }
      //	if (not iupac_string_match(pattern, seed))
	
      int j, a;
      BOOST_FOREACH(boost::tie(j,a), pairs) {
	result(a, j) += p;
      }
    }
  }
  
  printf("Total probability of hamming neighbourhood in background is %f\n", prob_sum);
  //TIME_PRINT("find_multinomial_n_background took %.2f seconds.\n", t);
  return boost::make_tuple(result, prob_sum);
} // find_multinomial_n_background




boost::tuple<dmatrix, unsigned long, unsigned long>
count_all_occurrences(const string_to_tuple_type& string_to_contributions, const std::string& seed, const suffix_array& sa,
		      const std::vector<std::string>& sequences)
{
  typedef string_to_tuple_type::const_iterator iterator;
  unsigned long seed_count = 0;
  unsigned long total_count = 0;
  int k = seed.length();
  int lines = sequences.size();
  int L = sequences[0].length();
  std::vector<std::vector<int> > hit_positions(lines);
  std::vector<std::vector<int> > hit_directions(lines);
  std::vector<std::vector<iterator> > hit_contributions(lines);
  int divisor = L + 1;    // Includes the separator '#'
  
  dmatrix result(4, k);
  std::vector<boost::tuple<int, int> > pairs;
  std::string pattern;
  if (print_alignment) 
    printf("#String\tColumn\tHamming distance\tPalindrome\tCount\tMatches at col\n");
  //  BOOST_FOREACH(boost::tie(pattern, pairs), string_to_contributions) {

  // All the following hassle is just to categorize the hits by the sequence they appear in,
  // and to make sure that palindromes are counted correctly.
  for (iterator it=string_to_contributions.begin(); it != string_to_contributions.end(); ++it) {
    pattern = it->first;
    pairs = it->second;
    
    std::vector<long int> positions;
    sa.locate_iupac(pattern, positions);
    BOOST_FOREACH(int pos, positions) {
      int i = pos / divisor;  // index of the read containing the pos
      int j = pos % divisor;
      int dir = 1;
      if (i >= lines) {       // handle reverse complement
	i = 2*lines - i - 1;
	j = L - (j + k - 1) - 1;
	dir = -1;
      }
      bool is_palindrome = is_palindromic(sequences[i].substr(j, k));
      if (not is_palindrome or count_palindromes_twice) { // This branch allows the same site to be counted twice, for both orientations
	hit_positions[i].push_back(j);
	hit_directions[i].push_back(dir);
	hit_contributions[i].push_back(it);
      }
      else {  // count palindrome only once
	if (std::find(hit_positions[i].begin(), hit_positions[i].end(), j) == hit_positions[i].end()) {
	  hit_positions[i].push_back(j);
	  hit_directions[i].push_back(dir);
	  hit_contributions[i].push_back(it);
	}
      }
    } // end foreach pos

  } // end for it
  
  for (int i=0; i < lines; ++i) {
    int hit_count = hit_positions[i].size();
    
    // This sorting is just for the printing of the alignment
    std::vector<int> hits(hit_count);
    for (int index=0; index < hit_count; ++index) {
      hits[index] = index;
    }
    std::sort(hits.begin(), hits.end(), [i, &hit_positions](int a, int b) { return hit_positions[i][a] < hit_positions[i][b];} );

    BOOST_FOREACH(int index, hits) {
      boost::tie(pattern, pairs) = *(hit_contributions[i][index]);
      //      int start_pos = hit_positions[i][index];
      //      int dir = hit_directions[i][index];
      std::string pal = is_palindromic(pattern) ? "Palindrome" : "-";
      int hd = hamming_distance(seed, pattern);
      //printf("*%i\t%i\t%i\t%i\t%s\t%s\t%i\n", i, start_pos, dir, k, pattern.c_str(), pal.c_str(), hd);
      

      if (iupac_string_match(pattern, seed))
	seed_count += 1;
      total_count += 1;
      int j, a;
      BOOST_FOREACH(boost::tie(j,a), pairs) {
	if (print_alignment) 
	  printf("#%s\t%i\t%i\t%s\t%zu\t%s\n", pattern.c_str(), j, hd,
		 yesno(is_palindromic(pattern)), (size_t)1, yesno(seed[j]==pattern[j]));

	result(a, j) += 1;
      }
    } // end for index
  } // end i

  return boost::make_tuple(result, seed_count, total_count);
};


// Note! This assumes that all sequences are of equal length, for efficiency
boost::tuple<dmatrix, unsigned long, unsigned long>
count_sequence_contains_one(const string_to_tuple_type& string_to_contributions, const std::string& seed, const suffix_array& sa,
			    const std::vector<std::string>& sequences)
{
  // allow reads that contain exactly one hit.
  int lines = sequences.size();

  int L = sequences[0].length();
  for (int i=0; i < lines; ++i)
    assert(sequences[i].length() == L);
  
  unsigned long seed_count = 0;
  unsigned long total_count = 0;
  int k = seed.length();
  
  dmatrix result(4, k);
  typedef string_to_tuple_type::const_iterator iterator;
  std::vector<std::set<int> > hit_positions(lines);
  std::vector<std::vector<iterator> > hit_patterns(lines);
  int divisor = L + 1;    // Includes the separator '#'
  std::vector<boost::tuple<int, int> > pairs;
  std::string pattern;
  for (iterator it=string_to_contributions.begin(); it != string_to_contributions.end(); ++it) {
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
  
		    
  return boost::make_tuple(result, seed_count, total_count);
};

boost::tuple<dmatrix, unsigned long, unsigned long>
count_choose_one_per_cluster(const string_to_tuple_type& string_to_contributions, const std::string& seed, const suffix_array& sa,
			     const std::vector<std::string>& sequences)
{
  typedef string_to_tuple_type::const_iterator iterator;
  int lines = sequences.size();
  int L = sequences[0].length();
  for (int i=0; i < lines; ++i)
    assert(sequences[i].length() == L);
  std::string seed_rev = reverse_complement(seed);
  
  unsigned long seed_count = 0;
  unsigned long total_count = 0;
  int k = seed.length();
  int epsilon = k - cluster_threshold;
  esko_distance_type esko_distance(k);
  
  dmatrix result(4, k);
  std::vector<std::vector<int> > hit_positions(lines);
  std::vector<std::vector<int> > hit_directions(lines);
  std::vector<std::vector<iterator> > hit_contributions(lines);
  int divisor = L + 1;    // Includes the separator '#'
  std::vector<boost::tuple<int, int> > pairs;
  std::string pattern;

  // Bin the occurrences according to the sequence they are in.
  for (iterator it=string_to_contributions.begin(); it != string_to_contributions.end(); ++it) {  // iterator through string in Hamming neighbourhood
    std::vector<long int> positions;
    pattern = it->first;
    bool is_palindrome = is_palindromic(pattern);
    sa.locate_iupac(pattern, positions);
    BOOST_FOREACH(int pos, positions) {
      int i = pos / divisor;  // index of the read containing the pos
      int j = pos % divisor;
      int dir = 1;
      if (i >= lines) {       // handle reverse complement
	i = 2*lines - i - 1;
	j = L - (j + k - 1) - 1;
	dir = -1;
      }
      if (not is_palindrome or count_palindromes_twice) { // This branch allows the same site to be counted twice, for both orientations
	hit_positions[i].push_back(j);
	hit_directions[i].push_back(dir);
	hit_contributions[i].push_back(it);
      }
      else {  // count palindrome only once
	if (std::find(hit_positions[i].begin(), hit_positions[i].end(), j) == hit_positions[i].end()) {
	  hit_positions[i].push_back(j);
	  hit_directions[i].push_back(dir);
	  hit_contributions[i].push_back(it);
	}
      }
    }
  }
  
  for (int i=0; i < lines; ++i) {
    int hit_count = hit_positions[i].size();
    std::vector<int> best_from_each_cluster;
    if (hit_count <= 1)
      best_from_each_cluster.push_back(0);
    else {

      ////////////////////////
      // Form the clusters

      int max_cluster_size=0;
      int max_cluster_length=0;
      std::vector<std::vector<int> > clusters;
      std::vector<int> hits;               // This will contain indices to hit_positions/patterns vector
      for (int index=0; index < hit_count; ++index)
	hits.push_back(index);
      std::sort(hits.begin(), hits.end(), [i, &hit_positions](int a, int b) { return hit_positions[i][a] < hit_positions[i][b];} );
      std::vector<int> new_cluster;
      int previous_end_position=sequences[i].length();  // Sure to give an overlap in the following test
      //int distance_threshold = -1000;  // each occurrence will be considered its own cluster
      for (int t=0; t < hit_count; ++t) {
	if (previous_end_position - hit_positions[i][hits[t]] + 1 >= epsilon) {   // do two consequent hits overlap?
	  new_cluster.push_back(hits[t]);           // extend the old cluster
	  previous_end_position = hit_positions[i][hits[t]] + k - 1;
	} else {
	  clusters.push_back(new_cluster);   // store the old cluster
	  new_cluster.clear();               // Start a new cluster
	  new_cluster.push_back(hits[t]);    // new cluster now contains the current hit
	  previous_end_position = hit_positions[i][hits[t]] + k - 1;
	}
      }
      clusters.push_back(new_cluster);  // last cluster

      std::vector<int> cluster;         // cluster contains indices to hit_contributions vector
      BOOST_FOREACH(cluster, clusters) {
	int size = cluster.size();
	if (size > max_cluster_size)
	  max_cluster_size = size;
	int len = hit_positions[i][cluster[size-1]] + k - hit_positions[i][cluster[0]];
	if (len > max_cluster_length)
	  max_cluster_length = len;
      }
	  
      ////////////////////////
      // Iterate the clusters

      std::vector<int> cluster_size_distribution(max_cluster_size+1);
      std::vector<int> cluster_length_distribution(max_cluster_length+1);
      BOOST_FOREACH(cluster, clusters) {
	int size = cluster.size();
	int len = hit_positions[i][cluster[size-1]] + k - hit_positions[i][cluster[0]];
	++cluster_size_distribution[cluster.size()];
	++cluster_length_distribution[len];
	if (cluster.size() == 1)
	  best_from_each_cluster.push_back(cluster[0]);
	else {
	  std::map<int, int> hds;  // map from index to Hamming distance
	  if (use_two_strands) 
	    BOOST_FOREACH(int index, cluster)
	      hds[index] = std::min(esko_distance.distance(seed, hit_contributions[i][index]->first),
				  esko_distance.distance(seed_rev, hit_contributions[i][index]->first));
	  else
	    BOOST_FOREACH(int index, cluster)
	      hds[index] = esko_distance.distance(seed, hit_contributions[i][index]->first);

	    
	  // sort cluster by Hamming distance of occurrences to the seed
	  std::sort(cluster.begin(), cluster.end(), [i, &hit_positions, &hds](int a, int b) { 
	      // a and b are indices to hit_positions vector

	      // Sort primarily by Hamming index and secondarily by start position
	      return hds[a] < hds[b] or (hds[a] == hds[b] and hit_positions[i][a] < hit_positions[i][b]);
	    } );
	    
	  int best_hd = hds[cluster[0]];
	  int best_hd_count = 0;

	  for (int t = 0; t < cluster.size() and hds[cluster[t]] == best_hd; ++t)
	    ++best_hd_count;

	  if (best_hd_count <= 1)
	    best_from_each_cluster.push_back(cluster[0]); // Add the element (index) of cluster with smallest Hamming distance
	  /*
	  else {
	    int best_theoretical_starting_point = floor((hit_positions[i][cluster[best_hd_count-1]] - hit_positions[i][cluster[0]])/2); // center point
	    int min_dist=std::numeric_limits<int>::max();
	    int min_arg=0;
	    for (int t=0; t < best_hd_count; ++t) {
	      int dist = abs(best_theoretical_starting_point - hit_positions[i][cluster[t]]);
	      if (dist < min_dist) {
		min_dist = dist;
		min_arg = t;
	      }
	    }
	    best_from_each_cluster.push_back(cluster[min_arg]); // Add the element of cluster with smallest Hamming distance, if not unique
	    // use the centermost of those having the best Hamming distance
	  }
	  */
	      
	}
      }  // end BOOST_FOREACH cluster
      printf("#Cluster size\tCount\n");
      for (int i=0; i<= max_cluster_size; ++i) {
	printf("#%i\t%i\n", i, cluster_size_distribution[i]);
      }
	  
      printf("$Cluster length\tCount\n");
      for (int i=0; i<= max_cluster_length; ++i) {
	printf("$%i\t%i\n", i, cluster_length_distribution[i]);
      }
    }
    for (int t=0; t < best_from_each_cluster.size(); ++t) {
      int index = best_from_each_cluster[t];
      boost::tie(pattern, pairs) = *(hit_contributions[i][index]);
      int start_pos = hit_positions[i][index];
      int dir = hit_directions[i][index];
      std::string pal = is_palindromic(pattern) ? "Palindrome" : "-";
      int hd = hamming_distance(seed, pattern);
      printf("*%i\t%i\t%i\t%i\t%s\t%s\t%i\n", i, start_pos, dir, k, pattern.c_str(), pal.c_str(), hd);

      if (not iupac_string_match(pattern, seed))
	total_count += 1;
      else
	++seed_count;
      int j, a;
      BOOST_FOREACH(boost::tie(j,a), pairs) {
	result(a, j) += 1;
      }
    }
  }  // end i
      
  return boost::make_tuple(result, seed_count, total_count);
};


boost::tuple<dmatrix, unsigned long, unsigned long>
count_neighbourhood_contains_one(const string_to_tuple_type& string_to_contributions, const std::string& seed, const suffix_array& sa,
				  const std::vector<std::string>& sequences)
{
  typedef string_to_tuple_type::const_iterator iterator;
  int lines = sequences.size();
  int L = sequences[0].length();
  for (int i=0; i < lines; ++i)
    assert(sequences[i].length() == L);

  unsigned long seed_count = 0;
  unsigned long total_count = 0;
  int k = seed.length();
  std::string seed_rev = reverse_complement(seed);
  //  int epsilon = k - cluster_threshold;
  esko_distance_type esko_distance(k);
  dmatrix result(4, k);
  std::vector<std::vector<int> > hit_positions(lines);
  std::vector<std::vector<int> > hit_directions(lines);
  std::vector<std::vector<iterator> > hit_contributions(lines);
  int divisor = L + 1;    // Includes the separator '#'
  std::vector<boost::tuple<int, int> > pairs;
  std::string pattern;

  // Bin the occurrences according to the sequence they are in.
  for (iterator it=string_to_contributions.begin(); it != string_to_contributions.end(); ++it) {  // iterator through strings in Hamming neighbourhood
    std::vector<long int> positions;
    pattern = it->first;
    bool is_palindrome = is_palindromic(pattern);
    sa.locate_iupac(pattern, positions);
    BOOST_FOREACH(int pos, positions) {
      int i = pos / divisor;  // index of the read containing the pos
      int j = pos % divisor;
      int dir = 1;
      if (i >= lines) {       // handle reverse complement
	i = 2*lines - i - 1;
	j = L - (j + k - 1) - 1;
	dir = -1;
      }
      if (not is_palindrome or count_palindromes_twice) {  // This branch allows a site to be counted twice, for both orientations
	hit_positions[i].push_back(j);
	hit_directions[i].push_back(dir);
	hit_contributions[i].push_back(it);
      }
      else {  // count palindrome only once
	if (std::find(hit_positions[i].begin(), hit_positions[i].end(), j) == hit_positions[i].end()) {
	  hit_positions[i].push_back(j);
	  hit_directions[i].push_back(dir);
	  hit_contributions[i].push_back(it);
	}
      }
    }
  }

  
  for (int i=0; i < lines; ++i) {
    int hit_count = hit_positions[i].size();
    std::vector<int> non_intersecting_occurrences;
    if (hit_count == 1)
      non_intersecting_occurrences.push_back(0);
    else {


      //std::vector<std::vector<int> > clusters;
      std::vector<int> hits(hit_count);               // This will contain indices to hit_positions/patterns vector
      std::vector<int> hamming_distances(hit_count);
      for (int index=0; index < hit_count; ++index) {
	hits[index] = index;
	if (use_two_strands)
	  hamming_distances[index] = std::min(esko_distance.distance(hit_contributions[i][index]->first, seed),
					      esko_distance.distance(hit_contributions[i][index]->first, seed_rev));
	else
	  hamming_distances[index] = esko_distance.distance(hit_contributions[i][index]->first, seed);
      }
      std::sort(hits.begin(), hits.end(), [i, &hit_positions](int a, int b) { return hit_positions[i][a] < hit_positions[i][b];} );

      for (int current=0; current < hit_count; ++current) {
	int current_hd = hamming_distances[hits[current]];
	int current_pos = hit_positions[i][hits[current]];
	int prev_index = current - 1;
	bool stop = false;
	while (prev_index >= 0 and current_pos - hit_positions[i][hits[prev_index]] <= cluster_threshold) {
	  if (hamming_distances[hits[prev_index]] <= current_hd and hit_positions[i][hits[prev_index]] != current_pos) {
	    stop = true;
	    break;
	  }
	  --prev_index;
	}
	if (stop)
	  continue;
	int next_index = current + 1;
	while (next_index < hit_count and hit_positions[i][hits[next_index]] - current_pos <= cluster_threshold) {
	  if (hamming_distances[hits[next_index]] <= current_hd and hit_positions[i][hits[next_index]] != current_pos) {
	    stop = true;
	    break;
	  }
	  ++next_index;
	}
	if (stop)
	  continue;
	non_intersecting_occurrences.push_back(hits[current]);
      }

      /*				       

      ////////////////////////
      // Form the clusters, and store clusters with size one to non_intersecting_occurrences container

      std::vector<int> new_cluster;
      int previous_end_position=sequences[i].length();  // Sure to give an overlap in the following test
      //int distance_threshold = -1000;  // each occurrence will be considered its own cluster
      for (int t=0; t < hit_count; ++t) {
	if (previous_end_position - hit_positions[i][hits[t]] + 1 >= epsilon) {   // do two consequent hits overlap?
	  new_cluster.push_back(hits[t]);           // extend the old cluster
	  previous_end_position = hit_positions[i][hits[t]] + k - 1;
	} else {
	  if (new_cluster.size() == 1)         // store only if cluster has one element
	    clusters.push_back(new_cluster);   // store the old cluster
	  new_cluster.clear();               // Start a new cluster
	  new_cluster.push_back(hits[t]);    // new cluster now contains the current hit
	  previous_end_position = hit_positions[i][hits[t]] + k - 1;
	}
      }
      if (new_cluster.size() == 1)
	clusters.push_back(new_cluster);  // last cluster


	  
      ////////////////////////
      // Iterate the clusters

      std::vector<int> cluster;         // cluster contains indices to hit_patterns vector
      BOOST_FOREACH(cluster, clusters) {
	non_intersecting_occurrences.push_back(cluster[0]);
      }  // end BOOST_FOREACH cluster

      */
 
    }

    for (int t=0; t < non_intersecting_occurrences.size(); ++t) {
      int index = non_intersecting_occurrences[t];
      boost::tie(pattern, pairs) = *(hit_contributions[i][index]);
      /*
      int start_pos = hit_positions[i][index];
      int dir = hit_directions[i][index];
      std::string pal = is_palindromic(pattern) ? "Palindrome" : "-";
      int hd = hamming_distance(seed, pattern);
      printf("*%i\t%i\t%i\t%i\t%s\t%s\t%i\n", i, start_pos, dir, k, pattern.c_str(), pal.c_str(), hd);
      */
      if (iupac_string_match(pattern, seed))
	++seed_count;
      
      ++total_count;
      int j, a;
      BOOST_FOREACH(boost::tie(j,a), pairs) {
	result(a, j) += 1;
      }
    }
  }  // end i
      
  return boost::make_tuple(result, seed_count, total_count);
};


boost::tuple<dmatrix,int>
find_multinomial_n_suffix_array(const std::string& seed, const std::vector<std::string>& sequences, const suffix_array& sa, int n, bool use_multimer)
{
  const int k = seed.length();
  //  const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  //char nucs[] = "ACGT";
  dmatrix result(4, k);
  //code_to_tuple_type code_to_tuple;

  string_to_tuple_type string_to_contributions;
  string_to_contributions = get_n_neighbourhood_mononucleotide_contributions(seed, n);

  printf("Number of patterns %lu\n", string_to_contributions.size());
  unsigned long seed_count=0;
  unsigned long total_count=0;
  //unsigned long number_of_sites=0;

  std::string seed_rev = reverse_complement(seed);
  switch (data_counting) {
  case all_occurrences:
    boost::tie(result, seed_count, total_count) =
      count_all_occurrences(string_to_contributions, seed, sa, sequences);
    break;
  case choose_one_per_cluster:  
    boost::tie(result, seed_count, total_count) = count_choose_one_per_cluster(string_to_contributions, seed, sa, sequences);
    break;
  case neighbourhood_contains_one:
    boost::tie(result, seed_count, total_count) = count_neighbourhood_contains_one(string_to_contributions, seed, sa, sequences);
    break;
  case sequence_contains_at_least_one:
    error(true, "Not implemented.");
    break;
  case sequence_contains_one:
    boost::tie(result, seed_count, total_count) = count_sequence_contains_one(string_to_contributions, seed, sa, sequences);
    break;
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
