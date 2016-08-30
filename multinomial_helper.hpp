#include "matrix.hpp"
#include "suffix_array_wrapper.hpp"

#include <boost/tuple/tuple.hpp>
#include <string>
#include <vector>
#include <boost/unordered_map.hpp>


extern int palindromic_index_limit;
extern bool use_palindromic_correction;


typedef boost::unordered_map<std::string, std::vector<boost::tuple<int, int> > > string_to_tuple_type;

boost::tuple<std::string,int>
most_common_pattern_multimer(const std::vector<std::string>& sequences, int k, std::string seed,
			     bool contains_N, int hamming_radius=0);

int
conflict_free_palindromic_index(int hamming_radius);


boost::tuple<std::string,int>
most_common_pattern_monomer(const std::vector<std::string>& sequences, int k, std::string seed = "",
			    int hamming_radius=0);


// Do not reject sequences with multiple occurrences of query strings
// compute the counts for the multinomial1 matrix
dmatrix
find_snips_multimer(const std::string& seed, const std::vector<std::string>& sequences, int hamming_distance);


boost::tuple<dmatrix,int,int>
find_snips_multimer_helper(const std::string& seed, const std::vector<std::string>& sequences);


string_to_tuple_type
get_n_neighbourhood(const std::string&seed, int n);

std::vector<std::pair<std::string, std::vector<boost::tuple<int, int> > > >
get_n_neighbourhood_in_vector(const std::string&seed, int n);

boost::tuple<dmatrix,int>
find_multinomial_n_suffix_array(const std::string& seed, const std::vector<std::string>& sequences, const suffix_array& sa, int n, bool use_multimer);

boost::tuple<dmatrix,int>
find_multinomial_n_background(const std::string& seed, const std::vector<std::string>& sequences, const std::vector<double>& bg,
			      int n, bool use_multimer);

dmatrix
align_all(const std::vector<std::string>& sequences);