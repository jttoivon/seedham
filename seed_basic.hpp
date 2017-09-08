#include "multinomial_helper.hpp"

#include <string>
#include <vector>

extern int palindromic_index_limit;
extern counting_type seed_counting;

int
conflict_free_palindromic_index(int hamming_radius);


boost::tuple<std::string,int>
most_common_pattern(const std::vector<std::string>& sequences, int k, std::string seed,
			     bool contains_N, int hamming_radius=0);
