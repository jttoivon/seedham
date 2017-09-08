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
#ifndef MULTINOMIAL_HELPER_HPP
#define MULTINOMIAL_HELPER_HPP

#include "matrix.hpp"
#include "suffix_array_wrapper.hpp"
#include "data.hpp"
#include "aho_corasick_wrapper.hpp"

#include <boost/tuple/tuple.hpp>
#include <string>
#include <vector>
#include <boost/unordered_map.hpp>
#include <boost/multi_array.hpp>


extern bool use_palindromic_correction;
extern int extra_flank;
extern int cluster_threshold;
extern bool use_weighted_hamming;
extern bool use_multimer;

enum counting_type {
  all_occurrences,                      // no restrictions
  sequence_contains_one,                // count only those sequences that contain exactly one occurrence
  sequence_contains_at_least_one,       // 
  neighbourhood_contains_one,           // count only those occurrences that have no neighbours (a neighbour is an intersecting site)
  choose_one_per_cluster     // count only one occurrence from a cluster
};

enum background_correction_type {
  no_correction,                                 // no correction
  uniform_correction,                              // Bernoulli model using uniform distribution
  data_correction,                                 // Bernoulli model using distribution learned from the sequnce data file
  file_correction                                  // Background model is the set of sequences in a given file
};

static std::map< counting_type, const char * > counting_type_to_string = {
   {all_occurrences, "all_occurrences"},
   {sequence_contains_one, "sequence_contains_one"},
   {neighbourhood_contains_one, "neighbourhood_contains_one"},
   {choose_one_per_cluster, "choose_one_per_cluster"}
};

static std::map< background_correction_type, const char * > background_correction_type_to_string = {
  {no_correction,      "no_correction"},    
  {uniform_correction,   "uniform_correction"}, 
  {data_correction,      "data_correction"},    
  {file_correction,       "file_correction"}     
};


extern background_correction_type background_correction;
extern counting_type data_counting;
extern counting_type background_counting;


typedef boost::unordered_map<int, int> id_to_count_type;
//typedef std::vector<int> id_to_count_type;
typedef boost::unordered_map<big_int, std::vector<boost::tuple<int, int> > > code_to_tuple_type;
// This actually maps a string to vector of (position,nucleotide) pairs.
typedef boost::unordered_map<std::string, std::vector<boost::tuple<int, int> > > string_to_tuple_type;



/*
boost::tuple<std::string,int>
most_common_pattern_monomer(const std::vector<std::string>& sequences, int k, std::string seed = "",
			    int hamming_radius=0);
*/

// Do not reject sequences with multiple occurrences of query strings
// compute the counts for the multinomial1 matrix
dmatrix
find_snips_multimer(const std::string& seed, const std::vector<std::string>& sequences, int hamming_distance);


boost::tuple<dmatrix,int,int>
find_snips_multimer_helper(const std::string& seed, const std::vector<std::string>& sequences);


string_to_tuple_type
get_n_neighbourhood_mononucleotide_contributions(const std::string&seed, int n);

std::vector<std::pair<std::string, std::vector<boost::tuple<int, int> > > >
get_n_neighbourhood_in_vector(const std::string&seed, int n);

boost::tuple<dmatrix,int>
find_multinomial_n_suffix_array(const std::string& seed, const std::vector<std::string>& sequences, const suffix_array& sa, int n, bool use_multimer);

boost::tuple<dmatrix,double>
neighbour_expected_pfm_in_pfm_distribution(const std::string& seed,
					   const dmatrix& m,
					   int hamming_radius);

boost::tuple<dmatrix,double>
find_multinomial_n_background(const std::string& seed, const std::vector<std::string>& sequences, const std::vector<double>& bg,
			      int n, bool use_multimer);

dmatrix
align_all(const std::vector<std::string>& sequences);

double
estimate_real_count(int observed_count, const std::string& seed, const dmatrix& m, int hamming_radius, int L,
		    const string_to_tuple_type& string_to_contributions, const std::vector<double>& bg);


double
estimate_real_count_version_2(int observed_count, const std::string& seed, const dmatrix& m, int hamming_radius, int L,
			      const string_to_tuple_type& string_to_contributions, const std::vector<double>& bg, double gamma=0.2);

dmatrix
get_shifted_window_pwm(const dmatrix& m, const dvector& q, int j);

void
print_ics(const dmatrix& m);

// this computes the multinomial-n matrix counts using suffix array (was AC-automaton)
boost::tuple<dmatrix,int>
find_multinomial_n(const std::string& seed, const std::vector<std::string>& sequences, int n);

// this computes the multinomial-n matrix counts by scanning all possible windows
boost::tuple<dmatrix,int>
find_multinomial_n_scan(const std::string& seed, const std::vector<std::string>& sequences, int n);

int match_handler_counts(MATCH* m, void* param);

double 
estimate_signal_fraction(const dmatrix& pfm, dmatrix bg_pfm, const boost::multi_array<dmatrix,1>& shifted_bg_pfms,
			 const std::vector<std::string>& sequences, const std::vector<double>& q);

#endif // MULTINOMIAL_HELPER_HPP
