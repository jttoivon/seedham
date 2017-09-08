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
#include "seed_basic.hpp"
#include "kmer_tools.hpp"
#include "parameters.hpp"
#include "common.hpp"
#include "huddinge.hpp"

#include <boost/tuple/tuple.hpp>
#include <boost/foreach.hpp>

int fixed_low_count_limit = 20;
int palindromic_index_limit=0;
counting_type seed_counting = all_occurrences;

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

// The triples are of (code, count, palindromic index) type 
// Primary sort key is the count, the secondary sort key is the palindromic index
bool
triple_comp(boost::tuple<big_int, int, int> a, boost::tuple<big_int, int, int> b)
{
  return a.get<1>() > b.get<1>() or (a.get<1>() == b.get<1>() and a.get<2>() > b.get<2>());
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


// A sequence must have at most occurrence per line to be considered
count_container_t
get_counts_of_unique(const std::vector<std::string>& sequences, int k)
{
  count_container_t number_of_occurrences;// items are (code,count) pairs
  bool count_palindromes_twice = use_two_strands;  // There is also global version of this

  int lines = sequences.size();
  //  int max_count=-1;
  //  big_int argmax=-1;
  big_int id;
  big_int id2;
  for  (int i=0; i < lines; ++i) {
    const std::string& line = sequences[i];
    //std::map<big_int, std::vector<int> >occurences;

    boost::unordered_map<big_int, std::vector<int> > occurrences;  // occurrences on this line
    std::set<big_int> ids;

    // find all occurrences in sequence
    for (int j=0; j < line.length()-k+1; ++j) {
      id = dna_to_number<big_int>(line.substr(j,k));
      occurrences[id].push_back(j);
      ids.insert(id);
      if (use_two_strands) {
	id2 = dna_to_number<big_int>(reverse_complement(line.substr(j,k)));
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
	/*	
	if (number_of_occurrences[id] > max_count) {
	  max_count = number_of_occurrences[id];
	  argmax = id;
	}
	if (use_two_strands && number_of_occurrences[id2] > max_count) {
	  max_count = number_of_occurrences[id2];
	  argmax = id2;
	}
	*/	  
      }

    }
  }

  return number_of_occurrences;
}

// A sequence gets count 1 if it appears at least once per line otherwise the count for that line is zero
count_container_t
get_counts_once_per_line(const std::vector<std::string>& sequences, int k)
{
  count_container_t number_of_occurrences;// items are (code,count) pairs
  //  bool count_palindromes_twice = use_two_strands;  // There is also global version of this

  int lines = sequences.size();
  //  int max_count=-1;
  //  big_int argmax=-1;
  big_int id;
  big_int id2;
  for  (int i=0; i < lines; ++i) {
    const std::string& line = sequences[i];

    boost::unordered_map<big_int, std::vector<int> > occurrences;  // occurrences on this line
    std::set<big_int> ids;

    // find all occurrences in sequence
    for (int j=0; j < line.length()-k+1; ++j) {
      id = dna_to_number<big_int>(line.substr(j,k));
      occurrences[id].push_back(j);
      ids.insert(id);
      if (use_two_strands) {
	id2 = dna_to_number<big_int>(reverse_complement(line.substr(j,k)));
	occurrences[id2].push_back(j);
      }
    }
    // accept only subsequences that appear only once per sequence, or is a palindromic occurrence
    for (std::set<big_int>::iterator it=ids.begin(); it != ids.end(); ++it) {
      big_int id = *it;
      std::vector<int>& r = occurrences[id];
      if (r.size() > 0) {
	big_int id2 = reverse_complement_2bitstring(id, k);
	++number_of_occurrences[id];
	if (use_two_strands)
	  ++number_of_occurrences[id2];
	
      }

    }
  }

  return number_of_occurrences;
}

boost::tuple<std::string,int>
most_common_pattern(const std::vector<std::string>& sequences, int k, std::string seed,
			     bool contains_N, int hamming_radius)
{
  assert(k > 1 && k <= max_matrix_len); 
  int lines = sequences.size();
  int L = sequences[0].length();
  
  //double low_count_limit = fixed_low_count_limit;
  double p = pow(4, -k);
  int sites = lines * (L-k+1) * (use_two_strands ? 2 : 1);
  
  double stddev = sqrt(sites*p*(1-p));
  double expected = sites*p;
  double low_count_limit = std::max(expected + 2*stddev, (double)fixed_low_count_limit);
  
  printf("PI optimization low count limit is %f\n", low_count_limit);
  
  count_container_t number_of_occurrences;// items are (code,count) pairs
  bool count_palindromes_twice = use_two_strands;  // There is also global version of this

  //  counting_type seed_counting = sequence_contains_at_least_one;
  if (seed_counting==all_occurrences) {
    if (contains_N)   // Data contains 'N's
      get_kmer_counts(remove_masked_areas(sequences, k), k, number_of_occurrences, use_two_strands, count_palindromes_twice);
    else
      get_kmer_counts(sequences, k, number_of_occurrences, use_two_strands, count_palindromes_twice);
  }
  else if (seed_counting == sequence_contains_one) 
    number_of_occurrences = get_counts_of_unique(sequences, k);
  else if (seed_counting == sequence_contains_at_least_one)
    number_of_occurrences = get_counts_once_per_line(sequences, k);
  else {
    error(true, "Option for seed_counting not supported");
  }
  
  code_t argmax;
  std::vector<boost::tuple<big_int, int, int> > v;
  code_t code;
  int count;
  BOOST_FOREACH(boost::tie(code, count), number_of_occurrences)
    v.push_back(boost::make_tuple(code, count, palindromic_index(number_to_dna(code, k))));
  std::sort(v.begin(), v.end(), triple_comp);   // Primary sort key is the count, the secondary sort key is the palindromic index

  // These are the count and pi of the most common kmer
  int top_count;
  int top_pi;
  code_t top_code;
  boost::tie(top_code, top_count, top_pi) = v[0];
  printf("!First candidate for seed: %s %i %i\n", number_to_dna(top_code, k).c_str(), top_count, top_pi);
  double ratio_cutoff = 0.25;   // Increase of two units in PI can be allowed to drop the corresponding count into one quarter
  if (palindromic_index_limit > 0) {
    // the palindromic index of a seed needs to be at least 'limit' for the n-Hamming-neighbourhood to be conflict-free
    code_t code=0;
    code_t max_code = top_code;
    int max_pi = top_pi;
    int max_count = top_count;
    std::vector<int> alignments;
    std::string topseed =  number_to_dna(v[0].get<0>(), k);   // Seed with highest count
    std::string topseed_revcomp = reverse_complement(topseed);
    for (int i=0; i < v.size(); ++i) {
      code = v[i].get<0>();
      int count = v[i].get<1>();
      int pi = v[i].get<2>();
      //      if (count < low_count_limit and max_pi >= 0)   // If we have already one candidate and the current kmer is too small, then quit.
      if (count < low_count_limit)   // If the current kmer is too small, then quit.
	break;
      std::string temp = number_to_dna(code, k);
      if (huddinge_distance(topseed, temp) <= huddinge_distance(topseed_revcomp, temp))
	alignments = huddinge_alignment(topseed, temp);
      else
	alignments = huddinge_alignment(topseed_revcomp, temp);
      
      // middle part of the condition checks that candidate is not a shift of topseed (or its reverse complement)
      //      if (pi > max_pi and alignments.size() == 1 and alignments[0] == 0 and (float)count/top_count >= pow(ratio_cutoff, (pi - top_pi)/2)) {   

      if (pi > max_pi and alignments.size() == 1 and alignments[0] == 0 and (float)count/max_count >= ratio_cutoff) {   
	max_pi = pi;
	max_code = code;
	max_count = count;
	printf("!New candidate for seed: %s %i %i\n", temp.c_str(), count, pi);
      }
      if (max_pi >= palindromic_index_limit)
	break;                   // found a good seed
    }
    argmax = max_code;
  }
  else {
    argmax = v[0].get<0>();
  }
  
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

  printf("Seed %s has %i occurrences\n",
	 result.c_str(),number_of_occurrences[dna_to_number<code_t>(result)]);

  return boost::make_tuple(result, number_of_occurrences[dna_to_number<code_t>(result)]);
  //return result;
}


/*
boost::tuple<std::string,int>
most_common_pattern_monomer(const std::vector<std::string>& sequences, int k, std::string seed,
			    bool contains_N, int hamming_radius)
{
  assert(k > 1 && k <= max_matrix_len); 
  count_container_t number_of_occurrences;// items are (code,count) pairs
  bool count_palindromes_twice = use_two_strands;  // There is also global version of this

  int lines = sequences.size();
  int max_count=-1;
  big_int argmax=-1;
  big_int id;
  big_int id2;
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
  code_t code;
  int count;
  BOOST_FOREACH(boost::tie(code, count), number_of_occurrences)
    v.push_back(boost::make_tuple(code, count, palindromic_index(number_to_dna(code, k))));
  std::sort(v.begin(), v.end(), triple_comp);   // compares according the second member of the pair: the count
  // These are the count and pi of the most common kmer
  int top_count;
  int top_pi;
  code_t top_code;
  boost::tie(top_code, top_count, top_pi) = v[0];

  double ratio_cutoff = 0.25;   // Increase of two units in PI can be allowed to drop the corresponding count into one quarter
  if (palindromic_index_limit > 0) {
    code_t code=0;
    code_t max_code = top_code;
    int max_pi = top_pi;
    int max_count = top_count;
    std::vector<int> alignments;
    std::string topseed =  number_to_dna(v[0].get<0>(), k);   // Seed with highest count
    std::string topseed_revcomp = reverse_complement(topseed);
    for (int i=0; i < v.size(); ++i) {
      code = v[i].get<0>();
      int count = v[i].get<1>();
      int pi = v[i].get<2>();
      if (count < low_count_limit)   // If we have already one candidate and the current kmer is too small, then quit.
	break;

      std::string temp = number_to_dna(code, k);
      if (huddinge_distance(topseed, temp) <= huddinge_distance(topseed_revcomp, temp))
	alignments = huddinge_alignment(topseed, temp);
      else
	alignments = huddinge_alignment(topseed_revcomp, temp);
      
      // latter part of the condition checks that candidate is not a shift of topseed (or its reverse complement)
      if (pi > max_pi and alignments.size() == 1 and alignments[0] == 0 and (float)count/max_count >= ratio_cutoff) {
	max_pi = pi;
	max_code = code;
	max_count = count;
	printf("!New candidate for seed: %s %i %i\n", temp.c_str(), count, pi);
      }
      if (max_pi >= palindromic_index_limit)
	break;                   // found a good seed
    }
    argmax = max_code;
  } else
    argmax = v[0].get<0>();

  // between string and its reverse complement, choose lexicographically smaller
  std::string result;
  if (seed.length() == k)
    result = seed;
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
}

*/

