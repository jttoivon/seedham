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
#include "kmer_tools.hpp"
#include "multinomial_helper.hpp"
#include "seed_basic.hpp"

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

//bool use_background_correction=false;
bool use_dinucleotide_model=false;
bool use_cell_correction=false;
bool correct_for_seed_bias = false;
bool contains_N = false;
bool lambda_given = false;   // whether lambda was given as command line parameter

std::string real_matrix_filename; // for testing purposes

prior<double> pseudo_counts;


std::vector<double> background_frequencies(4);
std::vector<double> background_probabilities(4);
matrix<double> background_frequency_matrix(4,4);   // for background noise
matrix<double> background_probability_matrix(4,4); // for background noise

// this accepts multiple occurrences per sequence
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



// print occurrence statistics for all possible k-mers
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
	if (d==l) {          // found occurrence
	  ++all[dna_to_number<size_t>(line.substr(j,k))];
	  seqs[dna_to_number<size_t>(line.substr(j,k))].insert(i);
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
    printf("On average, %lg occurrences per sequence\n", ratio/count);
    printf("%.3lg%% of the above strings occur only once per sequence\n",
	   (double)count2/count*100);
  }
}

// Find out number of sites in data within Hamming distance 'delta' from 'seed'
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

// print occurrence statistics for all Hamming neighbourhoods H(str, d) for 0<=d<=delta
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
    printf("\n%i occurrences of the Hamming neighbourhood H(str,%i) on %i sequences\n", 
	   all, l, seqs);
    int count2 = std::count(occs.begin(), occs.end(), 1);
//     for (int i=0; i < lines; ++i)
//       if (occs[i]==1)
// 	++count2;
    printf("On average, %lg occurrences per sequence\n", (double)all/seqs);
    printf("%.3lg%% of sequences have exactly one occurrence\n",
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



// Number of subsequences used to build the pwm
// NOTE! This only works when the Hamming radius is 1x
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







int 
myskip(int i, int skip)
{
  assert(i >= 0);
  assert(i <= 2);
  assert(skip >= 0);
  assert(skip < 4);

  return i >= skip ? i+1 : i;
}


// Convert memory area of length 'bytes' bytes to a string of chars '0' and '1'
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



// This computes the init_motif matrix counts.
// Chooses the first occurrence with lowest Hamming distance for each sequence
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




// Returns a 0-1-matrix corresponding to 'seed', which can be an IUPAC sequence
dmatrix
make_consensus_matrix(const std::string& seed, int count)
{
  int k = seed.length();
  dmatrix result(4, k);
  for (int i=0;i<k;++i)
    result.set_column(i, iupac_probability(seed[i]));
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

typedef boost::tuple<dmatrix,int> (*func_ptr_t)(const std::string&, const std::vector<std::string>&, int);
typedef dmatrix (*func_ptr2_t)(const std::string&, const std::vector<std::string>&, int);


// For testing purposes
char GATA[] =
  "0.468455	0.000000	1.000000	0.000000	0.776771	0.616743	0.144221	0.371156\n"
  "0.215102	0.000000	0.000000	0.000000	0.054366	0.115108	0.301166	0.209968\n"
  "0.071038	1.000000	0.000000	0.000000	0.038715	0.165468	0.452810	0.309650\n"
  "0.245405	0.000000	0.000000	1.000000	0.130148	0.102681	0.101803	0.109226\n";

// Creates a matrix of width k+2*g with another matrix 'm' in the middle.
// Rest of the columns are filled with q
dmatrix
make_helper_pfm(const dmatrix& m, const dvector& q, int g)
{
  int k = m.get_columns();
  int len = k + 2*g;
  dmatrix result(4, len);
  for (int i=0; i < g; ++i) {
    result.set_column(i, q);
    result.set_column(i+k+g, q);
  }
  result.inject(m, 0, g);
  
  return result;
}

dmatrix
find_model(func_ptr_t func_ptr, const std::string& name, const std::string& seed, int hamming_radius,
	   std::vector<std::string>& sequences, std::vector<std::string>& sequences_bg, const std::vector<double>& data_bg,
	   double lambda)
{
  TIME_START(t);
  int k = seed.length();
  int lines = sequences.size();
  int lines_bg = sequences_bg.size();
  int L = sequences[0].length();
  //  int N = lines*(L-k+1); // number of sites
  dmatrix motif;
  int total_count;
  TIME_START(t2);
  boost::tie(motif, total_count) = func_ptr(seed, sequences, hamming_radius);
  TIME_PRINT("Multinomial-n observed took %.2f seconds.\n", t2);

  TIME_START(t3);
  if (background_correction == file_correction) {   // This bg correction uses previous generation SELEX data
    printf("\n");
    write_matrix(stdout, motif, to_string("Signal %s counts before correction:\n", name.c_str()), "%.0f");

    dmatrix motif_bg;
    int bg_total_count;
    boost::tie(motif_bg, bg_total_count) = func_ptr(seed, sequences_bg, hamming_radius);
    write_matrix(stdout, motif_bg, to_string("Background %s counts before correction:\n", name.c_str()), "%.0f");
    printf("Lambda is %f\n", lambda);
    dmatrix m_new;
    if (use_cell_correction)
      m_new = (double)lines_bg/lines*motif - (1.0-lambda)*motif_bg;
    else
      m_new = motif - ((double)lines/lines_bg*(1.0-lambda))*motif_bg;

    m_new.apply(cut);
    motif = m_new;
    printf("\n");
  } else if ((background_correction == uniform_correction or background_correction == data_correction) and k <= 10) {   // NOTE! FIX THIS! MAGIC NUMBER AS WELL!

    // This bg correction do NOT use previous generation SELEX data
    
    //assert(use_multimer);
    write_matrix(stdout, motif, to_string("Signal %s counts before correction:\n", name.c_str()), "%.0f");
    dmatrix expected_matrix;
    double bg_total_prob;
    std::vector<double> uniform_bg(4, 0.25);
    const std::vector<double>& bg = background_correction == uniform_correction ? uniform_bg : data_bg;
    boost::tie(expected_matrix, bg_total_prob) =
      find_multinomial_n_background(seed, sequences, bg,
				    hamming_radius, use_multimer);
    dmatrix m_new;
    int sites = 0;
    if (data_counting == all_occurrences) 
      sites = lines * (L-k+1);     // For all occurrences method
    else if (data_counting == neighbourhood_contains_one) {
      //  int sites = lines * (L - (k+2*g) + 1);
      sites = lines * (L - k + 1);
    }
    printf("Total %s bg count is %f\n", name.c_str(), sites*bg_total_prob);
    dmatrix motif_bg;
    if (data_counting == all_occurrences) {
      if (not lambda_given)
	lambda = 0.0;
      printf("Lambda is %f\n", lambda);
      motif_bg = (1.0-lambda)*sites*expected_matrix;
    }
    else if (data_counting == neighbourhood_contains_one) {   // assumed to be neighbour method
      int g = cluster_threshold;
      double bg_total_prob2;
      //      dmatrix estimate = motif-motif_bg;
      //estimate.apply(cut);
      dmatrix estimate;
      dmatrix mstar;
      if (real_matrix_filename != "") {
	estimate = read_matrix_file(real_matrix_filename);
	mstar = make_helper_pfm(0.00001+estimate, bg, g);  // Note! Pseudo count added
      }
      else {
	estimate = motif;
	mstar = make_helper_pfm(1+estimate, bg, g);  // Note! Pseudo count added
      }
      normalize_matrix_columns(mstar);
      write_matrix(stdout, mstar, to_string("mstar:\n"), "%.2f");
      typedef boost::multi_array_types::extent_range range;
      boost::multi_array<dmatrix,1> shifted_bg_pfms(boost::extents[range(-(k+g-1), k+g)]);
      for (int j = -(k+g-1); j <= k+g-1; ++j) {
	if (j == 0)
	  continue;
	dmatrix mjstar = get_shifted_window_pwm(mstar, bg, j);
	dmatrix em;
	boost::tie(em, bg_total_prob2) =
	  neighbour_expected_pfm_in_pfm_distribution(seed, mjstar, hamming_radius);
	shifted_bg_pfms[j] = em;
      }
      
      if (not lambda_given)
	lambda = estimate_signal_fraction(motif, expected_matrix, shifted_bg_pfms, sequences, bg);
      printf("Lambda is %f\n", lambda);
      //      motif_bg = (1.0 - lambda*(2*g+1)) * sites * expected_matrix;
      motif_bg = (1.0 - lambda*(2*(g+k)-1)) * sites * expected_matrix;
      write_matrix(stdout, motif_bg, to_string("%s background matrix:\n", name.c_str()), "%.5f");
      for (int j = -(k+g-1); j <= k+g-1; ++j) {
	if (j == 0)
	  continue;
	dmatrix temp = lambda * sites * shifted_bg_pfms[j];
	write_matrix(stdout, temp, to_string("contribution j=%i:\n", j), "%f");
	motif_bg += temp;
	//printf("Total j=%i %s bg count is %f\n", j, name.c_str(), sites*bg_total_prob2);
      }
    }
    else {
      error(true, "Not implemented!");
    }
    
    write_matrix(stdout, motif_bg, to_string("%s background matrix:\n", name.c_str()), "%.0f");
    m_new = motif - motif_bg;
    m_new.apply(cut);
    motif = m_new;
    printf("\n");

  }
  TIME_PRINT("Multinomial-n background took %.2f seconds.\n", t3);

  write_matrix(stdout, motif, to_string("%s motif matrix counts:\n", name.c_str()), "%.0f");
  dmatrix norm=motif;
  normalize_matrix_columns(norm);
  
  write_matrix(stdout, norm, to_string("%s motif matrix:\n", name.c_str()), "%.6f"); 
  print_ics(norm);
  printf("\n");
  TIME_PRINT("Total multinomial-n took %.2f seconds.\n", t);
  return motif;
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
  use_multimer=false;   // allow multiple occurrences per sequence
  double lambda=0.0;
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  // 
  // Do the initialization
  //
  ///////////////////////////////////////////////////////////////////////////////////////////

  // Declare the supported options.
  po::options_description nonhidden("Allowed options");
  nonhidden.add_options()
    ("help", "produce help message")
    ("hamming-radius", po::value<int>(), m("Maximum Hamming radius", hamming_radius).c_str())
    //    ("stairs", m("Use weighted Hamming distance", use_weighted_hamming).c_str())
    ("statistics", "Print matching statistics")
    ("single-strand", m("Assume sequences can only come from forward strand", not use_two_strands).c_str())
    ("reverse-strand", m("Assume sequences come from reverse strand", false).c_str())
    ("count-palindromes-twice", m("Count palindromes twice", count_palindromes_twice).c_str())
    ("count", po::value<std::string>(), "Counting method: all, neighbour, or cluster, default: all")
    ("seed-count", po::value<std::string>(), "Counting method for finding seeds: all, unique, once, default: all")
    ("cluster-threshold", po::value<int>(), m("Cluster threshold", cluster_threshold).c_str())
    ("pseudo-counts", "Use pseudo counts, default: no")
    //    ("correct-for-seed-bias", m("This should be used for learning dinucleotide model.",
    //			      correct_for_seed_bias).c_str())
    //    ("use-bernoulli-background", m("Use bernoulli model to subtract background", use_bernoulli_background).c_str())
    // ("prior", po::value<std::string>(), 
    //  "Choose either addone or dirichlet prior")
    ("background", po::value<std::string>(),      
     "Either 'uniform', 'data', or the filename of the background sequences, enables also background correction. "
     "The options 'uniform' and 'data' use the Bernoulli model for background.")
    ("lambda", po::value<double>(), m("Proportion of signal in data. If not given and background correction is in use, the lambda will be estimated from the data.", lambda).c_str())
    ("output-matrix", po::value<std::string>(), "Name of the output matrixfile, default: none")
    //    ("real-matrix", po::value<std::string>(), "Name of the real matrixfile, default: none")
    //    ("multimer", m("Use multimer version for multinomial-1/2", use_multimer).c_str())
    //    ("palindromic-correction", m("Correct mixing of direction by Esko's method", use_palindromic_correction).c_str())
    ("palindromic-index",  po::value<std::string>(), m("Require that seed has high-enough palindromic index: automatic (depending on Hamming radius), 0, 1, ...", palindromic_index_limit).c_str())
    //    ("dependence-matrix", m("Find dependence matrix of result", find_dependence_matrix).c_str())
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

    if (vm.count("stairs")) {
      use_weighted_hamming=true;
    }   

    if (vm.count("single-strand"))
      use_two_strands = false;

    if (vm.count("count-palindromes-twice"))
      count_palindromes_twice = true;

    if (vm.count("reverse-strand")) {
      use_two_strands = false;
      use_reverse_strand = true;
    }


    if (vm.count("dinucleotide-model")) {
      use_dinucleotide_model = true;
      hamming_radius = 2;
    }
    if (vm.count("pseudo-counts"))
      use_pseudo_counts = true;

    // Note! This has to be below the handling of --dinucleotide-model
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

    
    if (vm.count("cluster-threshold")) {
      cluster_threshold = vm["cluster-threshold"].as< int >();
    }
    
    if (vm.count("count")) {
      std::string param = vm["count"].as<std::string>();
      if (param == "all") {
	data_counting = background_counting = all_occurrences;
      }
      else if (param == "cluster") {
	data_counting = background_counting = choose_one_per_cluster;
      }
      else if (param == "neighbour") {
	data_counting = background_counting = neighbourhood_contains_one;
      }
	
      use_multimer = true;
    }

    if (vm.count("seed-count")) {
      std::string param = vm["seed-count"].as<std::string>();
      if (param == "all") {
	seed_counting = all_occurrences;
      }
      else if (param == "unique") {
	seed_counting = sequence_contains_one;
      }
      else if (param == "once") {
	seed_counting = sequence_contains_at_least_one;
      }
	
    }

    if (vm.count("background")) {
      std::string t = vm["background"].as< string >();
      if (t == "uniform") {
	background_correction = uniform_correction;
      }
      else if (t == "data") {
	background_correction = data_correction;
      }
      else {
	background_correction = file_correction;
	background_file = vm["background"].as< string >();

      }
    }
    
    if (vm.count("use-cell-correction")) {
      error(vm.count("background") == 0, "The '--use-cell-correction' option must be used with the '--background filename option'.");
      use_cell_correction = true;
    }
    
    if (vm.count("dependence-matrix"))
      find_dependence_matrix = true;

    if (vm.count("print-neighbourhood"))
      print_neighbourhood = true;

    if (vm.count("print-alignment"))
      print_alignment = true;

    if (vm.count("correct-for-seed-bias"))
      correct_for_seed_bias = true;


    if (background_correction != file_correction and vm.count("lambda")) {
      lambda = vm["lambda"].as<double>();
      lambda_given = true;
    }
    
    // if (vm.count("prior")) {
    //   prior_parameter = vm["prior"].as< string >();
    //   if (prior_parameter == "addone" || prior_parameter == "dirichlet")
    // 	use_pseudo_counts=true;
    //   else {
    // 	fprintf(stderr, "Invalid prior: %s\n", prior_parameter.c_str());
    // 	exit(1);
    //   }
    // }   
    
    

    if (vm.count("output-matrix")) 
      matrixfile = vm["output-matrix"].as< string >();
    
    if (vm.count("real-matrix")) 
      real_matrix_filename = vm["real-matrix"].as< string >();
    
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


  //  boost::tie(lines, bad_lines) = read_sequences(seqsfile, sequences, true);
  boost::tie(lines, bad_lines) = read_sequences(seqsfile, sequences, false);   // Don't allow IUPAC
  



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

  printf("Background correction: %s\n", background_correction_type_to_string[background_correction]);
  printf("Counting type for data is %s\n", counting_type_to_string[data_counting]);
  if (background_correction != no_correction and background_correction != file_correction) {
    printf("Counting type for background is %s\n", counting_type_to_string[background_counting]);
  }
  
  printf("Cluster threshold is %i\n", cluster_threshold);
  printf("Hamming radius is %i\n", hamming_radius);
  printf("Use weighted Hamming distance: %s\n", yesno(use_weighted_hamming));
  printf("Use two dna strands: %s\n", use_two_strands ? "yes" : "no");
  printf("Count palindromes twice: %s\n", count_palindromes_twice ? "yes" : "no");
  printf("Use reverse strand: %s\n", use_reverse_strand ? "yes" : "no");
  //  printf("Use bernoulli background: %s\n", use_bernoulli_background ? "yes" : "no");
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
    boost::tie(seed,seed_count) = most_common_pattern(sequences, k, seed_param, contains_N, hamming_radius);
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

  if (false) {
    if (not use_multimer and L < 60)
      distribution_of_startpositions(sequences, seed);
  }
  
  int lines_bg=0, bad_lines_bg;
  std::vector<std::string> background_seqs;
  if (background_correction == file_correction) {
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
    lambda = 1.0 - lambda; // elsewhere signal lambda is used, hence the conversion

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
  func_ptr_t func_ptr;// = use_multimer ? find_snips_multimer : find_snips_monomer;
  dmatrix multinomial1_motif;
  /*

  if (k <= 64) 
    multinomial1_motif = find_model(func_ptr, "multinomial-1", seed, 1, sequences, background_seqs, background_probabilities,
					    lambda);
  */


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

  dmatrix multinomial_n_motif;
  if (not use_dinucleotide_model) {
    if (k >= 20 or hamming_radius >= 6) 
      func_ptr = find_multinomial_n_scan;
    else
      func_ptr = find_multinomial_n;
    if (hamming_radius < 1) {
      printf("Warning! Hamming radius of at least 1 should be used to get unbiased pfm model\n");
    }
  
    multinomial_n_motif= find_model(func_ptr, "multinomial-n", seed, hamming_radius, sequences, background_seqs, background_probabilities,
				    lambda);
  }
  

  /////////////////////////////////////
  //
  // compute the positional background
  //
  /////////////////////////////////////

  if (false) {
    if (not use_multimer and L < 60) {
      positional_background = count_positional_background(sequences);
    
      normalize_matrix_columns(positional_background);
      assert(is_column_stochastic_matrix(positional_background));
      write_matrix(stdout, positional_background, "Positional background probility matrix:\n", "%.6f");
    }
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




  if (matrixfile != "-") {
    write_matrix_file(matrixfile, multinomial_n_motif, "%.0f");
  }

  if (statistics) {
    print_hamming_statistics(sequences, seed, hamming_radius);
    print_hamming_statistics2(sequences, seed, hamming_radius);
  }

  return 0;
}


