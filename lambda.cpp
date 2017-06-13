#include "lambda.hpp"
#include "common.hpp"
#include "data.hpp"
#include "parameters.hpp"

#include <boost/tuple/tuple.hpp>

//#include <vector>
#include <set>
#include <cmath>



int myget1(boost::tuple<int, int> t)
{
  return boost::get<0>(t);
}

int myget2(boost::tuple<int, int> t)
{
  return boost::get<1>(t);
}

int
find_nonspecific_fraction(const std::vector<std::string>& sequences)
{
  int lines = sequences.size();
  int k=8;
  int L=sequences[0].length();
  int m=L-k+1;


  // kmers is a set of pairs (kmer_id,count)
  std::vector<boost::tuple<int, int> > kmers(pow(4,k), boost::make_tuple(0,0));
  for (int i=0; i < pow(4,k); ++i)
    boost::get<0>(kmers[i]) = i;

  // count the number of kmers, only single strand considered here
  for (int i=0; i < lines; ++i) {
    const std::string& line = sequences[i];
    for (int j=0; j < m; ++j) {
      boost::get<1>(kmers[dna_to_number<size_t>(line.substr(j,k))])++;
      boost::get<1>(kmers[dna_to_number<size_t>(reverse_complement(line.substr(j,k)))])++;
    }
  }

  key_sort(kmers.begin(), kmers.end(), myget2); // sort by count

  int a = ceil(0.25*pow(4,k));
  int b = floor(0.75*pow(4,k));

  int sum=0;
  for (int i=a; i <= b; ++i)    // number of kmers that have relative frequency in [0.25,0.75]
    sum += boost::get<1>(kmers[i]);
  
  return sum;
}




double
jussi_lambda(const std::vector<std::string>& s0, 
	    const std::vector<std::string>& s1)
{
  int L = s0[0].length();
  int k = 8;
  double f0=find_nonspecific_fraction(s0);
  double f1=find_nonspecific_fraction(s1);
  int size0 = s0.size()*(L-k+1);
  int size1 = s1.size()*(L-k+1);
  printf("\t\tTwo middle quads\tAll 8-mers\n");
  printf("\tDataset 1:\t%i\t%i\n", (int)f0, size0);
  printf("\tDataset 2:\t%i\t%i\n", (int)f1, size1);
  
  f0 /= size0;
  f1 /= size1;
  printf("\tFreqs 0: %f Freqs 1: %f\n", f0, f1);

  double lambda=f1/f0;
  
  return lambda;
}


// returns set of sequences of length L that do not contain a subsequence
// close to the consensus by Hamming distance,
// current limit is ???
std::set<std::string>
get_unspecific_set_helper(const std::vector<std::string>& dataset, const std::string& consensus)
{
  std::set<std::string> result;

  int limit= 3; //ceil(consensus.length()/1.5);

  std::string reverse_consensus = reverse_complement(consensus);

  // find all unspecific sequences in dataset
  for (int i=0; i < dataset.size(); ++i) {
    if (min_hamming_distance(dataset[i], consensus) >= limit && 
	(not use_two_strands || min_hamming_distance(dataset[i], reverse_consensus) >= limit))
      result.insert(dataset[i]);
  }

  return result;
}

std::set<std::string>
get_unspecific_set(const std::vector<std::string>& dataset1, 
		   const std::vector<std::string>& dataset2, 
		   const std::string& consensus)
{
  
  std::set<std::string> result1 = get_unspecific_set_helper(dataset1, consensus);
  std::set<std::string> result2 = get_unspecific_set_helper(dataset2, consensus);

  std::set<std::string> result;

  printf("Dataset1 distinct unspecific sequences %zu\n", result1.size());
  printf("Dataset2 distinct unspecific sequences %zu\n", result2.size());
  //set_intersection(result1.begin(), result1.end(), result2.begin(), result2.end(), inserter(result, result.begin()));
  set_union(result1.begin(), result1.end(), result2.begin(), result2.end(), inserter(result, result.begin()));
  //printf("Result dataset unspecific %zu\n", result.size());
  return result;
}

int
count_unspecific(const std::vector<std::string>& dataset, const  std::set<std::string>& unspecific)
{
  int count = 0;
  for (int i=0; i < dataset.size(); ++i) {
    if (unspecific.find(dataset[i]) != unspecific.end())
      ++count;
  }
  return count;
}

// compute lambda as the ratio of frequencies of an unspecific set
// between two generations, unspecific set is the union of
// two unspecific sets of the generations
double
hamming_lambda(const std::vector<std::string>& data1, 
		const std::vector<std::string>& data2,
		const std::string& consensus, int generation)
{
  //const std::vector<std::string>& data1=datasets[generation-1];
  //const std::vector<std::string>& data2=datasets[generation];

  std::set<std::string> unspecific_set = 
    get_unspecific_set(data1, data2, consensus);

  printf("Size of unspecific set of generations %i and %i is %zu\n", 
	 generation-1, generation, unspecific_set.size());


  printf("Dataset1: Unspecific %i  Total %zu\n", 
	 count_unspecific(data1, unspecific_set), data1.size());
  printf("Dataset2: Unspecific %i  Total %zu\n", 
	 count_unspecific(data2, unspecific_set), data2.size());

  double unspec1 = count_unspecific(data1, unspecific_set) / (double)data1.size();
  double unspec2 = count_unspecific(data2, unspecific_set) / (double)data2.size();
  printf("Unspecific frequencies, generation %i: %f   generation %i: %f   ratio(lambda): %f\n", 
	 generation-1, unspec1, generation, unspec2, unspec2/unspec1);


  return unspec2/unspec1;
}
