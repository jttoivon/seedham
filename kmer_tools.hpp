#ifndef KMER_TOOLS_HPP
#define KMER_TOOLS_HPP


#include "type.hpp"

#include <vector>
#include <string>

#include <boost/unordered_map.hpp>
#include <boost/multi_array.hpp>
#include "unordered_map.hpp"


extern bool use_middle_gap;




/*
namespace std {
  std::size_t 
  hash_value(const myuint128& input);
}
*/

// Hash function to make int128 work with boost::hash.
struct hash128
    : std::unary_function<myuint128, std::size_t>
{
  std::size_t 
  operator()(const myuint128& input) const
  {
    boost::hash<unsigned long long> hasher;
    unsigned long long hashVal = 0;
    unsigned long long mask = (1ull<<32) - 1;
    
    for(int i = 0; i < 3; ++i)
      {
    	hashVal *= 37;
    	hashVal += (mask & (input>>(32*i)));
    }

    return hasher(hashVal);
  }
};

#define USE_HASH 1


//typedef unsigned long long int code_t;   // These code bit representations of DNA sequences
typedef myuint128 code_t;
//typedef uint64_t code_t;

typedef unsigned int count_t;


#ifdef USE_HASH
//typedef boost::unordered_map<code_t, count_t> count_container_t;
typedef my_unordered_map<code_t, count_t, hash128> count_container_t;
#else
typedef std::vector<count_t> count_container_t;   // index type of std::vector is size_t
#endif

typedef boost::multi_array<count_container_t, 2> gapped_count_container_t; 


void
get_kmer_counts(const std::vector<std::string>& sequences, int k,
		count_container_t& count, bool use_two_strands = true, bool count_palindromes_twice = false);

void
get_gapped_kmer_counts(const std::vector<std::string>& sequences, int k, int gaps,
		       gapped_count_container_t& count);

void
get_gapped_kmer_counts_fast(const std::vector<std::string>& sequences, int k, int max_gap,
			    gapped_count_container_t& count);



class dummy  // this functor reverses the 4-mer stored in a byte
{
public:

  dummy() : v(256) {        // initialise the v vector
    unsigned char mask[4] = {0xc0, 0x30, 0x0c, 0x03};
    for (int i=0; i < 256; ++i) {
      unsigned char x0=(i & mask[0])>>6;
      unsigned char x1=(i & mask[1])>>2;
      unsigned char x2=(i & mask[2])<<2;
      unsigned char x3=(i & mask[3])<<6;
      v[i] = x0 | x1 | x2 | x3;
    }
  }

  unsigned char
  operator()(unsigned char i) { return v[i]; }

private:
  std::vector<unsigned char> v;
};

extern dummy my_reverse_4_2bit_string;

uint32_t
reverse32_2bitstring(uint32_t u, int k);

uint64_t
reverse64_2bitstring(uint64_t u, int k);

template <typename T>
T
reverse_2bitstring(T u, int k)
{
  unsigned char* p=reinterpret_cast<unsigned char*>(&u);
  int bytes = sizeof(T);
  int bits = bytes * 8;
  
  for (int i=0; i < bytes; ++i)
    p[i]=my_reverse_4_2bit_string(p[i]);
  for (int i=0; i < bytes/2; ++i)
    std::swap(p[i], p[bytes-i-1]);
    
  return u>>(bits-2*k);
}

big_int
reverse_complement_2bitstring(big_int c, int l);



big_int
remove_gap_from_bitstring(big_int c, int k, int gap, int pos);

class kmers
{
public:

  kmers(const std::vector<std::string>& sequences, int maxk_);

  unsigned
  count(int k, big_int code) const ;


  unsigned
  count(int k, const std::string& str) const ;

  double 
  probability(int k, big_int code) const ;


  double
  probability(int k, const std::string& str) const;

  int get_maxk() const 
  { return maxk; }

private:
  int maxk;
  std::vector<count_container_t> counts;
  std::vector<unsigned int> total_counts;
};

#endif // KMER_TOOLS_HPP
