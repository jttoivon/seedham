#include "matrix.hpp"
#include <string>
#include <vector>


// Position dependent first order model for a site of length k
class dinuc_model
{
public:


  dinuc_model(const std::vector<std::string>& sequences, const std::string& seed, int n=2);
  dinuc_model(const std::string& filename);

  dinuc_model() { }

  void
  init(const dmatrix& dm_);     // initialize the position dependent first order model using dinucleotide count array

  double
  cond(int i, int c1, int c2) const;  // If character in position i is c1, the what is the probability of getting c2 in the next pos 

  int
  length() const;   // returns k, the width of the binding site

  void
  print() const;    // Prints the counts in a 16 x (k-1) matrix

  double
  score(const std::string& s, int start_pos = 0) const ;

  dmatrix dm;  // 16 x (k-1)
  dmatrix ip;  // 4 x k
  dmatrix cp;  // 16 x (k-1)
};

dinuc_model
reverse_complement(const dinuc_model& dm);

dinuc_model
pwm_to_dinucleotide(const dmatrix& pwm);

dmatrix
dinucleotide_counts_suffix_array(const std::string& seed, const std::vector<std::string>& sequences, int n);

dmatrix
dinucleotide_counts_scan(const std::string& seed, const std::vector<std::string>& sequences, int n);
