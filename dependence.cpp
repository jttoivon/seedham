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
#include "dependence.hpp"

#include "common.hpp"
#include "parameters.hpp"
#include "data.hpp"
#include "kmp.hpp"
#include "matrix_tools.hpp"
#include "suffix_array_wrapper.hpp"

#include <vector>
#include <cmath>
#include <cassert>



// boost::math::chi_squared_distribution<> ch(16-1);

// double
// chi_sq(double param)
// {
//   return 1.0 - cdf(ch, param);
// }





// double
// chi_sq_statistic2(const std::vector<double>& counts, 
// 		 const std::vector<double>& expectations)
// {
//   assert(counts.size() == expectations.size());
//   double sum=0;
//   for (int i = 0; i < counts.size(); ++i)
//     sum += chi_sq_statistic_helper(counts[i], expectations[i]);

//   //return chi_sq(sum);
//   return sum;
// }

double
matrix_chi_sq_statistic(const matrix<double>& counts, 
			const matrix<double>& expectations)
{
  
  assert(counts.get_columns() == expectations.get_columns());
  assert(counts.get_rows() == expectations.get_rows());
  int rows = counts.get_rows();
  int columns = counts.get_columns();
  double sum=0;
  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < columns; ++j)
      sum += chi_sq_statistic_helper(counts(i,j), expectations(i,j));

  return sum;
}

double
matrix_r_sq_statistic(const std::vector<double>& a, const std::vector<double>& b, 
		      const matrix<double>& ab)
{
  assert(a.size() == b.size());
  assert(ab.get_columns() == ab.get_rows());
  assert(ab.get_rows() == a.size());
  int rows = ab.get_rows();
  int columns = ab.get_columns();
  double sum=0;
  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < columns; ++j)
      sum += (ab(i,j) - a[i]*b[j]) / ((1-a[i])*(1-b[j]));

  return sum;
}

double
find_tnips(const std::string& consensus, int x, int y, const dmatrix& m,
	   const suffix_array& sa1, const suffix_array& sa2)
{
  char nuc[] = "ACGT";
  std::string temp = consensus;
  double total_c=0;
  double total_e=0;
  matrix<double> counts(4,4);
  matrix<double> expectations(4,4);
  std::vector<double> result(4*4,0);

  for (int i=0; i < 4; ++i) {      // iterate through all characters in pos x
    for (int j=0; j < 4; ++j) {      // iterate through all characters in pos y
      temp[x]=nuc[i];
      temp[y]=nuc[j];
      counts(i,j) += sa1.count_iupac(temp);
      if (use_two_strands)
	counts(i,j) += sa2.count_iupac(temp);
      //counts(i, j) += 0.01; // pseudo count
      total_c += counts(i,j);
    }
  }

  if (true) {
    printf("===================================\n");
    printf("X=%i  Y=%i\n", x, y);
    printf("Count matrix:\n");
    counts.print(10, 0);
    printf("Total of counts = %.0f\n", total_c);
  }
  double row_sums[4];
  double column_sums[4];


  //compute row and column sums
  for (int i=0; i < 4; ++i) {
    row_sums[i]=sum(counts.row(i));
    column_sums[i]=sum(counts.column(i));
  }


  assert(fabs(sum(row_sums) - total_c) < 0.0001);
  assert(fabs(sum(column_sums) - total_c) < 0.0001);


  // compute observed and expected frequencies for chi² statistic:
  // sum_i=1^n (O_i - E_i)^2/E_i
  double p_sum = 0;
  for (int i=0; i < 4; ++i) {      // iterate through all characters in pos x
    for (int j=0; j < 4; ++j) {      // iterate through all characters in pos y
      // double px=row_sums[i]/total_c;
      // double py=column_sums[j]/total_c;
      double px=m(i, x);
      double py=m(j, y);
      assert(px <= 1);
      assert(px >= 0);
      assert(py <= 1);
      assert(py >= 0);
      p_sum += (px*py);
      expectations(i,j) = px * py * total_c;
      total_e += expectations(i,j);
   }
  }
  if (true) {
    printf("p_sum = %f\n", p_sum);
    printf("Expectation matrix:\n");
    expectations.print(10,1);
    printf("Total of expectations = %.0f\n", total_e);
    

  }

  add_pseudo_counts(counts, 5.0);
  add_pseudo_counts(expectations, 5.0);

  dmatrix stat(4, 4);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      stat(i, j) = chi_sq_statistic_helper(counts(i,j), expectations(i,j));
  printf("Chi statistics matrix:\n");
  stat.print(10,1);
  printf("Total of stats = %.0f\n", sum(stat));
  printf("===================================\n");

  //  return  chi_sq_statistic(counts.to_vector(), expectations.to_vector());
  return matrix_chi_sq_statistic(counts, expectations);
}

void
find_tnips2(const std::string& consensus, int x, int y, 
	    const std::string& str1, const std::string& str2, 
	    const matrix<double>& m, matrix<double>& result)
{
  char nuc[] = "ACGT";
  std::string temp = consensus;
  double total_c=0;
  //double total_e=0;
  matrix<double> counts(4,4);
  //std::vector<double> result(4*4,0);
  for (int i=0; i < 4; ++i) {      // iterate through all characters in pos x
    for (int j=0; j < 4; ++j) {      // iterate through all characters in pos y
      temp[x]=nuc[i];
      temp[y]=nuc[j];
      counts(i,j)=KMP(str1,temp);
      //products[4*i+j]=m(i,x)*m(j,y);
      if (use_two_strands)
	counts(i,j) += KMP(str2,temp);
      //counts2(i,j)=counts[4*i+j];
      total_c += counts(i,j);
    }
  }

  if (false) {
    printf("===================================\n");
    printf("X=%i  Y=%i\n", x, y);
    printf("Count matrix:\n");
    counts.print(10, 0);
    printf("Total of counts = %.0f\n", total_c);
  }
  double row_sums[4];
  double column_sums[4];


  //compute row and column sums
  for (int i=0; i < 4; ++i) {
    row_sums[i]=sum(counts.row(i));
    column_sums[i]=sum(counts.column(i));
  }


  assert(sum(row_sums) == total_c);
  assert(sum(column_sums) == total_c);

  //  sub_matrix<double> subi(result, x*4, y*4, 4, 4);

  // compute expected counts for chi² statistic:
  // sum_i=1^n (O_i - E_i)^2/E_i
  // double p_sum = 0;
  for (int i=0; i < 4; ++i) {      // iterate through all characters in pos x
    for (int j=0; j < 4; ++j) {      // iterate through all characters in pos y
      assert(m(i,x) != 0);
      assert(m(j,y) != 0);
      matrix<double> sub_counts(2,2);
      sub_counts(0,0) = counts(i,j);
      sub_counts(1,0) = column_sums[i] - counts(i,j);
      sub_counts(0,1) = row_sums[j] - counts(i,j);
      sub_counts(1,1) = total_c - sub_counts(0,0) - sub_counts(1,0) - sub_counts(0,1);

      matrix<double> expectations(2,2);
      expectations(0,0) = m(i,x)     * m(j,y)     * total_c;   // i and j
      expectations(1,0) = (1-m(i,x)) * m(j,y)     * total_c;   // not i but j
      expectations(0,1) = m(i,x)     * (1-m(j,y)) * total_c;   // i but not j
      expectations(1,1) = (1-m(i,x)) * (1-m(j,y)) * total_c;   // not i nor j
      //total_e += expectations(i,j);

      std::vector<double> a(2); a[0] = m(i,x); a[1] = 1-m(i,x);
      std::vector<double> b(2); b[0] = m(j,y); b[1] = 1-m(j,y);

      //subi(i,j) = chi_sq_statistic_helper(counts(i,j), expectations(i,j));
      //result(x*4+i, y*4+j) = matrix_chi_sq_statistic(sub_counts, expectations);
      result(x*4+i, y*4+j) = matrix_r_sq_statistic(a, b, expectations);
   }
  }
  if (false) {
    printf("===================================\n");
  }


  return;
}

// Test of goodness of fit.
//
// Dependent position need not be adjacent
// compute the kxk dependence matrix
matrix<double>
get_dependence_matrix(const std::vector<std::string>& seqs, const std::string& seed, const matrix<double>& m)
{
  
  std::string str1;
  std::string str2;

  str1=join(seqs, '#');
  str2=join_rev(seqs, '#');
  suffix_array sa1(str1);
  suffix_array sa2(str2);

  int k = seed.length();

  matrix<double> result(k,k);
  for (int x=0; x < k; ++x) {
    for (int y=0; y < x; ++y) {
      result(x,y) = find_tnips(seed, x, y, m, sa1, sa2);
    }
  }

  return result;
}

// Test of independence.
//
// Dependent position need not be adjacent
// compute the (kx4)x(kx4) dependence matrix
matrix<double>
get_dependence_matrix2(const std::vector<std::string>& seqs, const std::string& str, matrix<double> m)
{
  
  std::string str1;
  std::string str2;

  str1=join(seqs, '#');
  str2=join_rev(seqs, '#');

  int k = str.length();

  matrix<double> result(4*k,4*k);
  for (int x=0; x < k; ++x) {
    for (int y=0; y < x; ++y) {
      find_tnips2(str, x, y, str1, str2, m, result);
    }
  }

  return result;
}
