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
#include "matrix.hpp"

#include <string>
#include <boost/math/distributions/chi_squared.hpp>

matrix<double>
get_dependence_matrix(const std::vector<std::string>& seqs, const std::string& str, const matrix<double>& m);

matrix<double>
get_dependence_matrix2(const std::vector<std::string>& seqs, const std::string& str, matrix<double> m);

//double
//chi_sq(double param);

inline
double
chi_sq_statistic_helper(double observed, double expected)
{
  assert(expected != 0);
  return pow((observed - expected), 2) / expected;
}

inline
double
G_statistic_helper(double observed, double expected)
{
  assert(expected != 0);
  if (observed == 0)
    return 0;
  else
    return observed * log(observed / expected);
}

class chi_square
{
public:

  chi_square(int d) : df(d), ch(df) {}  // defines the degrees of freedom

  double
  operator()(double param) const {
    return 1.0 - cdf(ch, param);   // area of the right tail
  }

private:
  int df;    // degrees of freedom
  boost::math::chi_squared_distribution<> ch;
};
