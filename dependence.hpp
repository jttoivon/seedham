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
