#include "packed_string.hpp"

#include <cassert>

packed_string::packed_string(const std::string& s) : bits(0), k(0) 
{
  k = s.length();
  assert(k <= max_len);
  bits = dna_to_number(s);
}

int 
packed_string::length() const
{
  return k;
}

char
packed_string::operator[](size_t i) const 
{
  assert(i < k);
  return (bits >> (2*(k-i-1))) & 3;
}

packed_string::my_wrapper
packed_string::operator[](size_t i)
{
  assert(i < k);
  return my_wrapper(bits, 2*(k-i-1));
}

packed_string&
packed_string::erase(int pos, int dummy_count)
{
  assert(pos < k);
  big_int suffix_mask = 1;
  suffix_mask = (suffix_mask << (2*(k-pos-1))) - 1;
  big_int prefix_mask = ~suffix_mask;
  bits = ((bits >> 2) & prefix_mask) | (bits & suffix_mask);
  --k;
  return *this;
}

packed_string&
packed_string::insert(int pos, int dummy_count, char ch)
{
  assert(pos <= k);
  big_int suffix_mask = 1;
  suffix_mask = (suffix_mask << (2*(k-pos))) - 1;
  big_int prefix_mask = ~suffix_mask;
  big_int c = ch & 3;
  c <<= (2*(k-pos));
  bits = ((bits & prefix_mask) << 2) | c | (bits & suffix_mask);
  ++k;
  assert(k <= max_len);

  return *this;
}
