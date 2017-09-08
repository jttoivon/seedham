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
#include "data.hpp"
#include <cstdio>


class packed_string
{
public:

  class my_wrapper
  {
  public:
    my_wrapper(big_int& b, int s) : bits(b), shift(s) {}

    // for reading the char
    operator char() 
    {
      return (bits >> shift) & 3;
    }

    // for setting the char
    my_wrapper&
    operator=(char ch)
    {
      big_int c = (ch & 3) << shift;
      big_int mask = 3;
      mask = ~(mask << shift);

      bits = (bits & mask) | c;

      return *this;
    }


  private:
    big_int& bits;
    int shift;
  };


  packed_string() : bits(0), k(0) {}

  packed_string(big_int b, int k_) : bits(b), k(k_) {}

  packed_string(const std::string& s);

 
  int 
  length() const;


  bool
  operator==(const packed_string& rhs) const
  {
    return k == rhs.k and bits == rhs.bits;
  }

  bool
  operator!=(const packed_string& rhs) const
  {
    return not (*this == rhs);
  }

  char
  operator[](size_t i) const ;

  my_wrapper
  operator[](size_t i);

  std::string 
  to_string() const
  {
    return number_to_dna(bits, k);
  }

  packed_string&
  erase(int pos, int dummy_count);

  void
  debug()
  {
    int number_of_bytes = sizeof(bits);
    int number_of_bits = 8*number_of_bytes;
    std::string bit_repr(number_of_bits, '0');
    int pos=0;
    unsigned char* start = (unsigned char*)&bits;
    for (int i = number_of_bytes-1; i >= 0; --i) {
      unsigned char c =start[i];
      for (int j = 7; j >= 0; --j) {
	unsigned char mask = 1 << j;
	bit_repr[pos] = c & mask ? '1' : '0';
	++pos;
	//printf("pos=%i\n", pos);
      }
    }
    printf("bits: %s\n", bit_repr.c_str());
    printf("k: %i\n", k);
    printf("max_len: %i\n", max_len);
  }

  packed_string&
  insert(int pos, int dummy_count, char ch);

  big_int
  get_bits() const
  { 
    return bits;
  }

private:
  big_int bits;
  int k;
  static const int max_len = sizeof(big_int) * 4;


};





