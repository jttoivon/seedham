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
#include "multinomial_helper.hpp"

#include <string>
#include <vector>

extern int palindromic_index_limit;
extern counting_type seed_counting;

int
conflict_free_palindromic_index(int hamming_radius);


boost::tuple<std::string,int>
most_common_pattern(const std::vector<std::string>& sequences, int k, std::string seed,
			     bool contains_N, int hamming_radius=0);
