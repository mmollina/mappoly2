/*
 MAPpoly: a package to construct genetic maps in autopolyploids
 Copyright (C) 2014-2023 Marcelo Mollinari

 This file is part of MAPpoly.

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

/*
 File: combinatorics.cpp

 Description: The functions involve various calculations related to genetics
 and combinatorics. The code also includes the use of Rcpp, which allows for
 integration between R and C++.

 Functions Written by Marcelo Mollinari.

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version:       2022
 Last update: May 25, 2023
 */

// [[Rcpp::plugins(cpp11)]]
#include <R.h>
#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include "combinatorics.h"
#include <math.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <unordered_set>
#define THRESH 200.0
using namespace Rcpp;
using namespace std;

/*
This function calculates the binomial coefficient of `n` and `k`.
*/
int nChoosek(int n, int k)
{
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;
  int result = n;
  for( int i = 2; i <= k ; ++i )
  {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}

/*
Given the ploidy level and two indices representing two genotypes,
this function calculates the number of recombinant events between
the two genotypes.
*/
int n_rec_given_genk_and_k1(int ploidy, int index1, int index2)
{
  int i, result = 0;
  std::vector<bool> vec1(ploidy), vec2(ploidy);
  std::fill(vec1.begin(), vec1.end()-ploidy/2, false);
  std::fill(vec2.begin(), vec2.end()-ploidy/2, false);
  vec1=get_boolean_vec_from_lexicographical_index(ploidy, index1);
  vec2=get_boolean_vec_from_lexicographical_index(ploidy, index2);
  for(i=0; i < ploidy; i++)
  {
    if((vec1[i]+vec2[i]) == 2)
    {
      result++;
    }
  }
  result = ploidy/2 - result;
  return result;
}

/*
This function generates a Boolean vector representing a lexicographical
combination of a given index at a certain ploidy level.
*/
std::vector <bool> get_boolean_vec_from_lexicographical_index(int ploidy, int index)
{
  int i, j, increment, sentinel;
  std::vector<bool> vec(ploidy+1);
  i=0;
  j=1;
  increment=0;
  sentinel=0;
  std::fill(vec.begin(), vec.end(), 0);
  while(sentinel < ploidy/2)
  {
    if(index > nChoosek((ploidy-j), (ploidy/2 - (i+1))) + increment)
    {
      vec[j-1]=0;
      increment += nChoosek((ploidy-j), (ploidy/2 - (i+1)));
    }
    else
    {
      vec[j-1]=1;
      i++;
    }
    sentinel += vec[j-1];
    j++;
  }
  return vec;
}

/*
This function checks if a permutation vector `v` is valid, given matrix `H`
and vector `d`.
*/
bool valid_permutation(NumericVector v, NumericMatrix H, NumericVector d) {
  for (int i = 0; i < H.nrow(); ++i) {
    if (!R_IsNA(d[i])) {
      double dot_product = 0.0;
      for (int j = 0; j < H.ncol(); ++j) {
        dot_product += v[j] * H(i, j);
      }
      if (dot_product != d[i]) {
        return false;
      }
    }
  }
  return true;
}


/*
A recursive function used to find all valid permutations of a given vector `v`.
*/
void find_permutations(NumericVector v, NumericMatrix H, NumericVector d, int start_idx, int num_ones, NumericMatrix &results, int &resultIdx) {
  if (num_ones == 0) {
    if (valid_permutation(v, H, d)) {
      results(resultIdx, _) = v;
      resultIdx++;
    }
    return;
  }
  for (int i = start_idx; i < v.length(); ++i) {
    v[i] = 1;
    find_permutations(v, H, d, i + 1, num_ones - 1, results, resultIdx);
    v[i] = 0;
  }
}

/*
This function uses `find_permutations()` to find all valid permutations, then
returns them as a `NumericMatrix`.
*/
// [[Rcpp::export]]
NumericMatrix find_valid_permutations(NumericMatrix H, NumericVector d, int x) {
  NumericVector v(H.ncol(), 0.0);
  int maxPermutations = nChoosek(H.ncol(), x);
  NumericMatrix results(maxPermutations, H.ncol());
  int resultIdx = 0;

  find_permutations(v, H, d, 0, x, results, resultIdx);

  if (resultIdx == 0) {
    return NumericMatrix(0, H.ncol());
  }
  else if (resultIdx < maxPermutations) {
    results = results(Range(0, resultIdx-1), _);
  }

  return results;
}

/*
This function converts a matrix into a string representation.
*/
string mat_to_string(IntegerMatrix mat) {
  string str = "";
  for(int i = 0; i < mat.nrow(); i++) {
    for(int j = 0; j < mat.ncol(); j++) {
      str += to_string(mat(i,j));
    }
  }
  return str;
}

/*
This function generates all unique permutations of columns for a given
 binary matrix.
*/
unordered_set<string> get_all_permutations(IntegerMatrix mat) {
  unordered_set<string> permutations;
  int ncols = mat.ncol();

  // Generating initial state
  IntegerVector column_order(ncols);
  for(int i = 0; i < ncols; i++) {
    column_order[i] = i;
  }

  do {
    // Creating permuted matrix
    IntegerMatrix perm_mat(mat.nrow(), ncols);
    for(int i = 0; i < ncols; i++) {
      perm_mat(_,i) = mat(_,column_order[i]);
    }
    permutations.insert(mat_to_string(perm_mat));
  } while(std::next_permutation(column_order.begin(), column_order.end()));

  return permutations;
}

/*
This function removes any matrices from a list that have the same permutations
as another matrix in the list.
*/
// [[Rcpp::export]]
List filter_matrices(List mat_list) {
  int n = mat_list.size();
  vector<bool> is_unique(n, true);

  for(int i = 0; i < n; i++) {
    if(!is_unique[i]) continue;
    IntegerMatrix mat_i = as<IntegerMatrix>(mat_list[i]);
    unordered_set<string> perms_i = get_all_permutations(mat_i);
    for(int j = i+1; j < n; j++) {
      if(!is_unique[j]) continue;
      IntegerMatrix mat_j = as<IntegerMatrix>(mat_list[j]);
      if(perms_i.count(mat_to_string(mat_j))) {
        is_unique[j] = false;
      }
    }
  }

  // Creating list of unique matrices
  List unique_matrices;
  for(int i = 0; i < n; i++) {
    if(is_unique[i]) {
      unique_matrices.push_back(mat_list[i]);
    }
  }

  return unique_matrices;
}

//end of file
