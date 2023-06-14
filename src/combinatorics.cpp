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

//This function calculates the binomial coefficient of `n` and `k`.
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

//Given the ploidy level and two indices representing two genotypes, this function calculates the number of recombinant events between the two genotypes.
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

//This function generates a Boolean vector representing a lexicographical combination of a given index at a certain ploidy level.
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

//This function checks if a permutation vector `v` is valid, given matrix `H` and vector `d`.
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

//A recursive function used to find all valid permutations of a given vector `v`.
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

// Finds all valid permutations, then returns them as a `NumericMatrix`.
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

// This function converts a matrix into a string representation.
string mat_to_string(IntegerMatrix mat) {
  string str = "";
  for(int i = 0; i < mat.nrow(); i++) {
    for(int j = 0; j < mat.ncol(); j++) {
      str += to_string(mat(i,j));
    }
  }
  return str;
}

// This function generates all unique permutations of columns for a given binary matrix.
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

// This function removes any matrices from a list that have the same permutations as another matrix in the list.
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

// Helper function to calculate the combinations
IntegerMatrix combn(NumericVector x, int m) {
  int n = x.size();
  if (m > n) return IntegerMatrix(0);
  if (m == n) return IntegerMatrix(x);

  int count = Rf_choose(n, m);
  IntegerMatrix result(m, count);

  std::vector<int> indices(m);
  std::iota(indices.begin(), indices.end(), 0);

  for (int col = 0; col < count; ++col) {
    for (int row = 0; row < m; ++row) {
      result(row, col) = x[indices[row]];
    }

    for (int i = m - 1; i >= 0; --i) {
      if (++indices[i] <= n - m + i) {
        while (++i < m) {
          indices[i] = indices[i - 1] + 1;
        }
        break;
      }
    }
  }

  return result;
}

// Creates a matrix of dimensions (nrow x ncol) populated by double x
NumericMatrix make_mat(double x, int nrow, int ncol) {
  NumericMatrix res(nrow,ncol);
  for(int i = 0; i < nrow; i++)
    for(int j = 0; j < ncol; j++)
      res(i,j) = x;
  return res;
}

// Function for rep_each
IntegerVector rep_each(IntegerVector x, int n) {
  int sz = x.size() * n;
  IntegerVector res(sz);
  for (int i = 0; i < x.size(); ++i) {
    std::fill_n(res.begin() + i * n, n, x[i]);
  }
  return res;
}

// Function for rep_len
IntegerVector rep_len(IntegerVector x, int n) {
  IntegerVector res(n);
  std::copy_n(x.begin(), std::min<int>(n, x.size()), res.begin());
  for (int i = x.size(); i < n; i += x.size()) {
    std::copy_n(x.begin(), std::min<int>(n - i, x.size()), res.begin() + i);
  }
  return res;
}

// Helper function for expand_grid
IntegerMatrix expand_grid(IntegerVector v1, IntegerVector v2) {
  int n_v1 = v1.size();
  int n_v2 = v2.size();
  IntegerVector res_v1 = rep_len(v1, n_v1 * n_v2); // Adjust the rep_len() usage
  IntegerVector res_v2 = rep_each(v2, n_v1); // Adjust the rep_each() usage
  return cbind(res_v2, res_v1);
}

// Helper function for which
IntegerVector which(LogicalVector x) {
  IntegerVector idx = seq_len(x.size());
  IntegerVector out;
  for (int i = 0; i < x.size(); ++i) {
    if (x[i]) {
      out.push_back(idx[i]);
    }
  }
  return out;
}

// Helper function for concatenating vectors
IntegerVector concatenate_vectors(IntegerVector x1, IntegerVector x2) {
  int n1 = x1.size();
  int n2 = x2.size();
  IntegerVector result(n1 + n2);
  std::copy(x1.begin(), x1.end(), result.begin());
  std::copy(x2.begin(), x2.end(), result.begin() + n1);
  return result;
}

// Helper function for calculate L and initialize H
List calculate_L_and_initialize_H(int n_fullsib_pop,
                                  int n_mrk,
                                  int n_ind,
                                  NumericVector ploidy_p1,
                                  NumericVector ploidy_p2) {
  List L(n_fullsib_pop);
  for (int i = 0; i < n_fullsib_pop; i++) {
    int ngam1 = R::choose(ploidy_p1[i], ploidy_p1[i] / 2);
    int ngam2 = R::choose(ploidy_p2[i], ploidy_p2[i] / 2);
    IntegerMatrix S = expand_grid(seq_len(ngam2) - 1, seq_len(ngam1) - 1);
    L[i] = S;
  }
  List H(n_mrk); //H[[n_mrk]][[n_ind]]
  List E(n_mrk);
  for (int k = 0; k < n_mrk; k++) {
    H[k] = List(n_ind);
    E[k] = List(n_ind);
  }
  return List::create(Named("L") = L,
                      Named("H") = H,
                      Named("E") = E);
}

// Helper function for unique rows
NumericMatrix retainUniqueAndSortByLastColumn(NumericMatrix x) {
  std::set<std::vector<double>> uniqueRows;
  int numRows = x.nrow();
  int numCols = x.ncol();

  // Iterate over each row
  for (int i = 0; i < numRows; i++) {
    std::vector<double> row;

    // Extract values from the current row
    for (int j = 0; j < numCols; j++) {
      row.push_back(x(i, j));
    }

    // Check if the row is already present
    if (uniqueRows.find(row) == uniqueRows.end()) {
      uniqueRows.insert(row); // Insert the row into the set
    }
  }

  // Sort the unique rows based on the last column value
  std::vector<std::vector<double>> sortedRows(uniqueRows.begin(), uniqueRows.end());
  std::sort(sortedRows.begin(), sortedRows.end(), [numCols](const std::vector<double>& a, const std::vector<double>& b) {
    return a[numCols - 1] < b[numCols - 1];
  });

  // Create a new matrix to store the sorted and unique rows
  NumericMatrix result(sortedRows.size(), numCols);
  for (int i = 0; i < sortedRows.size(); i++) {
    const std::vector<double>& row = sortedRows[i];
    for (int j = 0; j < numCols; j++) {
      result(i, j) = row[j];
    }
  }

  return result;
}

// Helper function to calculate genotype probabilities given homolog probabilities (h)
NumericVector calculate_hmm_combinatorial_products(NumericVector h,
                                                   int ploidy_p1,
                                                   int ploidy_p2){
  NumericVector v1(ploidy_p1),  v2(ploidy_p2);
  for(int i=0; i < ploidy_p1; i++)
    v1[i] = i;
  for(int i=0; i < ploidy_p2; i++)
    v2[i] = i + ploidy_p1;
  IntegerMatrix u1 = combn(v1, ploidy_p1/2);
  IntegerMatrix u2 = combn(v2, ploidy_p2/2);
  NumericVector res(u1.ncol() * u2.ncol());
  NumericVector w1(u1.ncol());
  for(int i=0; i < u1.ncol(); i++){
    double temp = 1.0;
    for(int j=0; j < u1.nrow(); j++){
      temp *= h[u1(j,i)];
    }
    w1[i] = temp;
  }
  NumericVector w2(u2.ncol());
  for(int i=0; i < u2.ncol(); i++){
    double temp = 1.0;
    for(int j=0; j < u2.nrow(); j++){
      temp *= h[u2(j,i)];
    }
    w2[i] = temp;
  }
  for(int i=0; i < u1.ncol(); i++)
    for(int j=0; j < u2.ncol(); j++)
      res[j + (i*u2.ncol())] = w1[i] * w2[j];
  return(res);
}

//end of file
