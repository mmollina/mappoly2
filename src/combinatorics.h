#ifndef COMBINATORICS_H
#define COMBINATORICS_H

#include <R.h>
#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <unordered_set>

#define THRESH 200.0

using namespace Rcpp;
using namespace std;

int nChoosek(int n, int k);

int n_rec_given_genk_and_k1(int ploidy, int index1, int index2);

std::vector <bool> get_boolean_vec_from_lexicographical_index(int ploidy, int index);

bool valid_permutation(NumericVector v, NumericMatrix H, NumericVector d);

void find_permutations(NumericVector v, NumericMatrix H, NumericVector d, int start_idx, int num_ones, NumericMatrix &results, int &resultIdx);

NumericMatrix find_valid_permutations(NumericMatrix H, NumericVector d, int x);

string mat_to_string(IntegerMatrix mat);

unordered_set<string> get_all_permutations(IntegerMatrix mat);

List filter_matrices(List mat_list);

#endif /* COMBINATORICS_H */
