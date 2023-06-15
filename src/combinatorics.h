#ifndef COMBINATORICS_H
#define COMBINATORICS_H

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

int nChoosek(int n, int k);
int n_rec_given_genk_and_k1(int ploidy, int index1, int index2);
std::vector <bool> get_boolean_vec_from_lexicographical_index(int ploidy, int index);
bool valid_permutation(NumericVector v, NumericMatrix H, NumericVector d);
void find_permutations(NumericVector v, NumericMatrix H, NumericVector d, int start_idx, int num_ones, NumericMatrix &results, int &resultIdx);
NumericMatrix find_valid_permutations(NumericMatrix H, NumericVector d, int x);
string mat_to_string(IntegerMatrix mat);
unordered_set<string> get_all_permutations(IntegerMatrix mat);
List filter_matrices(List mat_list);
IntegerMatrix combn(NumericVector x, int m);
IntegerMatrix get_all_combinations(int ploidy, int dose);
NumericMatrix make_mat(double x, int nrow, int ncol);
IntegerVector rep_each(IntegerVector x, int n);
IntegerVector rep_len(IntegerVector x, int n);
IntegerMatrix expand_grid(IntegerVector v1, IntegerVector v2);
IntegerVector which(LogicalVector x);
IntegerVector concatenate_vectors(IntegerVector x1, IntegerVector x2);
List calculate_L_and_initialize_H(int n_fullsib_pop, int n_mrk, int n_ind, NumericVector ploidy_p1, NumericVector ploidy_p2);
NumericMatrix retainUniqueAndSortByLastColumn(NumericMatrix x);
NumericVector calculate_hmm_combinatorial_products(NumericVector h,int ploidy_p1,int ploidy_p2);
#endif // COMBINATORICS_H
