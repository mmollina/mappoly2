#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <string_view>

using namespace std;
using namespace Rcpp;

IntegerMatrix combn(NumericVector x, int m);
IntegerVector rep_each(IntegerVector x, int n);
IntegerVector rep_len(IntegerVector x, int n);
IntegerMatrix expand_grid(IntegerVector v1, IntegerVector v2);
IntegerVector which(LogicalVector x);
IntegerVector concatenate_vectors(IntegerVector x1, IntegerVector x2);
List calculate_L_and_initialize_H(int n_fullsib_pop,
                                  int n_mrk,
                                  int n_ind,
                                  NumericVector ploidy_p1,
                                  NumericVector ploidy_p2);
NumericMatrix retainUniqueAndSortByLastColumn(NumericMatrix x);
List vs_multiallelic_Rcpp(List PH,
                          List GENO,
                          NumericMatrix pedigree);
List vs_biallelic_Rcpp(List PH,
                       IntegerMatrix G,
                       NumericMatrix pedigree);

List vs_biallelic_single_Rcpp(NumericMatrix PH,
                              IntegerMatrix G);

List hmm_vectors(List input_list);
#endif /* UTILS_H */
