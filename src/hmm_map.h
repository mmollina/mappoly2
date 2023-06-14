#ifndef HMM_MAP_H
#define HMM_MAP_H

#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <string_view>
#include "combinatorics.h"
#include "hmm_elements.h"
#include "utils.h"

using namespace std;
using namespace Rcpp;

double calc_loglike(std::vector<std::vector<std::vector<int> > > v,
                    std::vector<std::vector<std::vector<double> > > emit,
                    std::vector<double> rf_vec,
                    NumericVector ploidy_p1,
                    NumericVector ploidy_p2);

double calc_loglike_single(std::vector<std::vector<std::vector<int> > > v,
                           std::vector<std::vector<std::vector<double> > > emit,
                           std::vector<double> rf_vec,
                           int ploidy);

List main_hmm_full(NumericVector ploidy_p1,
                   NumericVector ploidy_p2,
                   std::vector<std::vector<std::vector<int> > > v,
                   std::vector<std::vector<std::vector<double> > > e,
                   NumericVector init_rf,
                   double tol,
                   bool verbose,
                   bool detailed_verbose,
                   bool ret_H0);

List main_hmm_full_single(int ploidy,
                          std::vector<std::vector<std::vector<int> > > v,
                          std::vector<std::vector<std::vector<double> > > e,
                          NumericVector init_rf,
                          double tol,
                          bool verbose,
                          bool detailed_verbose,
                          bool ret_H0);

List est_hmm_map_biallelic(List PH,
                           IntegerMatrix G,
                           NumericMatrix pedigree,
                           NumericVector rf,
                           double err,
                           bool verbose,
                           bool detailed_verbose,
                           double tol,
                           bool ret_H0);

List est_hmm_map_biallelic_single(NumericMatrix PH,
                                  IntegerMatrix G,
                                  NumericVector rf,
                                  double err,
                                  bool verbose,
                                  bool detailed_verbose,
                                  double tol,
                                  bool ret_H0);

#endif // HMM_MAP_H
