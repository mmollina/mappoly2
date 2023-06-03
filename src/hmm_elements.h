#ifndef HMM_FUNCTIONS_H
#define HMM_FUNCTIONS_H

#include <algorithm>
#include <iostream>
#include <vector>
#include "combinatorics.h"
#include "hmm_elements.h"
#include <math.h>
#include <Rmath.h>
#include <Rcpp.h>
#include <R_ext/PrtUtil.h>

#define THRESH 200.0

using namespace std;
using namespace Rcpp;

double prob_k1_given_k_l_m(int ploidy, int l, double rf);
double log_prob_k1_given_k_l_m(int ploidy, int l, double rf);
std::vector<std::vector<double> > rec_num(int ploidy);
std::vector<std::vector<double> > transition(int ploidy, double rf);
std::vector<std::vector<double> > log_transition(int ploidy, double rf);
std::vector<double> forward_emit(std::vector<double>& fk,
                                 std::vector<int>& ik,
                                 std::vector<int>& ik1,
                                 std::vector<double>& emit,
                                 std::vector<std::vector<double> >& T1,
                                 std::vector<std::vector<double> >& T2);
std::vector<double> forward_emit_single_parent(int m,
                                               std::vector<double>& fk,
                                               std::vector<int>& ik,
                                               std::vector<int>& ik1,
                                               std::vector<double>& emit,
                                               std::vector<std::vector<double> >& T);
std::vector<double> backward_emit(std::vector<double>& fk1,
                                  std::vector<int>& ik,
                                  std::vector<int>& ik1,
                                  std::vector<double>& emit,
                                  std::vector<std::vector<double> >& T1,
                                  std::vector<std::vector<double> >& T2);
std::vector<double> backward_emit_single_parent(int m,
                                                std::vector<double>& fk1,
                                                std::vector<int>& ik,
                                                std::vector<int>& ik1,
                                                std::vector<double>& emit,
                                                std::vector<std::vector<double> >& T);
double addlog(double a, double b);

#endif // HMM_FUNCTIONS_H
