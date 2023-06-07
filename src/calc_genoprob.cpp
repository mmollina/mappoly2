/*
 MAPpoly-2.0: a package to construct genetic maps in autopolyploids
 Copyright (C) 2014-2023 Marcelo Mollinari

 This file is part of MAPpoly-2.0.

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
 File: calc_genoprob.cpp

 Description:

 Functions Written by Marcelo Mollinari.

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version:       2022
 Last update: Jun 6, 2023
 */

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


// [[Rcpp::export]]
List calc_genoprob_biallelic(List PH,
                             IntegerMatrix G,
                             NumericMatrix pedigree,
                             NumericVector rf,
                             double err) {
  NumericVector ploidy_p1 = pedigree(_,2)/2 - 1;
  NumericVector ploidy_p2 = pedigree(_,3)/2 - 1;
  int n_ind = pedigree.nrow();
  NumericMatrix temp_phase_mat = PH[0];
  int n_mrk = temp_phase_mat.nrow();
  std::vector<int> pl{2,4,6};
  int k, k1;

  // Getting the maximum ploidy level among founders
  int mpi1 = max(ploidy_p1);
  int mpi2 = max(ploidy_p2);
  int max_ploidy_id;
  if(mpi1 >= mpi2)
    max_ploidy_id = mpi1;
  else
    max_ploidy_id = mpi2;

  // HMM states that should be visited given the phase of
  // the founders, genotype of the offspring and pedigree
  List result = vs_biallelic_error(PH, G, pedigree, err, 0);

  List ve = hmm_vectors(result);
  std::vector<std::vector<std::vector<int> > > v = ve["v"];
  std::vector<std::vector<std::vector<double> > > e = ve["e"];

  //Initializing alpha, beta, and gamma
  std::vector<std::vector<std::vector<double> > > alpha(n_ind);
  std::vector<std::vector<std::vector<double> > > beta(n_ind);
  std::vector<std::vector<std::vector<double> > > gamma(n_ind);
  for(int ind=0; ind < n_ind; ind++)
  {
    for(int i=0; i < n_mrk; i++)
    {
      std::vector<double> temp3(v[i][ind].size()/2);
      alpha[ind].push_back(temp3);
      beta[ind].push_back(temp3);
      gamma[ind].push_back(temp3);
    }
  }

  //Initializing transition matrices
  std::vector< std::vector< std::vector< std::vector<double> > > >T;
  for(int j=0; j <= max_ploidy_id; j++)
  {
    std::vector< std::vector< std::vector<double> > > Ttemp;
    for(int i=0; i < n_mrk-1; i++)
    {
      Ttemp.push_back(transition(pl[j], rf[i]));
    }
    T.push_back(Ttemp);
  }
  //Loop over all individuals
  for(int ind=0; ind < n_ind; ind++)
  {
    R_CheckUserInterrupt();
    for(int j=0; (unsigned)j < e[0][ind].size(); j++)
      alpha[ind][0][j] = e[0][ind][j];

    std::fill(beta[ind][n_mrk-1].begin(), beta[ind][n_mrk-1].end(), 1.0);
    //forward-backward
    for(k=1,k1=n_mrk-2; k < n_mrk; k++, k1--)
    {
      std::vector<double> temp4 (v[k][ind].size()/2);
      temp4 = forward_emit(alpha[ind][k-1],
                           v[k-1][ind],
                                 v[k][ind],
                                     e[k][ind],
                                         T[ploidy_p1[ind]][k-1],
                                                          T[ploidy_p2[ind]][k-1]);
      // Normalization to avoid underflow
      // NOTE: The LogSumExp (LSE) method is not used here for efficiency reasons,
      // as it has been observed that this normalization technique performs adequately.
      double zeta = 0;
      for(int j=0; (unsigned)j < temp4.size(); j++)
        zeta = zeta + temp4[j];
      for(int j=0; (unsigned)j < temp4.size(); j++)
      {
        alpha[ind][k][j]=temp4[j]/zeta;
      }
      std::vector<double> temp5 (v[k1][ind].size()/2);
      temp5=backward_emit(beta[ind][k1+1],
                          v[k1][ind],
                               v[k1+1][ind],
                                      e[k1+1][ind],
                                             T[ploidy_p1[ind]][k1],
                                                              T[ploidy_p2[ind]][k1]);
      // Normalization to avoid underflow
      zeta = 0;
      for(int j=0; (unsigned)j < temp5.size(); j++)
        zeta = zeta + temp5[j];
      for(int j=0; (unsigned)j < temp5.size(); j++)
      {
        beta[ind][k1][j]=temp5[j]/zeta;
      }
    }
  } // End loop over individuals

  //Computing gamma
  for(int ind=0; ind < n_ind; ind++)
  {
    for(int k=0; k < n_mrk; k++)
    {
      double w = 0.0;
      for(int j=0; (unsigned)j < alpha[ind][k].size(); j++)
        w += alpha[ind][k][j]*beta[ind][k][j];
      for(int j=0; (unsigned)j < alpha[ind][k].size(); j++)
        gamma[ind][k][j] = alpha[ind][k][j]*beta[ind][k][j]/w;
    }
  }

  int n1 = gamma.size();
  int n2 = gamma[0].size();
  int n3 = gamma[0][0].size();

  List L(n1);

  for(int i = 0; i < n1; i++) {
    NumericMatrix M(n3, n2);
    for(int j = 0; j < n2; j++) {
      for(int k = 0; k < n3; k++) {
        M(k,j) = gamma[i][j][k];
      }
    }
    L[i] = M;
  }
  return L ;
}

// [[Rcpp::export]]
List calc_genoprob_biallelic_single(NumericMatrix PH,
                                    IntegerMatrix G,
                                    NumericVector rf,
                                    double err){
  int n_mrk = G.nrow();
  int n_ind = G.ncol();
  int ploidy = PH.ncol();
  int k, k1;

  // HMM states that should be visited given the phase of
  // the founders, genotype of the offspring and pedigree
  List result = vs_biallelic_error_single(PH, G, err, 0);
  List ve = hmm_vectors(result);
  std::vector<std::vector<std::vector<int> > > v = ve["v"];
  std::vector<std::vector<std::vector<double> > > e = ve["e"];

  //Initializing alpha, beta, and gamma
  std::vector<std::vector<std::vector<double> > > alpha(n_ind);
  std::vector<std::vector<std::vector<double> > > beta(n_ind);
  std::vector<std::vector<std::vector<double> > > gamma(n_ind);
  for(int ind=0; ind < n_ind; ind++)
  {
    for(int i=0; i < n_mrk; i++)
    {
      std::vector<double> temp3(v[i][ind].size());
      alpha[ind].push_back(temp3);
      beta[ind].push_back(temp3);
      gamma[ind].push_back(temp3);
    }
  }

  //Initializing transition matrices
  std::vector< std::vector< std::vector<double> > > T;
  for(int i=0; i < n_mrk-1; i++)
    T.push_back(transition(ploidy, rf[i]));
  //Loop over all individuals
  for(int ind=0; ind < n_ind; ind++)
  {
    R_CheckUserInterrupt();
    for(int j=0; (unsigned)j < e[0][ind].size(); j++)
      alpha[ind][0][j] = e[0][ind][j];
    std::fill(beta[ind][n_mrk-1].begin(), beta[ind][n_mrk-1].end(), 1.0);
    //forward-backward
    for(k=1,k1=n_mrk-2; k < n_mrk; k++, k1--)
    {
      std::vector<double> temp4 (v[k][ind].size());
      temp4 = forward_emit_single_parent(ploidy, alpha[ind][k-1], v[k-1][ind], v[k][ind], e[k][ind], T[k-1]);
      // Normalization to avoid underflow
      // NOTE: The LogSumExp (LSE) method is not used here for efficiency reasons,
      // as it has been observed that this normalization technique performs adequately.
      double zeta = 0;
      for(int j=0; (unsigned)j < temp4.size(); j++)
        zeta = zeta + temp4[j];
      for(int j=0; (unsigned)j < temp4.size(); j++)
        alpha[ind][k][j]=temp4[j]/zeta;
      std::vector<double> temp5 (v[k1][ind].size());
      temp5=backward_emit_single_parent(ploidy, beta[ind][k1+1], v[k1][ind], v[k1+1][ind], e[k1+1][ind], T[k1]);
      // Normalization to avoid underflow
      zeta = 0;
      for(int j=0; (unsigned)j < temp5.size(); j++)
        zeta = zeta + temp5[j];
      for(int j=0; (unsigned)j < temp5.size(); j++)
        beta[ind][k1][j]=temp5[j]/zeta;
    }
  }
  //Computing gamma
  for(int ind=0; ind < n_ind; ind++)
  {
    for(int k=0; k < n_mrk; k++)
    {
      double w = 0.0;
      for(int j=0; (unsigned)j < alpha[ind][k].size(); j++)
        w += alpha[ind][k][j]*beta[ind][k][j];
      for(int j=0; (unsigned)j < alpha[ind][k].size(); j++)
        gamma[ind][k][j] = alpha[ind][k][j]*beta[ind][k][j]/w;
    }
  }
  int n1 = gamma.size();
  int n2 = gamma[0].size();
  int n3 = gamma[0][0].size();

  List L(n1);

  for(int i = 0; i < n1; i++) {
    NumericMatrix M(n3, n2);
    for(int j = 0; j < n2; j++) {
      for(int k = 0; k < n3; k++) {
        M(k,j) = gamma[i][j][k];
      }
    }
    L[i] = M;
  }
  return L ;
}
