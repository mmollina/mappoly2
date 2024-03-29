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
 File: est_hmm_map.cpp

 Description:

 Functions Written by Marcelo Mollinari.

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version:       2022
 Last update: Jun 3, 2023
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
#define SparseThreshold 0.01
using namespace std;
using namespace Rcpp;

double calc_loglike_log(std::vector<std::vector<std::vector<int> > > v,
                        std::vector<std::vector<std::vector<double> > > emit,
                        std::vector<double> rf_vec,
                        NumericVector ploidy_p1,
                        NumericVector ploidy_p2) {

  std::vector<int> pl{2,4,6};
  double loglike = 0.0;

  // Getting the maximum ploidy level among founders
  int mpi1 = max(ploidy_p1);
  int mpi2 = max(ploidy_p2);
  int max_ploidy_id;
  if(mpi1 >= mpi2)
    max_ploidy_id = mpi1;
  else
    max_ploidy_id = mpi2;

  int n_mrk  = v.size();
  int n_ind  = v[0].size();

  //Initializing alpha and beta
  std::vector<std::vector<std::vector<double> > > alpha(n_ind);
  for(int ind=0; ind < n_ind; ind++)
  {
    for(int k=0; k < n_mrk; k++)
    {
      std::vector<double> temp(emit[k][ind].size());
      alpha[ind].push_back(temp);
    }
  }

  //Initializing transition matrices
  std::vector< std::vector< std::vector< std::vector<double> > > >T;
  for(int j=0; j <= max_ploidy_id; j++)
  {
    std::vector< std::vector< std::vector<double> > > Ttemp;
    for(int i=0; i < n_mrk-1; i++)
    {
      Ttemp.push_back(log_transition(pl[j], rf_vec[i]));
    }
    T.push_back(Ttemp);
  }
  //Loop over all individuals
  for(int i=0; i < n_ind; i++)
  {
    R_CheckUserInterrupt();
    for(int j=0; (unsigned)j < emit[0][i].size(); j++)
      alpha[i][0][j] = emit[0][i][j];

    //forward
    for(int k=1; k < n_mrk; k++)
    {
      int ngen_k0 = v[k-1][i].size()/2;
      int ngen_k1 = v[k][i].size()/2;
      for(int s1 = 0; s1 < ngen_k1; s1++)
      {
        alpha[i][k][s1] = alpha[i][k-1][0] +
          T[ploidy_p1[i]][k-1][0][v[k][i][s1]] +
          T[ploidy_p2[i]][k-1][0][v[k][i][s1+ngen_k1]];
        for(int s0 = 0; s0 < ngen_k0; s0++){
          alpha[i][k][s1] =  addlog(alpha[i][k][s1], alpha[i][k-1][s0] +
            T[ploidy_p1[i]][k-1][v[k-1][i][s0]][v[k][i][s1]] +
            T[ploidy_p2[i]][k-1][v[k-1][i][s0+ngen_k0]][v[k][i][s1+ngen_k1]]);
        }
        alpha[i][k][s1] += emit[k][i][s1];
      }
    }
    double term = alpha[i][n_mrk-1][0];
    for(int j=1; j < alpha[i][n_mrk-1].size(); j++)
      term = addlog(term, alpha[i][n_mrk-1][j]);
    loglike += term;
  }
  return(loglike);
}

// [[Rcpp::export]]
List est_hmm_map_biallelic_log_implementation(List PH,
                                              IntegerMatrix G,
                                              NumericMatrix pedigree,
                                              NumericVector rf,
                                              double err,
                                              bool verbose,
                                              bool detailed_verbose,
                                              double tol,
                                              bool ret_H0) {
  NumericVector ploidy_p1 = pedigree(_,2)/2 - 1;
  NumericVector ploidy_p2 = pedigree(_,3)/2 - 1;
  int n_ind = pedigree.nrow();
  NumericMatrix temp_phase_mat = PH[0];
  int n_mrk = temp_phase_mat.nrow();
  std::vector<int> pl{2,4,6};
  int k, k1,  maxit = 1000, flag=0;
  double s, loglike = 0.0, nr=0.0;

  // Getting the maximum ploidy level among founders
  int mpi1 = max(ploidy_p1);
  int mpi2 = max(ploidy_p2);
  int max_ploidy_id;
  if(mpi1 >= mpi2)
    max_ploidy_id = mpi1;
  else
    max_ploidy_id = mpi2;
  std::vector<double> rf_cur(rf.size());

  // HMM states that should be visited given the phase of
  // the founders, genotype of the offspring and pedigree
  List result = vs_biallelic(PH, G, pedigree);
  List ve = hmm_vectors(result);
  std::vector<std::vector<std::vector<int> > > v = ve["v"];
  std::vector<std::vector<std::vector<double> > > log_emit = ve["e"];

  //Initializing alpha and beta
  std::vector<std::vector<std::vector<double> > > alpha(n_ind);
  std::vector<std::vector<std::vector<double> > > beta(n_ind);
  for(int ind=0; ind < n_ind; ind++)
  {
    for(int i=0; i < n_mrk; i++)
    {
      std::vector<double> temp3(v[i][ind].size()/2);
      alpha[ind].push_back(temp3);
      beta[ind].push_back(temp3);
    }
  }
  //Initializing recombination number matrix
  std::vector< std::vector< std::vector<double> > > R;
  for(int j=0; j <= max_ploidy_id; j++)
    R.push_back(rec_num(pl[j]));

  for(int it=0; it < maxit ; it++)
  {
    //Initializing recombination fraction vector for Baum-Welch
    for(int j=0; j<n_mrk-1; j++)
    {
      rf_cur[j] = rf[j];
      rf[j] = 0.0;
    }
    //Initializing transition matrices
    std::vector< std::vector< std::vector< std::vector<double> > > >T;
    for(int j=0; j <= max_ploidy_id; j++)
    {
      std::vector< std::vector< std::vector<double> > > Ttemp;
      for(int i=0; i < n_mrk-1; i++)
      {
        Ttemp.push_back(log_transition(pl[j], rf_cur[i]));
      }
      T.push_back(Ttemp);
    }
    //Loop over all individuals
    for(int ind=0; ind < n_ind; ind++)
    {
      R_CheckUserInterrupt();
      for(int j=0; (unsigned)j < log_emit[0][ind].size(); j++)
        alpha[ind][0][j] = log_emit[0][ind][j];
      //forward
      for(k=1; k < n_mrk; k++)
      {
        int ngen_k0 = v[k-1][ind].size()/2;
        int ngen_k1 = v[k][ind].size()/2;
        for(int s1 = 0; s1 < ngen_k1; s1++)
        {
          alpha[ind][k][s1] = alpha[ind][k-1][0] +
            T[ploidy_p1[ind]][k-1][v[k-1][ind][0]][v[k][ind][s1]] +
            T[ploidy_p2[ind]][k-1][v[k-1][ind][ngen_k0]][v[k][ind][s1+ngen_k1]];
          for(int s0 = 1; s0 < ngen_k0; s0++){
            alpha[ind][k][s1] =  addlog(alpha[ind][k][s1], alpha[ind][k-1][s0] +
              T[ploidy_p1[ind]][k-1][v[k-1][ind][s0]][v[k][ind][s1]] +
              T[ploidy_p2[ind]][k-1][v[k-1][ind][s0+ngen_k0]][v[k][ind][s1+ngen_k1]]);
          }
          alpha[ind][k][s1] += log_emit[k][ind][s1];
        }
      }
      //backward
      std::fill(beta[ind][n_mrk-1].begin(), beta[ind][n_mrk-1].end(), 1.0);
      for(k1 = n_mrk-2; k1 >=0; k1--)
      {
        int ngen_k1 = v[k1][ind].size()/2;
        int ngen_k2 = v[k1+1][ind].size()/2;
        for(int s1 = 0; s1 < ngen_k1; s1++)
        {
          beta[ind][k1][s1] = beta[ind][k1+1][0] + T[ploidy_p1[ind]][k1][v[k1][ind][s1]][v[k1+1][ind][0]] +
            T[ploidy_p2[ind]][k1][v[k1][ind][s1+ngen_k1]][v[k1+1][ind][ngen_k2]] +
            log_emit[k1+1][ind][0];
          for(int s2 = 1; s2 < ngen_k2; s2++){
            beta[ind][k1][s1] = addlog(beta[ind][k1][s1], beta[ind][k1+1][s2] +
              T[ploidy_p1[ind]][k1][v[k1][ind][s1]][v[k1+1][ind][s2]] +
              T[ploidy_p2[ind]][k1][v[k1][ind][s1+ngen_k1]][v[k1+1][ind][s2+ngen_k2]] +
              log_emit[k1+1][ind][s2]);
          }
        }
      }
      if(ret_H0 == 0)
      {
        //Updating recombination fraction
        for(int k = 0; k < n_mrk-1; k++)
        {
          int ngeni = alpha[ind][k].size();
          int ngenj = beta[ind][k+1].size();
          vector<vector<double> > gamma(ngeni, vector<double>(ngenj));
          s=0.0;
          for(int i = 0; i < ngeni; i++)
          {
            for(int j = 0; j < ngenj; j++)
            {
              gamma[i][j] = alpha[ind][k][i] + beta[ind][k+1][j] +
                T[ploidy_p1[ind]][k][v[k][ind][i]][v[k+1][ind][j]] +
                T[ploidy_p2[ind]][k][v[k][ind][i+ngeni]][v[k+1][ind][j+ngenj]] +
                log_emit[k+1][ind][j];
              if(i==0 && j==0) s = gamma[i][j];
              else s = addlog(s, gamma[i][j]);
            }
          }
          for(int i=0; i < ngeni; i++)
          {
            for(int j=0; j < ngenj; j++)
            {
              nr=R[ploidy_p1[ind]][v[k][ind][i]][v[k+1][ind][j]] +
                R[ploidy_p2[ind]][v[k][ind][i+ngeni]][v[k+1][ind][j+ngenj]];
              if(s != 0)
                rf[k] +=  nr * exp(gamma[i][j] - s);
            }
          }
        }
      }
    } // loop over individuals

    //Likelihood using a specific recombination fraction vector
    //Usually, this is used to compute LOD Score under H0: rf=0.5
    if(ret_H0 == 1)
    {
      loglike = calc_loglike_log(v, log_emit, rf_cur, ploidy_p1, ploidy_p2);
      if(verbose)
        Rcpp::Rcout << "   \n";
      List z = List::create(wrap(loglike), rf_cur);
      return(z);
    }
    // re-scale
    for(int j=0; j<n_mrk-1; j++)
    {
      rf[j] /= (double)n_ind;
      if(rf[j] < tol/100.0) rf[j] = tol/100.0;
      else if(rf[j] > 0.5-tol/100.0) rf[j] = 0.5-tol/100.0;
    }
    // check convergence
    flag=0;
    for(int j=0; j < n_mrk-1; j++)
    {
      if(fabs(rf[j] - rf_cur[j]) > tol*(rf_cur[j]+tol*100.0))
      {
        flag = 1;
        break;
      }
    }
    if(verbose & !detailed_verbose){
      if(it%30 == 0 && it!=0) Rcout << "\n            ";
      Rcout << "." ;
    }
    if(detailed_verbose)
    {
      Rcpp::Rcout << "\t\n Iter: " << it+1 << "\t";
      for(int j = 0; j < n_mrk-1; j++)
      {
        Rcpp::Rcout.precision(3);
        Rcpp::Rcout << std::fixed << rf[j] << " ";
      }
    }
    if(!flag) break;
  }//end of EM algorithm
  if(flag && verbose) Rcpp::Rcout << "Didn't converge!\n";

  for(int j=0; j<n_mrk-1; j++)
    rf_cur[j] = rf[j];

  //Loglike computation
  loglike = calc_loglike_log(v, log_emit, rf_cur, ploidy_p1, ploidy_p2);
  if(verbose)
    Rcpp::Rcout << "\n";
  List z = List::create(wrap(loglike), wrap(rf_cur));
  return(z);
}








