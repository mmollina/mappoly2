/*
 MAPpoly: a package to construct genetic maps in autopolyploids
 Copyright (C) 2014-2022 Marcelo Mollinari

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
 File: hmm_elements.cpp

 Description:

 Functions Written by Marcelo Mollinari.

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version:       2022
 Last update: Oct 06, 2022
 */


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

/* FUNCTION: prob_k1_given_k_l_m
 This is equation 5 on the paper
 -----------------------------------------------------
 Calculates the genotypic transition probability based
 on l, which denotes the number of recombinant bivalents
 between loci k and k + 1 in one parent.
 */
double prob_k1_given_k_l_m(int ploidy, int l, double rf)
{
  return ((pow((1-rf),(ploidy/2-l))*pow(rf,l))/nChoosek(ploidy/2, l));
}

/* FUNCTION: log_prob_k1_given_k_l_m
 * -----------------------------------------------------
 */
double log_prob_k1_given_k_l_m(int ploidy, int l, double rf)
{
  return ((ploidy/2-l) * log(1-rf) + l * log(rf) - log(nChoosek(ploidy/2, l)));
}

/* FUNCTION: rec_num
 -----------------------------------------------------
 Returns a matrix containing the number of recombination
 events between loci k and k + 1 in one parent, given the
 ploidy level.
 */
std::vector<std::vector<double> > rec_num(int ploidy)
{
  int g = nChoosek(ploidy, ploidy/2);
  std::vector<std::vector<double> > R(g);
  for(int i = 0; (unsigned)i < R.size(); ++i)
  {
    for(int j = 0; j < g; ++j)
    {
      R[i].push_back(n_rec_given_genk_and_k1(ploidy, i+1, j+1)/(double)ploidy);
    }
  }
  return(R);
}

/* FUNCTION: transition
 -----------------------------------------------------
 Returns a transition matrix between loci k and k + 1 in
 one parent i.e. Prop(p_{k+1}|p_k), given the ploidy level
 and the recombination fraction.
 */
std::vector<std::vector<double> > transition(int ploidy, double rf)
{
  int g = nChoosek(ploidy, ploidy/2);
  std::vector<std::vector<double> > T(g);
  for(int i = 0; i < g; ++i)
  {
    for(int j = 0; j < g; ++j)
    {
      T[i].push_back(prob_k1_given_k_l_m(ploidy,
                                         n_rec_given_genk_and_k1(ploidy, i+1, j+1),
                                         rf));
    }
  }
  return(T);
}

/* FUNCTION: transition
 -----------------------------------------------------
 Returns a transition matrix between loci k and k + 1 in
 one parent i.e. Prop(p_{k+1}|p_k), given the ploidy level
 and the recombination fraction.
 */
std::vector<std::vector<double> > log_transition(int ploidy, double rf)
{
  int g = nChoosek(ploidy, ploidy/2);
  std::vector<std::vector<double> > T(g);
  for(int i = 0; i < g; ++i)
  {
    for(int j = 0; j < g; ++j)
    {
      T[i].push_back(log_prob_k1_given_k_l_m(ploidy,
                                             n_rec_given_genk_and_k1(ploidy, i+1, j+1),
                                             rf));
    }
  }
  return(T);
}

/* FUNCTION: forward_emit (with both informative parents)
 -----------------------------------------------------
 Classical forward equation presented in Rabiner 1989.
 */
std::vector<double> forward_emit(std::vector<double>& fk,
                                 std::vector<int>& ik,
                                 std::vector<int>& ik1,
                                 std::vector<double>& emit,
                                 std::vector<std::vector<double> >& T1,
                                 std::vector<std::vector<double> >& T2)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<double> fk1(ngenk1);
  std::fill(fk1.begin(), fk1.end(), 0.0);
  for(int k1 = 0; k1 < ngenk1; k1++ )
  {
    for(int k = 0; k < ngenk; k++ )
    {
      fk1[k1] = fk1[k1] + fk[k] * T1[ik[k]][ik1[k1]] * T2[ik[k+ngenk]][ik1[k1+ngenk1]];
    }
    fk1[k1] = fk1[k1] * emit[k1];
  }
  return(fk1);
}

/* FUNCTION: forward_emit (with both informative parents)
 -----------------------------------------------------
 Classical forward equation presented in Rabiner 1989.
 */
std::vector<long double> forward_emit_highprec(std::vector<long double>& fk,
                                               std::vector<int>& ik,
                                               std::vector<int>& ik1,
                                               std::vector<double>& emit,
                                               std::vector<std::vector<double> >& T1,
                                               std::vector<std::vector<double> >& T2)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<long double> fk1(ngenk1);
  std::fill(fk1.begin(), fk1.end(), 0.0);
  for(int k1 = 0; k1 < ngenk1; k1++ )
  {
    for(int k = 0; k < ngenk; k++ )
    {
      fk1[k1] = fk1[k1] + fk[k] * T1[ik[k]][ik1[k1]] * T2[ik[k+ngenk]][ik1[k1+ngenk1]];
    }
    fk1[k1] = fk1[k1] * emit[k1];
  }
  return(fk1);
}

/* FUNCTION: forward (with one informative parent)
 -----------------------------------------------------
 Classical forward equation presented in Rabiner 1989.
 */
std::vector<double> forward_emit_single_parent(int m,
                                               std::vector<double>& fk,
                                               std::vector<int>& ik,
                                               std::vector<int>& ik1,
                                               std::vector<double>& emit,
                                               std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size();
  int ngenk1 = ik1.size();
  std::vector<double> fk1(ngenk1);
  std::fill(fk1.begin(), fk1.end(), 0.0);
  for(int k1 = 0; k1 < ngenk1; k1++ )
  {
    for(int k = 0; k < ngenk; k++ )
    {
      fk1[k1] = fk1[k1] + fk[k] * T[ik[k]][ik1[k1]];
    }
    fk1[k1] = fk1[k1] * emit[k1];
  }
  return(fk1);
}


/* FUNCTION: forward_emit (with both informative parents)
 -----------------------------------------------------
 Classical forward equation presented in Rabiner 1989.
 */
std::vector<double> log_forward_emit(std::vector<double>& fk,
                                     std::vector<int>& ik,
                                     std::vector<int>& ik1,
                                     std::vector<double>& emit,
                                     std::vector<std::vector<double> >& T1,
                                     std::vector<std::vector<double> >& T2)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<double> fk1(ngenk1);
  std::fill(fk1.begin(), fk1.end(), 1.0);
  for(int k1 = 0; k1 < ngenk1; k1++ )
  {
    for(int k = 0; k < ngenk; k++ )
    {
      fk1[k1] = addlog(fk1[k1], fk[k] + T1[ik[k]][ik1[k1]] + T2[ik[k+ngenk]][ik1[k1+ngenk1]]);
    }
    fk1[k1] = fk1[k1] + emit[k1];
  }
  return(fk1);
}


/* FUNCTION: backward (with both informative parents)
 -----------------------------------------------------
 Classical backward equation presented in Rabiner 1989.
 */
std::vector<double> backward_emit(std::vector<double>& fk1,
                                  std::vector<int>& ik,
                                  std::vector<int>& ik1,
                                  std::vector<double>& emit,
                                  std::vector<std::vector<double> >& T1,
                                  std::vector<std::vector<double> >& T2)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<double> fk(ngenk);
  std::fill(fk.begin(), fk.end(), 0.0);
  for(int k = 0; k < ngenk; k++ )
  {
    for(int k1 = 0; k1 < ngenk1; k1++ )
    {
      fk[k] =  fk[k] + fk1[k1] * T1[ik[k]][ik1[k1]] * T2[ik[k+ngenk]][ik1[k1+ngenk1]] * emit[k1];
    }
  }
  return(fk);
}
/* FUNCTION: backward (with both informative parents)
 -----------------------------------------------------
 Classical backward equation presented in Rabiner 1989.
 */
std::vector<long double> backward_emit_highprec(std::vector<long double>& fk1,
                                                std::vector<int>& ik,
                                                std::vector<int>& ik1,
                                                std::vector<double>& emit,
                                                std::vector<std::vector<double> >& T1,
                                                std::vector<std::vector<double> >& T2)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<long double> fk(ngenk);
  std::fill(fk.begin(), fk.end(), 0.0);
  for(int k = 0; k < ngenk; k++ )
  {
    for(int k1 = 0; k1 < ngenk1; k1++ )
    {
      fk[k] =  fk[k] + fk1[k1] * T1[ik[k]][ik1[k1]] * T2[ik[k+ngenk]][ik1[k1+ngenk1]] * emit[k1];
    }
  }
  return(fk);
}

/* FUNCTION: backward (with one informative parent)
 -----------------------------------------------------
 Classical backward equation presented in Rabiner 1989.
 */
std::vector<double> backward_emit_single_parent(int m,
                                                std::vector<double>& fk1,
                                                std::vector<int>& ik,
                                                std::vector<int>& ik1,
                                                std::vector<double>& emit,
                                                std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size();
  int ngenk1 = ik1.size();
  std::vector<double> fk(ngenk);
  std::fill(fk.begin(), fk.end(), 0.0);
  for(int k = 0; k < ngenk; k++ )
  {
    for(int k1 = 0; k1 < ngenk1; k1++ )
    {
      fk[k] =  fk[k] + fk1[k1] * T[ik[k]][ik1[k1]] * emit[k1];
    }
  }
  return(fk);
}

/* FUNCTION: backward (with both informative parents)
 -----------------------------------------------------
 Classical backward equation presented in Rabiner 1989.
 */
std::vector<double> log_backward_emit(std::vector<double>& fk1,
                                      std::vector<int>& ik,
                                      std::vector<int>& ik1,
                                      std::vector<double>& emit,
                                      std::vector<std::vector<double> >& T1,
                                      std::vector<std::vector<double> >& T2)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<double> fk(ngenk);
  std::fill(fk.begin(), fk.end(), 0.0);
  for(int k = 0; k < ngenk; k++ )
  {
    for(int k1 = 0; k1 < ngenk1; k1++ )
    {
      fk[k] =  addlog(fk[k], fk1[k1] +
                             T1[ik[k]][ik1[k1]] +
                             T2[ik[k+ngenk]][ik1[k1+ngenk1]] +
                             emit[k1]);
    }
  }
  return(fk);
}

/* FUNCTION: addlog
 * -----------------------------------------------------
 */
double addlog(double a, double b)
{
  if(b > a + THRESH) return(b);
  else if(a > b + THRESH) return(a);
  else return(a + log1p(exp(b-a)));
}
