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
 File: twopt_phasing.cpp

 Description:

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
#include "utils.h"
#include <math.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <unordered_set>
#include <progress.hpp>
#define THRESH 200.0
using namespace Rcpp;
using namespace std;

// Helper function for list concatenation
List cLists(List x, List y) {
  int nsize = x.size();
  int msize = y.size();
  List out(nsize + msize);

  for(int i = 0; i < nsize; i++) {
    out[i] = x[i];
  }
  for(int i = 0; i < msize; i++) {
    out[nsize+i] = y[i];
  }
  return(out);
}

// Helper function for matrix binding
NumericMatrix rbind_cpp(NumericMatrix a, NumericVector b){
  int nra = a.nrow(), nca = a.ncol();
  NumericMatrix out(nra+1, nca);

  // Fill the matrix
  for(int i=0; i < nra; i++)
    out(i, _) = a(i, _);
  out(nra, _) = b;

  return out;
}

// Helper function to initiate phasing matrix
NumericMatrix init_phase_mat(int d, int p){
  NumericMatrix v(1, p);

  if(d == 0)
    return v;
  else {
    for(int i = 0; i < d; i++)
      v(0,i) = 1.0;
    return v;
  }
}

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::export]]
List twopt_phasing_cpp(CharacterVector mrk_id,
                       int ploidy,
                       IntegerVector dose_vec,
                       NumericMatrix S,
                       int max_conf_number,
                       bool verbose){

  IntegerVector idx = 1;
  List H(1);
  H[0] = init_phase_mat(dose_vec[0],ploidy);

  int n = mrk_id.size();
  Progress p(n, verbose); //Progress bar

  for(int i = 1; i < mrk_id.size(); ++i){

    if (Progress::check_abort())
      return -1.0;

    int x = dose_vec[i];

    if(x == ploidy){
      for(int j = 0; j < H.size(); ++j)
        H[j] = rbind_cpp(H[j], NumericVector(ploidy, 1.0));
      idx.push_back(i);
    }

    else if(x == 0){
      for(int j = 0; j < H.size(); ++j)
        H[j] = rbind_cpp(H[j], NumericVector(ploidy, 0.0));
      idx.push_back(i);
    }

    else {
      NumericVector d(idx.size());
      for(int j = 0; j < idx.size(); ++j){
        d[j] = S(i, idx(j));
      }
      List vtemp(H.size());
      List Hres;
      //Rcpp::Rcout << "mrk: " << i << " --> Configs: " << H.size() << "\n";
      for(int j = 0; j < H.size(); ++j){
        vtemp[j] = find_valid_permutations(H[j], d, x);
        NumericMatrix z = vtemp[j];
        if(z.nrow() == 0) continue;

        List Htemp(z.nrow());

        for(int k = 0; k < z.nrow(); ++k){
          Htemp[k] = rbind_cpp(H[j], z(k, _));
        }
        Hres = cLists(Hres, filter_matrices(Htemp));
      }
      if(Hres.size() > 1)
        Hres = filter_matrices(Hres);

      if(Hres.size() > max_conf_number || Hres.size() == 0){
        //Rcpp::Rcout << "Skiping marker " << i << "\n";
        continue;
      }
      H = Hres;
      idx.push_back(i);
    }
    p.increment();
  }
  CharacterVector mrk_names(idx.size());
  for(int i  = 0; i < idx.size(); i++)
    mrk_names[i] = mrk_id[idx[i]];
  List z = List::create(Named("marker_names") = mrk_names,
                        Named("phase_configs") = H);
  return z;
}

// [[Rcpp::export]]
List phasing_one(CharacterVector mrk_id,
                 IntegerVector dose_vec,
                 NumericMatrix S,
                 NumericMatrix InitPh,
                 bool verbose){

  int ploidy = InitPh.ncol();//ploidy
  int n = mrk_id.size();//markers to be positioned
  List H(n);
  Progress p(n, verbose); //Progress bar
  for(int i = 0; i < n; ++i){
    if (Progress::check_abort())
      return -1.0;
    int x = dose_vec[i];
    if(x == ploidy)
      H[i] = make_mat(1.0, 1, ploidy);
    else if(x == 0)
      H[i] = make_mat(0.0, 1, ploidy);
    else {
      NumericMatrix z = find_valid_permutations(InitPh, S(i,_), x);
      if(z.nrow() == 0)
        H[i] = make_mat(-1.0, 1, ploidy); //FIXME
      else{
        List Htemp(z.nrow());
        for(int k = 0; k < z.nrow(); ++k)
          Htemp[k] = rbind_cpp(InitPh, z(k, _));
        List Hres = filter_matrices(Htemp);
        NumericMatrix w(Hres.size(), ploidy);
        for(int k = 0; k < Hres.size(); k++){
          NumericMatrix y = Hres[k];
          w(k,_) = y(y.nrow() - 1,_);
        }
        H[i] = w;
      }
    }
  }
  return H;
}




