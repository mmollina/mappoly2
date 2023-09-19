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

 Functions Written by Marcelo Mollinari.

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version:       2022
 Last update: Jun 14, 2023
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
#define SparseThreshold 0.005
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

/* FUNCTION: log transition
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
 Forward equation presented in Rabiner 1989.
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
    fk1[k1] *= emit[k1];
  }
  return(fk1);
}

/* FUNCTION: forward (with one informative parent)
 -----------------------------------------------------
 Forward equation presented in Rabiner 1989.
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

/* FUNCTION: addlog
 * -----------------------------------------------------
 */
double addlog(double a, double b)
{
  if(b > a + THRESH) return(b);
  else if(a > b + THRESH) return(a);
  else return(a + log1p(exp(b-a)));
}


// Helper function for states to visit - multiallelic
List vs_multiallelic_Rcpp(List PH,
                          List GENO,
                          NumericMatrix pedigree) {
  NumericMatrix unique_pop_mat = retainUniqueAndSortByLastColumn(pedigree);
  int n_fullsib_pop = unique_pop_mat.nrow();
  int n_ind = pedigree.nrow();
  NumericMatrix temp_phase_mat = PH[0];
  int n_mrk = temp_phase_mat.nrow();
  NumericVector ploidy_p1 = unique_pop_mat(_,2);
  NumericVector ploidy_p2 = unique_pop_mat(_,3);
  List result = calculate_L_and_initialize_H(n_fullsib_pop, n_mrk, n_ind, ploidy_p1, ploidy_p2);
  List L = result["L"];
  List H = result["H"];
  List E = result["E"]; // Emission probabilities to be implemented
  for(int k = 0; k < n_mrk; k++) { // ***************************************** Markers
    List H_k(n_ind);
    List E_k(n_ind);// Emission probabilities to be implemented
    for(int pop_id = 0; pop_id < n_fullsib_pop; pop_id ++){ // ************ Pop id
      // Emission function
      NumericMatrix L_mat = as<NumericMatrix>(L[pop_id]);
      NumericMatrix temp_emit(L_mat.nrow(), 1);
      std::fill(temp_emit.begin(), temp_emit.end(), 1);
      // std::fill(temp_emit.begin(), temp_emit.end(), 1/temp_emit.size());
      // States to visit
      IntegerVector ind_id = which(pedigree(_,4) == pop_id + 1) - 1;
      NumericMatrix matrix_PH1 =  PH[unique_pop_mat(pop_id,0) - 1];
      NumericMatrix matrix_PH2 = PH[unique_pop_mat(pop_id,1) - 1];
      NumericVector v1 = matrix_PH1(k, _);
      NumericVector v2 = matrix_PH2(k, _);
      if (any(is_na(v1)).is_true() || any(is_na(v2)).is_true()) {
        for(int j = 0; j < ind_id.size(); j++) { // ***************************** Ind Parents NA
          H_k[ind_id[j]] = L[pop_id];
          E_k[ind_id[j]] = temp_emit;
        } // end individual loop when NA
      } else {
        IntegerMatrix x1 = combn(v1, v1.size() / 2);
        IntegerMatrix x2 = combn(v2, v2.size() / 2);
        IntegerMatrix x(x1.nrow() + x2.nrow(), x1.ncol() * x2.ncol());
        for (int i = 0; i < x1.ncol(); i++) {
          for (int j = 0; j < x2.ncol(); j++) {
            IntegerVector temp = concatenate_vectors(x1(_, i), x2(_, j));
            temp.sort();
            x(_, (i * x2.ncol()) + j) = temp;
          }
        }
        for(int j = 0; j < ind_id.size(); j++) { // ***************************** Ind ALL
          IntegerMatrix matrix_GENO = as<IntegerMatrix>(GENO[ind_id[j]]);
          IntegerVector a = matrix_GENO(k, _);
          if (is_true(any(is_na(a)))) {  // ************************************* Ind NA
            H_k[ind_id[j]] = L[pop_id];
            E_k[ind_id[j]] = temp_emit;
          } else{  // *********************************************************** Ind Visit
            a.sort();
            IntegerVector y;
            for (int i = 0; i < x.ncol(); i++) {
              bool all_equal = true;
              for (int l = 0; l < x.nrow(); l++) {
                if (x(l, i) != a[l]) {
                  all_equal = false;
                  break;
                }
              }
              if (all_equal) {
                y.push_back(i);
              }
            }
            IntegerMatrix L_pop = L[pop_id];
            IntegerMatrix subset_L_pop(y.size(), L_pop.ncol());
            for (int row = 0; row < y.size(); ++row) {
              subset_L_pop(row, _) = L_pop(y[row], _);
            }
            H_k[ind_id[j]] = subset_L_pop;
            NumericMatrix temp_emit2(y.size(), 1);
            std::fill(temp_emit2.begin(), temp_emit2.end(), 1);
            E_k[ind_id[j]] = temp_emit2;
          }
        }
      }
    } // end population loop
    H[k] = H_k;
    E[k] = E_k;
  } // end marker loop
  return List::create(Named("states") = H,
                      Named("emit") = E);
}

// Helper function for states to visit - marker in PH and G will be inserted between M(,1) and M(,2)
// [[Rcpp::export]]
List vs_inserted_mrk(List PH, //list of vectors
                     IntegerVector G, //vector of dosages for inserted marker
                     NumericMatrix pedigree,
                     NumericMatrix  M,
                     IntegerVector idx) {
  NumericMatrix unique_pop_mat = retainUniqueAndSortByLastColumn(pedigree);
  int n_fullsib_pop = unique_pop_mat.nrow();
  int n_ind = pedigree.nrow();
  int n_mrk = 3;
  NumericVector ploidy_p1 = unique_pop_mat(_,2);
  NumericVector ploidy_p2 = unique_pop_mat(_,3);
  List result = calculate_L_and_initialize_H(n_fullsib_pop, n_mrk, n_ind,
                                             ploidy_p1,ploidy_p2);
  List L = result["L"];
  List H = result["H"];
  List E = result["E"];
  List H_k(n_ind);
  List E_k(n_ind);
  // For loop to retrieve the states associated with markers that have not yet been mapped
  for(int pop_id = 0; pop_id < n_fullsib_pop; pop_id ++){ // ************ Pop id
    // Emission function
    NumericMatrix L_mat = as<NumericMatrix>(L[pop_id]);
    NumericMatrix temp_emit(L_mat.nrow(), 1);
    std::fill(temp_emit.begin(), temp_emit.end(), 1.0);
    // States to visit
    IntegerVector ind_id = which(pedigree(_,4) == pop_id + 1) - 1;
    NumericVector v1 = PH[unique_pop_mat(pop_id,0) - 1];
    NumericVector v2 = PH[unique_pop_mat(pop_id,1) - 1];
    if (any(is_na(v1)).is_true() || any(is_na(v2)).is_true()) {
      for(int j = 0; j < ind_id.size(); j++) { // ***************************** Ind Parents NA
        H_k[ind_id[j]] = L[pop_id];
        E_k[ind_id[j]] = temp_emit;
      } // end individual loop when NA
    } else {
      IntegerMatrix x1 = combn(v1, v1.size() / 2);
      IntegerMatrix x2 = combn(v2, v2.size() / 2);
      IntegerVector x(x1.ncol() * x2.ncol());
      for (int i = 0; i < x1.ncol(); i++) {
        for (int j = 0; j < x2.ncol(); j++) {
          x((i * x2.ncol()) + j) = sum(x1(_ ,i)) + sum(x2(_ ,j));
        }
      }
      for(int j = 0; j < ind_id.size(); j++) { // ***************************** Ind ALL
        if (G(ind_id[j]) < 0) {  // ************************************* Ind NA
          H_k[ind_id[j]] = L[pop_id];
          E_k[ind_id[j]] = temp_emit;
        } else{  // *********************************************************** Ind Visit
          IntegerVector y;
          for (int i = 0; i < x.size(); i++) {
            if (x(i) == G(ind_id[j])) {
              y.push_back(i);
            }
          }
          IntegerMatrix L_pop = L[pop_id];
          IntegerMatrix subset_L_pop(y.size(), L_pop.ncol());
          for (int row = 0; row < y.size(); ++row) {
            subset_L_pop(row, _) = L_pop(y[row], _);
          }
          H_k[ind_id[j]] = subset_L_pop;
          NumericMatrix temp_emit2(y.size(), 1);
          std::fill(temp_emit2.begin(), temp_emit2.end(), 1.0);
          E_k[ind_id[j]] = temp_emit2;
        }
      }
    }
  } // end population loop
  H[idx(1)] = H_k;
  E[idx(1)] = E_k;
  // For loop to obtain the states of markers that are positioned adjacently to the un-mapped marker
  NumericMatrix L_mat = as<NumericMatrix>(L[0]);   // FIXME for multi-population: zero here refers to the first (and only bi-parental population)
  for(int k = 0; k < 2; k++){
    List H_k(n_ind);
    List E_k(n_ind);
    int cte = 0;
    for(int i = 0; i < n_ind; i++){
      NumericVector v(ploidy_p1[0]+ploidy_p2[0]); // FIXME for multi-population
      for(int j = 0; j < v.size(); j++) {
        v[j] = M(j + cte, k);
      }
      cte = cte + ploidy_p1[0]+ploidy_p2[0];      // FIXME for multi-population
      NumericVector w = calculate_hmm_combinatorial_products(v,ploidy_p1[0],ploidy_p2[0]);
      IntegerVector y;
      for (int j = 0; j < w.size(); j++) {
        if(w[j] > SparseThreshold){
          y.push_back(j);
        }
      }
      IntegerMatrix subset_L_pop(y.size(), L_mat.ncol());
      NumericMatrix subset_E_pop(y.size(), 1);
      for (int row = 0; row < y.size(); ++row) {
        subset_L_pop(row, _) = L_mat(y[row], _);
        subset_E_pop(row, 0) = w(y[row]);
      }
      H_k[i] = subset_L_pop;
      E_k[i] = subset_E_pop;
    }
    H[idx(k * 2)] = H_k;
    E[idx(k * 2)] = E_k;
  }
  return List::create(Named("states") = H,
                      Named("emit") = E);
}

// Helper function for states to visit - biallelic
List vs_biallelic(List PH,
                  IntegerMatrix G,
                  NumericMatrix pedigree) {
  NumericMatrix unique_pop_mat = retainUniqueAndSortByLastColumn(pedigree);
  int n_fullsib_pop = unique_pop_mat.nrow();
  int n_ind = pedigree.nrow();
  NumericMatrix temp_phase_mat = PH[0];
  int n_mrk = temp_phase_mat.nrow();
  NumericVector ploidy_p1 = unique_pop_mat(_,2);
  NumericVector ploidy_p2 = unique_pop_mat(_,3);
  List result = calculate_L_and_initialize_H(n_fullsib_pop, n_mrk, n_ind, ploidy_p1, ploidy_p2);
  List L = result["L"];
  List H = result["H"];
  List E = result["E"]; // Emission probabilities to be implemented
  for(int k = 0; k < n_mrk; k++) { // ***************************************** Markers
    List H_k(n_ind);
    List E_k(n_ind);
    for(int pop_id = 0; pop_id < n_fullsib_pop; pop_id ++){ // ************ Pop id
      // Emission function
      NumericMatrix L_mat = as<NumericMatrix>(L[pop_id]);
      NumericMatrix temp_emit(L_mat.nrow(), 1);
      std::fill(temp_emit.begin(), temp_emit.end(), 1.0);
      // States to visit
      IntegerVector ind_id = which(pedigree(_,4) == pop_id + 1) - 1;
      NumericMatrix matrix_PH1 =  PH[unique_pop_mat(pop_id,0) - 1];
      NumericMatrix matrix_PH2 = PH[unique_pop_mat(pop_id,1) - 1];
      NumericVector v1 = matrix_PH1(k, _);
      NumericVector v2 = matrix_PH2(k, _);
      if (any(is_na(v1)).is_true() || any(is_na(v2)).is_true()) {
        for(int j = 0; j < ind_id.size(); j++) { // ***************************** Ind Parents NA
          H_k[ind_id[j]] = L[pop_id];
          E_k[ind_id[j]] = temp_emit;
        } // end individual loop when NA
      } else {
        IntegerMatrix x1 = combn(v1, v1.size() / 2);
        IntegerMatrix x2 = combn(v2, v2.size() / 2);
        IntegerVector x(x1.ncol() * x2.ncol());
        for (int i = 0; i < x1.ncol(); i++) {
          for (int j = 0; j < x2.ncol(); j++) {
            x((i * x2.ncol()) + j) = sum(x1(_ ,i)) + sum(x2(_ ,j));
          }
        }
        for(int j = 0; j < ind_id.size(); j++) { // ***************************** Ind ALL
          if (G(k, ind_id[j]) < 0) {  // ************************************* Ind NA
            H_k[ind_id[j]] = L[pop_id];
            E_k[ind_id[j]] = temp_emit;
            // Rcout << "Here!!!\n" << "\n";
          } else{  // *********************************************************** Ind Visit
            IntegerVector y;
            for (int i = 0; i < x.size(); i++) {
              if (x(i) == G(k, ind_id[j])) {
                y.push_back(i);
              }
            }
            IntegerMatrix L_pop = L[pop_id];
            IntegerMatrix subset_L_pop(y.size(), L_pop.ncol());
            for (int row = 0; row < y.size(); ++row) {
              subset_L_pop(row, _) = L_pop(y[row], _);
            }
            H_k[ind_id[j]] = subset_L_pop;
            NumericMatrix temp_emit2(y.size(), 1);
            std::fill(temp_emit2.begin(), temp_emit2.end(), 1.0);
            E_k[ind_id[j]] = temp_emit2;
          }
        }
      }
    } // end population loop
    H[k] = H_k;
    E[k] = E_k;
  } // end marker loop
  return List::create(Named("states") = H,
                      Named("emit") = E);
}

// Helper function for states to visit - biallelic - single parent
List vs_biallelic_single(NumericMatrix PH,
                         IntegerMatrix G) {
  int ploidy = PH.ncol();
  int n_ind = G.ncol();
  int n_mrk = G.nrow();
  List H(n_mrk); //H[[n_mrk]][[n_ind]]
  List E(n_mrk);
  for (int k = 0; k < n_mrk; k++) {
    H[k] = List(n_ind);
    E[k] = List(n_ind);// Emission probabilities to be implemented
  }
  for(int k = 0; k < n_mrk; k++) { // ***************************************** Markers
    List H_k(n_ind);
    List E_k(n_ind);
    int ngam = R::choose(ploidy, ploidy / 2);
    NumericMatrix L_mat(ngam, 1);
    NumericMatrix temp_emit(ngam, 1);
    for (int i = 0; i < ngam; i++) {
      L_mat(i, 0) = i;
      temp_emit(i,0) = 1.0;    // Emission [to be implemented]
    }
    // States to visit
    // "ind_id" will be 0:n_ind
    NumericVector v1 = PH(k, _);
    if (any(is_na(v1)).is_true()) {
      for(int j = 0; j < n_ind; j++) { // ***************************** Ind Parents NA
        H_k[j] = L_mat;
        E_k[j] = temp_emit;
      } // end individual loop when NA
    } else {
      IntegerMatrix x1 = combn(v1, v1.size()/2);
      IntegerVector x(x1.ncol());
      for (int i = 0; i < x1.ncol(); i++) {
        x(i) = sum(x1(_ ,i));
      }
      for(int j = 0; j < n_ind; j++) { // ***************************** Ind ALL
        if (G(k, j) < 0) {  // ************************************* Ind NA
          H_k[j] = L_mat;
          E_k[j] = temp_emit;
        } else{  // *********************************************************** Ind Visit
          IntegerVector y;
          for (int i = 0; i < x.size(); i++) {
            if (x(i) == G(k, j)) {
              y.push_back(i);
            }
          }
          IntegerMatrix subset_L_pop(y.size(), L_mat.ncol());
          for (int row = 0; row < y.size(); ++row) {
            subset_L_pop(row, _) = L_mat(y[row], _);
          }
          H_k[j] = subset_L_pop;
          NumericMatrix temp_emit2(y.size(), 1);
          std::fill(temp_emit2.begin(), temp_emit2.end(), 1.0);
          E_k[j] = temp_emit2;
        }
      }
    }
    H[k] = H_k;
    E[k] = E_k;
  } // end marker loop
  return List::create(Named("states") = H,
                      Named("emit") = E);
}

// Helper function for states to visit - biallelic
// using emission to model error
// [[Rcpp::export]]
List vs_biallelic_error(List PH,
                        IntegerMatrix G,
                        NumericMatrix pedigree,
                        double err,
                        bool logatithm) {
  if(err<1e-50) err = 1e-50;
  NumericMatrix unique_pop_mat = retainUniqueAndSortByLastColumn(pedigree);
  int n_fullsib_pop = unique_pop_mat.nrow();
  int n_ind = pedigree.nrow();
  NumericMatrix temp_phase_mat = PH[0];
  int n_mrk = temp_phase_mat.nrow();
  NumericVector ploidy_p1 = unique_pop_mat(_,2);
  NumericVector ploidy_p2 = unique_pop_mat(_,3);
  List result = calculate_L_and_initialize_H(n_fullsib_pop, n_mrk, n_ind, ploidy_p1, ploidy_p2);
  List L = result["L"];
  List H = result["H"];
  List E = result["E"];
  for(int k = 0; k < n_mrk; k++) { // ***************************************** Markers
    List H_k(n_ind);
    List E_k(n_ind);
    for(int pop_id = 0; pop_id < n_fullsib_pop; pop_id ++){ // ************ Pop id
      // Emission function
      NumericMatrix L_mat = as<NumericMatrix>(L[pop_id]);
      NumericMatrix temp_emit(L_mat.nrow(), 1);
      if(logatithm)
        std::fill(temp_emit.begin(), temp_emit.end(), -log(temp_emit.size()));
      else
        std::fill(temp_emit.begin(), temp_emit.end(), (1.0/temp_emit.size()));
      // States to visit
      IntegerVector ind_id = which(pedigree(_,4) == pop_id + 1) - 1;
      NumericMatrix matrix_PH1 =  PH[unique_pop_mat(pop_id,0) - 1];
      NumericMatrix matrix_PH2 = PH[unique_pop_mat(pop_id,1) - 1];
      NumericVector v1 = matrix_PH1(k, _);
      NumericVector v2 = matrix_PH2(k, _);
      if (any(is_na(v1)).is_true() || any(is_na(v2)).is_true()) {
        for(int j = 0; j < ind_id.size(); j++) { // ***************************** Ind Parents NA
          H_k[ind_id[j]] = L[pop_id];
          E_k[ind_id[j]] = temp_emit;
        } // end individual loop when NA
      } else {
        IntegerMatrix x1 = combn(v1, v1.size() / 2);
        IntegerMatrix x2 = combn(v2, v2.size() / 2);
        IntegerVector x(x1.ncol() * x2.ncol());
        for (int i = 0; i < x1.ncol(); i++) {
          for (int j = 0; j < x2.ncol(); j++) {
            x((i * x2.ncol()) + j) = sum(x1(_ ,i)) + sum(x2(_ ,j));
          }
        }
        for(int j = 0; j < ind_id.size(); j++) { // ***************************** Ind ALL
          if (G(k, ind_id[j]) < 0) {  // ************************************* Ind NA
            H_k[ind_id[j]] = L[pop_id];
            E_k[ind_id[j]] = temp_emit;
            //Rcout << "Here!!!\n" << "\n";
          } else{  // *********************************************************** Ind Visit
            IntegerVector y;
            for (int i = 0; i < x.size(); i++) {
              if (x(i) == G(k, ind_id[j])) {
                y.push_back(i);
              }
            }
            H_k[ind_id[j]] = L[pop_id];
            NumericMatrix temp_emit2(L_mat.nrow(), 1);
            if(logatithm)
              std::fill(temp_emit2.begin(), temp_emit2.end(), log(err) - log((x.size() - y.size())));
            else
              std::fill(temp_emit2.begin(), temp_emit2.end(), err/(x.size() - y.size()));

            for (int row = 0; row < y.size(); ++row){
              if(logatithm)
                temp_emit2(y[row], 0) = log(1.0 - err)-log(y.size());
              else
                temp_emit2(y[row], 0) = (1.0 - err)/y.size();
            }
            E_k[ind_id[j]] = temp_emit2;
          }
        }
      }
    } // end population loop
    H[k] = H_k;
    E[k] = E_k;
  } // end marker loop
  return List::create(Named("states") = H,
                      Named("emit") = E);
}

// Helper function for states to visit - biallelic - single parent
// using emission to model error
List vs_biallelic_error_single(NumericMatrix PH,
                               IntegerMatrix G,
                               double err,
                               bool logatithm) {
  int ploidy = PH.ncol();
  int n_ind = G.ncol();
  int n_mrk = G.nrow();
  List H(n_mrk); //H[[n_mrk]][[n_ind]]
  List E(n_mrk);
  for (int k = 0; k < n_mrk; k++) {
    H[k] = List(n_ind);
    E[k] = List(n_ind);// Emission probabilities to be implemented
  }
  for(int k = 0; k < n_mrk; k++) { // ***************************************** Markers
    List H_k(n_ind);
    List E_k(n_ind);
    int ngam = R::choose(ploidy, ploidy / 2);
    NumericMatrix L_mat(ngam, 1);
    NumericMatrix temp_emit(ngam, 1);
    for (int i = 0; i < ngam; i++) {
      L_mat(i, 0) = i;
      if(logatithm)
        temp_emit(i,0) = -log(ngam);
      else
        temp_emit(i,0) = 1.0/ngam;
    }
    // States to visit
    // "ind_id" will be 0:n_ind
    NumericVector v1 = PH(k, _);
    if (any(is_na(v1)).is_true()) {
      for(int j = 0; j < n_ind; j++) { // ***************************** Ind Parents NA
        H_k[j] = L_mat;
        E_k[j] = temp_emit;
      } // end individual loop when NA
    } else {
      IntegerMatrix x1 = combn(v1, v1.size()/2);
      IntegerVector x(x1.ncol());
      for (int i = 0; i < x1.ncol(); i++) {
        x(i) = sum(x1(_ ,i));
      }
      for(int j = 0; j < n_ind; j++) { // ***************************** Ind ALL
        if (G(k, j) < 0) {  // ************************************* Ind NA
          H_k[j] = L_mat;
          E_k[j] = temp_emit;
        } else{  // *********************************************************** Ind Visit
          IntegerVector y;
          for (int i = 0; i < x.size(); i++) {
            if (x(i) == G(k, j)) {
              y.push_back(i);
            }
          }
          IntegerMatrix subset_L_pop(y.size(), L_mat.ncol());
          H_k[j] = L_mat;
          NumericMatrix temp_emit2(L_mat.nrow(), 1);
          if(logatithm)
            std::fill(temp_emit2.begin(), temp_emit2.end(), log(err/(x.size() - y.size())));
          else
            std::fill(temp_emit2.begin(), temp_emit2.end(), (err/(x.size() - y.size())));

          for (int row = 0; row < y.size(); ++row){
            if(logatithm)
              temp_emit2(y[row], 0) = log((1.0 - err)/y.size());
            else
              temp_emit2(y[row], 0) = ((1.0 - err)/y.size());
          }
          E_k[j] = temp_emit2;
        }
      }
    }
    H[k] = H_k;
    E[k] = E_k;
  } // end marker loop
  return List::create(Named("states") = H,
                      Named("emit") = E);
}

// Wrapper function - states to visit - both parents
List visit_states_biallelic(List PH,
                            IntegerMatrix G,
                            NumericMatrix pedigree,
                            double err){
  if(err <= 0.001){
    return(vs_biallelic(PH, G, pedigree));
  } else {
    return(vs_biallelic_error(PH, G, pedigree, err, 0));
  }
}

// Wrapper function - states to visit - single parents
List visit_states_biallelic_single(NumericMatrix PH,
                                   IntegerMatrix G,
                                   double err){
  if(err <= 0.001){
    return(vs_biallelic_single(PH, G));
  } else {
    return(vs_biallelic_error_single(PH, G, err, 0));
  }
}

List hmm_vectors(List input_list) {
  // Getting the input lists
  List haplo = input_list["states"];
  List emit = input_list["emit"];

  // Initializing v: states hmm should visit for each marker
  // Initializing e: emission probabilities associated to the states hmm should visit for each marker
  std::vector<std::vector<std::vector<int> > > v;
  std::vector<std::vector<std::vector<double> > > e;
  for(int i=0; i < haplo.size(); i++) // i: number of markers
  {
    Rcpp::List haplo_temp(haplo[i]); //states hmm should visit for marker i
    Rcpp::List emit_temp(emit[i]); //emission probs. for states hmm should visit for marker i
    std::vector<std::vector<int> > v1;
    std::vector<std::vector<double> > e1;
    for(int j=0; j < haplo_temp.size(); j++) //iterate for all j individuals
    {
      Rcpp::NumericMatrix M_temp = haplo_temp[j];
      Rcpp::NumericVector E_temp = emit_temp[j];
      std::vector<int> v2 = as<std::vector<int> >(M_temp);
      std::vector<double> e2 = as<std::vector<double> >(E_temp);
      v1.push_back(v2);
      e1.push_back(e2);
    }
    v.push_back(v1);
    e.push_back(e1);
  }
  return List::create(Named("v") = v, Named("e") = e);
}
