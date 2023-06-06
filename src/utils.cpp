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
 File: utils.cpp

 Description:

 Functions Written by Marcelo Mollinari.

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version:       2022
 Last update: May 24, 2023
 */

#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <string_view>
#include "combinatorics.h"
#include "hmm_elements.h"

using namespace std;
using namespace Rcpp;

// Helper function to calculate the combinations
IntegerMatrix combn(NumericVector x, int m) {
  int n = x.size();
  if (m > n) return IntegerMatrix(0);
  if (m == n) return IntegerMatrix(x);

  int count = Rf_choose(n, m);
  IntegerMatrix result(m, count);

  std::vector<int> indices(m);
  std::iota(indices.begin(), indices.end(), 0);

  for (int col = 0; col < count; ++col) {
    for (int row = 0; row < m; ++row) {
      result(row, col) = x[indices[row]];
    }

    for (int i = m - 1; i >= 0; --i) {
      if (++indices[i] <= n - m + i) {
        while (++i < m) {
          indices[i] = indices[i - 1] + 1;
        }
        break;
      }
    }
  }

  return result;
}

// Function for rep_each
IntegerVector rep_each(IntegerVector x, int n) {
  int sz = x.size() * n;
  IntegerVector res(sz);
  for (int i = 0; i < x.size(); ++i) {
    std::fill_n(res.begin() + i * n, n, x[i]);
  }
  return res;
}

// Function for rep_len
IntegerVector rep_len(IntegerVector x, int n) {
  IntegerVector res(n);
  std::copy_n(x.begin(), std::min<int>(n, x.size()), res.begin());
  for (int i = x.size(); i < n; i += x.size()) {
    std::copy_n(x.begin(), std::min<int>(n - i, x.size()), res.begin() + i);
  }
  return res;
}

// Helper function for expand_grid
IntegerMatrix expand_grid(IntegerVector v1, IntegerVector v2) {
  int n_v1 = v1.size();
  int n_v2 = v2.size();
  IntegerVector res_v1 = rep_len(v1, n_v1 * n_v2); // Adjust the rep_len() usage
  IntegerVector res_v2 = rep_each(v2, n_v1); // Adjust the rep_each() usage
  return cbind(res_v2, res_v1);
}

// Helper function for which
IntegerVector which(LogicalVector x) {
  IntegerVector idx = seq_len(x.size());
  IntegerVector out;
  for (int i = 0; i < x.size(); ++i) {
    if (x[i]) {
      out.push_back(idx[i]);
    }
  }
  return out;
}

// Helper function for concatenating vectors
IntegerVector concatenate_vectors(IntegerVector x1, IntegerVector x2) {
  int n1 = x1.size();
  int n2 = x2.size();
  IntegerVector result(n1 + n2);
  std::copy(x1.begin(), x1.end(), result.begin());
  std::copy(x2.begin(), x2.end(), result.begin() + n1);
  return result;
}

// Helper function for calculate L and initialize H
List calculate_L_and_initialize_H(int n_fullsib_pop,
                                  int n_mrk,
                                  int n_ind,
                                  NumericVector ploidy_p1,
                                  NumericVector ploidy_p2) {
  List L(n_fullsib_pop);
  for (int i = 0; i < n_fullsib_pop; i++) {
    int ngam1 = R::choose(ploidy_p1[i], ploidy_p1[i] / 2);
    int ngam2 = R::choose(ploidy_p2[i], ploidy_p2[i] / 2);
    IntegerMatrix S = expand_grid(seq_len(ngam2) - 1, seq_len(ngam1) - 1);
    L[i] = S;
  }
  List H(n_mrk); //H[[n_mrk]][[n_ind]]
  List E(n_mrk);
  for (int k = 0; k < n_mrk; k++) {
    H[k] = List(n_ind);
    E[k] = List(n_ind);
  }
  return List::create(Named("L") = L,
                      Named("H") = H,
                      Named("E") = E);
}

// Helper function for unique rows
NumericMatrix retainUniqueAndSortByLastColumn(NumericMatrix x) {
  std::set<std::vector<double>> uniqueRows;
  int numRows = x.nrow();
  int numCols = x.ncol();

  // Iterate over each row
  for (int i = 0; i < numRows; i++) {
    std::vector<double> row;

    // Extract values from the current row
    for (int j = 0; j < numCols; j++) {
      row.push_back(x(i, j));
    }

    // Check if the row is already present
    if (uniqueRows.find(row) == uniqueRows.end()) {
      uniqueRows.insert(row); // Insert the row into the set
    }
  }

  // Sort the unique rows based on the last column value
  std::vector<std::vector<double>> sortedRows(uniqueRows.begin(), uniqueRows.end());
  std::sort(sortedRows.begin(), sortedRows.end(), [numCols](const std::vector<double>& a, const std::vector<double>& b) {
    return a[numCols - 1] < b[numCols - 1];
  });

  // Create a new matrix to store the sorted and unique rows
  NumericMatrix result(sortedRows.size(), numCols);
  for (int i = 0; i < sortedRows.size(); i++) {
    const std::vector<double>& row = sortedRows[i];
    for (int j = 0; j < numCols; j++) {
      result(i, j) = row[j];
    }
  }

  return result;
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

// Helper function for states to visit - biallelic
List vs_biallelic(List PH,
                  IntegerMatrix G,
                  NumericMatrix pedigree) {
  NumericMatrix unique_pop_mat = retainUniqueAndSortByLastColumn(pedigree);
  //Rcpp::Rcout << unique_pop_mat << "\n";
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
      //TEST
      //std::fill(temp_emit.begin(), temp_emit.end(), 1.0);
      std::fill(temp_emit.begin(), temp_emit.end(), 1.0/temp_emit.size());
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
          if (R_IsNA(G(k, ind_id[j]))) {  // ************************************* Ind NA
            H_k[ind_id[j]] = L[pop_id];
            E_k[ind_id[j]] = temp_emit;
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
            //TEST
            //std::fill(temp_emit2.begin(), temp_emit2.end(), 1.0);
            std::fill(temp_emit2.begin(), temp_emit2.end(), 1.0/temp_emit2.size());
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
        if (R_IsNA(G(k, j))) {  // ************************************* Ind NA
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
  List E = result["E"]; // Emission probabilities to be implemented
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
          if (R_IsNA(G(k, ind_id[j]))) {  // ************************************* Ind NA
            H_k[ind_id[j]] = L[pop_id];
            E_k[ind_id[j]] = temp_emit;
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
        if (R_IsNA(G(k, j))) {  // ************************************* Ind NA
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

// [[Rcpp::export]]
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

