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
 Last update: May 24, 2023
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
                                  NumericVector pl1,
                                  NumericVector pl2) {
  List L(n_fullsib_pop);
  for (int i = 0; i < n_fullsib_pop; i++) {
    int ngam1 = R::choose(pl1[i], pl1[i] / 2);
    int ngam2 = R::choose(pl2[i], pl2[i] / 2);
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
  NumericMatrix temp_mat = PH[0];
  int n_mrk = temp_mat.nrow();
  NumericVector pl1 = unique_pop_mat(_,2);
  NumericVector pl2 = unique_pop_mat(_,3);
  List result = calculate_L_and_initialize_H(n_fullsib_pop, n_mrk, n_ind, pl1, pl2);
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
List vs_biallelic_Rcpp(List PH,
                       IntegerMatrix G,
                       NumericMatrix pedigree) {
  NumericMatrix unique_pop_mat = retainUniqueAndSortByLastColumn(pedigree);
  //Rcpp::Rcout << unique_pop_mat << "\n";
  int n_fullsib_pop = unique_pop_mat.nrow();
  int n_ind = pedigree.nrow();
  NumericMatrix temp_mat = PH[0];
  int n_mrk = temp_mat.nrow();
  NumericVector pl1 = unique_pop_mat(_,2);
  NumericVector pl2 = unique_pop_mat(_,3);
  List result = calculate_L_and_initialize_H(n_fullsib_pop, n_mrk, n_ind, pl1, pl2);
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
      std::fill(temp_emit.begin(), temp_emit.end(), 1);
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


// [[Rcpp::export]]
List est_hmm_map_biallelic(List PH,
                           IntegerMatrix G,
                           NumericMatrix pedigree,
                           NumericVector rf,
                           bool verbose,
                           double tol,
                           bool ret_H0) {
  NumericVector p1 = pedigree(_,2)/2 - 1;
  NumericVector p2 = pedigree(_,3)/2 - 1;
  int n_ind = pedigree.nrow();
  NumericMatrix temp_mat = PH[0];
  int n_mrk = temp_mat.nrow();
  List result = vs_biallelic_Rcpp(PH, G, pedigree);
  List haplo = result["states"];
  List emit = result["emit"]; // Emission probabilities to be implemented
  std::vector<int> pl{2,4,6};
  //Initializing some variables
  int k, k1,  maxit = 1000, flag=0;
  double s, loglike=0.0, nr=0.0, temp=0.0;
  int mpi1 = max(p1);
  int mpi2 = max(p2);
  int max_ploidy_id;
  if(mpi1 >= mpi2)
    max_ploidy_id = mpi1;
  else
    max_ploidy_id = mpi2;
  std::vector<double> rf_cur(rf.size());
  std::vector<double> term(n_ind);
  std::fill(term.begin(), term.end(), 0.0);

  //Initializing v: states hmm should visit for each marker
  //Initializing e: emission probabilities associated to the states hmm should visit for each marker
  std::vector<std::vector<std::vector<int> > > v;
  std::vector<std::vector<std::vector<double> > > e;
  for(int i=0; i < haplo.size(); i++) // i: number of markers
  {
    Rcpp::List haplo_temp(haplo(i)); //states hmm should visit for marker i
    Rcpp::List emit_temp(emit(i)); //emission probs. for states hmm should visit for marker i
    std::vector<std::vector<int> > v1;
    std::vector<std::vector<double> > e1;
    for(int j=0; j < haplo_temp.size(); j++) //iterate for all j individuals
    {
      Rcpp::NumericMatrix M_temp = haplo_temp(j);
      Rcpp::NumericVector E_temp = emit_temp(j);
      std::vector<int> v2 = Rcpp::as<std::vector<int> >(M_temp);
      std::vector<double> e2 = Rcpp::as<std::vector<double> >(E_temp);
      v1.push_back(v2);
      e1.push_back(e2);
    }
    v.push_back(v1);
    e.push_back(e1);
  }
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
        Ttemp.push_back(transition(pl[j], rf_cur[i]));
      }
      T.push_back(Ttemp);
    }
    //Loop over all individuals
    for(int ind=0; ind < n_ind; ind++)
    {
      R_CheckUserInterrupt();
      for(int j=0; (unsigned)j < e[0][ind].size(); j++)
      {
        alpha[ind][0][j] = e[0][ind][j];
      }
      std::fill(beta[ind][n_mrk-1].begin(), beta[ind][n_mrk-1].end(), 1);
      //forward-backward
      for(k=1,k1=n_mrk-2; k < n_mrk; k++, k1--)
      {
        std::vector<double> temp4 (v[k][ind].size()/2);
        temp4 = forward_emit(alpha[ind][k-1],
                             v[k-1][ind],
                                   v[k][ind],
                                       e[k][ind],
                                           T[p1[ind]][k-1],
                                                     T[p2[ind]][k-1]);
        for(int j=0; (unsigned)j < temp4.size(); j++)
        {
          alpha[ind][k][j]=temp4[j];
        }
        std::vector<double> temp5 (v[k1][ind].size()/2);
        temp5=backward_emit(beta[ind][k1+1],
                            v[k1][ind],
                                 v[k1+1][ind],
                                        e[k1+1][ind],
                                               T[p1[ind]][k1],
                                                         T[p2[ind]][k1]);
        for(int j=0; (unsigned)j < temp5.size(); j++)
        {
          beta[ind][k1][j]=temp5[j];
        }
      }
      if(ret_H0 == 0)
      {
        //Updating recombination fraction
        for(int k = 0; k < n_mrk-1; k++)
        {
          vector<vector<double> > gamma(alpha[ind][k].size(), vector<double>(beta[ind][k+1].size()));
          s=0.0;
          int ngeni = alpha[ind][k].size();
          int ngenj = beta[ind][k+1].size();
          for(int i = 0; i < ngeni; i++)
          {
            for(int j = 0; j < ngenj; j++)
            {
              gamma[i][j] = alpha[ind][k][i] * beta[ind][k+1][j] *
                T[p1[ind]][k][v[k][ind][i]][v[k+1][ind][j]] *
                T[p2[ind]][k][v[k][ind][i+ngeni]][v[k+1][ind][j+ngenj]];
              if(i==0 && j==0) s = gamma[i][j];
              else s += gamma[i][j];
            }
          }
          for(int i=0; i < ngeni; i++)
          {
            for(int j=0; j < ngenj; j++)
            {
              nr=R[p1[ind]][v[k][ind][i]][v[k+1][ind][j]] +
                R[p2[ind]][v[k][ind][i+ngeni]][v[k+1][ind][j+ngenj]];
              if(s > 0) // Verify theoretical implications of this condition
                rf[k] +=  nr * gamma[i][j]/s;
            }
          }
        }
      }
      //Termination
      for(int j=0; (unsigned)j < alpha[ind][n_mrk-1].size(); j++)
      {
        term[ind] +=  alpha[ind][n_mrk-1][j];
      }
    } // loop over individuals
    //Likelihood using a specific recombination fraction vector
    //Usually, this is used to compute LOD Score under H0: rf=0.5
    if(ret_H0 == 1)
    {
      //Loglike computation
      for(int i=0; (unsigned)i < alpha.size(); i++)
      {
        temp=0.0;
        for(int j=0; (unsigned)j < alpha[i][n_mrk-1].size(); j++)
          temp += alpha[i][n_mrk-1][j];
        if(temp > 0)
          loglike += log10(temp);
      }
      if(verbose)
        Rcpp::Rcout << "\n";
      List z = List::create(wrap(loglike), rf_cur);
      return(z);
    }
    // rescale
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
    if(verbose)
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

  //Loglike computation
  for(int i=0; (unsigned)i < alpha.size(); i++)
  {
    temp=0.0;
    for(int j=0; (unsigned)j < alpha[i][n_mrk-1].size(); j++)
      temp += alpha[i][n_mrk-1][j];
    if(temp > 0)
      loglike += log10(temp);
  }
  if(verbose)
    Rcpp::Rcout << "\n";
  List z = List::create(wrap(loglike), rf);
  return(z);
}

