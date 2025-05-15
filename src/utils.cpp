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
#include <string>
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

/*
//' Mendelian segregation
 //'
 //' Computes the Mendelian segregation frequencies given the ploidy level
 //' of two parents and the dosage of the locus in both parents. It does
 //' not consider double reduction.
 //'
 //' @name segreg_poly
 //'
 //' @param ploidy.p1 ploidy level of parent 1
 //'
 //' @param ploidy.p2 ploidy level of parent 1
 //'
 //' @param d.p1 the dosage in parent 1
 //'
 //' @param d.p2 the dosage in parent 2
 //'
 //' @return a vector containing the expected segregation frequency for
 //'         genotypic classes.
 //'
 //' @examples
 //' seg1 <- segreg_poly(ploidy_p1 = 6, ploidy_p2 = 6, d_p1 = 3, d_p2 = 3)
 //' barplot(seg1)
 //' seg2 <- segreg_poly(ploidy_p1 = 2, ploidy_p2 = 4, d_p1 = 2, d_p2 = 3)
 //' barplot(seg2)
 //'
 //' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
 //'
 //' @export segreg_poly
 //' @importFrom Rcpp evalCpp
*/
 // [[Rcpp::export]]
 NumericVector segreg_poly(int ploidy_p1, int ploidy_p2, int d_p1, int d_p2) {

   if (d_p1 > ploidy_p1 || d_p2 > ploidy_p2)
     Rcpp::stop("dose should be smalled than ploidy");

   if (ploidy_p1 % 2 != 0 || ploidy_p2 % 2 != 0)
     Rcpp::stop("ploidy_p1 and ploidy_p2 must be even numbers");

   NumericVector p_dose(ploidy_p1 / 2 + ploidy_p2 / 2 + 1);
   NumericVector seg_p1 = dhyper(seq(0, ploidy_p1 + 1), d_p1, ploidy_p1 - d_p1, ploidy_p1 / 2);
   NumericVector seg_p2 = dhyper(seq(0, ploidy_p2 + 1), d_p2, ploidy_p2 - d_p2, ploidy_p2 / 2);
   int nrows = seg_p1.size();
   int ncols = seg_p2.size();
   NumericMatrix M(nrows, ncols);

   for (int i = 0; i < nrows; i++) {
     for (int j = 0; j < ncols; j++) {
       M(i, j) = seg_p1[i] * seg_p2[j];
     }
   }

   for (int i = 0; i < nrows; i++) {
     for (int j = 0; j < ncols; j++) {
       p_dose[i + j] += M(i, j);
     }
   }
   p_dose.names() = seq(0, p_dose.size() - 1);
   return p_dose;
 }

 // [[Rcpp::export]]
 double cpp_chisq_test(NumericVector observed, NumericVector expected_probs) {
   // Ensure expected_probs sums to 1
   double total_prob = sum(expected_probs);
   if (std::abs(total_prob - 1.0) > 1e-8) {
     stop("Expected probabilities do not sum to 1.");
   }

   // Compute total count from observed
   double total_count = sum(observed);

   // Scale expected probabilities to expected counts
   NumericVector expected_counts = expected_probs * total_count;

   // Check for zero expected counts to avoid division by zero
   if (is_true(any(expected_counts == 0))) {
     stop("Expected counts contain zeros, invalid for chi-squared test.");
   }

   // Compute chi-squared statistic
   double chi_squared = 0.0;
   for (int i = 0; i < observed.size(); ++i) {
     double obs = observed[i];
     double exp = expected_counts[i];
     chi_squared += std::pow(obs - exp, 2) / exp;
   }

   // Degrees of freedom: number of categories minus 1
   int df = observed.size() - 1;

   // Compute p-value using the chi-squared distribution CDF
   double p_value = R::pchisq(chi_squared, df, false, false); // lower.tail = false, log.p = false

   return p_value;
 }

 // [[Rcpp::export]]
 NumericVector mappoly_chisq_test(List input_data) {
   int ploidy_p1 = input_data["ploidy.p1"];
   int ploidy_p2 = input_data["ploidy.p2"];
   int n_ind = input_data["n.ind"];
   int n_mrk = input_data["n.mrk"];
   IntegerVector d_p1 = input_data["dosage.p1"];
   IntegerVector d_p2 = input_data["dosage.p2"];
   NumericMatrix geno_dose = input_data["geno.dose"];
   CharacterVector mrk_names = input_data["mrk.names"];
   int ploidy_pr = (ploidy_p1 + ploidy_p2) / 2;

   NumericVector chisq_p_out(n_mrk);
   NumericVector y(ploidy_pr + 1);

   for (int i = 0; i < n_mrk; ++i) {
     std::fill(y.begin(), y.end(), 0); // Reset y for the current marker

     // Count genotype doses, skipping NA values
     for (int j = 0; j < n_ind; ++j) {
       if (!R_IsNA(geno_dose(i, j))) { // Check if the value is not NA
         y[geno_dose(i, j)] += 1;
       }
     }

     NumericVector exp_seg = segreg_poly(ploidy_p1, ploidy_p2, d_p1[i], d_p2[i]);
     LogicalVector mask = exp_seg != 0;
     NumericVector y_filtered = y[mask];

     if (Rcpp::sum(y_filtered) == 0) {
       chisq_p_out[i] = 1e-50; // Assign a small p-value for invalid cases
     } else {
       NumericVector exp_seg_filtered = exp_seg[mask];
       exp_seg_filtered = exp_seg_filtered / Rcpp::sum(exp_seg_filtered); // Normalize

       // Use the cpp_chisq_test function
       chisq_p_out[i] = cpp_chisq_test(y_filtered, exp_seg_filtered);
     }
   }

   chisq_p_out.names() = mrk_names;
   return chisq_p_out;
 }

 // [[Rcpp::export]]
 List filter_non_conforming_classes(List input_data) {
   int ploidy_p1 = input_data["ploidy.p1"];
   int ploidy_p2 = input_data["ploidy.p2"];
   IntegerVector d_p1 = input_data["dosage.p1"];
   IntegerVector d_p2 = input_data["dosage.p2"];
   NumericMatrix geno_dose = input_data["geno.dose"];
   CharacterVector mrk_names = input_data["mrk.names"];
   int ploidy_pr = (ploidy_p1 + ploidy_p2) / 2;  // Progeny ploidy

   for (int i = 0; i < geno_dose.nrow(); ++i) {
     int d1 = d_p1[i];
     int d2 = d_p2[i];

     // Valid dosage range using the no double reduction rule
     int min_valid = std::max(0, d1 + d2 - ploidy_pr);
     int max_valid = std::min(d1 + d2, ploidy_pr);

     for (int j = 0; j < geno_dose.ncol(); ++j) {
       double x = geno_dose(i, j);
       if (NumericVector::is_na(x)) continue;

       // x must be an integer and within the valid range
       if (x < min_valid || x > max_valid || floor(x) != x) {
         geno_dose(i, j) = NA_REAL;
       }
     }
   }

   input_data["geno.dose"] = geno_dose;
   return input_data;
 }



 // Function to print a NumericMatrix or IntegerMatrix
 void print_matrix(NumericMatrix mat, std::string name) {
   Rcpp::Rcout << name << " (" << mat.nrow() << " x " << mat.ncol() << "):\n";
   for (int i = 0; i < mat.nrow(); ++i) {
     for (int j = 0; j < mat.ncol(); ++j) {
       Rcpp::Rcout << mat(i, j) << "";
     }
     Rcpp::Rcout << "\n";
   }
   Rcpp::Rcout << "\n";
 }

 void print_matrix(IntegerMatrix mat, std::string name) {

   Rcpp::Rcout << name << " (" << mat.nrow() << " x " << mat.ncol() << "):\n";
   for (int i = 0; i < mat.nrow(); ++i) {
     for (int j = 0; j < mat.ncol(); ++j) {
       Rcpp::Rcout << mat(i, j) << "";
     }
     Rcpp::Rcout << "\n";
   }
   Rcpp::Rcout << "\n";
 }


 NumericVector imf_h(NumericVector r) {
   NumericVector d = clone(r); // Clone the input vector to avoid modifying the original
   for(int i = 0; i < d.size(); ++i) {
     if(d[i] >= 0.5) {
       d[i] = 0.5 - 1e-14;
     }
     d[i] = -50 * std::log(1 - 2 * d[i]);
   }
   return d;
 }
