#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
#include <string>

// Function to compute Mendelian segregation frequencies
Rcpp::NumericVector segreg_poly(int ploidy_p1, int ploidy_p2, int d_p1, int d_p2);

// Function to compute chi-squared test p-value
double cpp_chisq_test(Rcpp::NumericVector observed, Rcpp::NumericVector expected_probs);

// Function to compute chi-squared test for all markers
Rcpp::NumericVector mappoly_chisq_test(Rcpp::List input_data);

// Function to filter non-conforming genotypic classes
Rcpp::List filter_non_conforming_classes(Rcpp::List input_data);

// Utility functions to print matrices
void print_matrix(Rcpp::IntegerMatrix mat, std::string name);
void print_matrix(Rcpp::NumericMatrix mat, std::string name);

// Haldane map function
NumericVector imf_h(NumericVector r);

#endif // UTILS_H
