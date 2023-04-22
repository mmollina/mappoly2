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

#include <Rcpp.h>
using namespace Rcpp;

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


