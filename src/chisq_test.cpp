#include <Rcpp.h>
#include  "segreg_poly.h"
 using namespace Rcpp;

 // [[Rcpp::export]]
 List mappoly_chisq_test(List input_data) {
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
   NumericVector::iterator it;
   NumericVector y(ploidy_pr + 1);

   for (int i = 0; i < n_mrk; ++i) {
     std::fill(y.begin(), y.end(), 0);
     for (int j = 0; j < n_ind; ++j) {
       y[geno_dose(i, j)] += 1;
     }
     it = std::remove(y.begin(), y.end(), NA_REAL);
     y.erase(it, y.end());

     NumericVector exp_seg = segreg_poly(ploidy_p1, ploidy_p2, d_p1[i], d_p2[i]);
     NumericVector y_filtered = y[exp_seg != 0];
     if(Rcpp::sum(y_filtered) == 0){
       chisq_p_out[i] = 10e-50;
     } else{
       NumericVector exp_seg_filtered = exp_seg[exp_seg != 0];
       Rcpp::Environment stats("package:stats");
       Rcpp::Function chisq_test = stats["chisq.test"];
       List x = chisq_test(Named("x", y_filtered), Named("p", exp_seg_filtered));
       chisq_p_out[i] = x["p.value"];
     }
   }
   chisq_p_out.names() = mrk_names;
   input_data["chisq.pval"] = chisq_p_out;
   return input_data;
 }

 // [[Rcpp::export]]
 List filter_non_conforming_classes(List input_data) {
   int ploidy_p1 = input_data["ploidy.p1"];
   int ploidy_p2 = input_data["ploidy.p2"];
   IntegerVector d_p1 = input_data["dosage.p1"];
   IntegerVector d_p2 = input_data["dosage.p2"];
   NumericMatrix geno_dose = input_data["geno.dose"];
   CharacterVector mrk_names = input_data["mrk.names"];
   int ploidy_pr = (ploidy_p1 + ploidy_p2) / 2;

   IntegerMatrix Dpop = cbind(d_p1, d_p2);
   NumericMatrix M(Dpop.nrow(), ploidy_pr + 1);
   for (int i = 0; i < Dpop.nrow(); ++i) {
     M(i, _) = segreg_poly(ploidy_p1, ploidy_p2, d_p1[i], d_p2[i]);
   }

   for (int i = 0; i < geno_dose.nrow(); ++i) {
     for (int j = 0; j < geno_dose.ncol(); ++j) {
       if (M(i, geno_dose(i, j)) == 0 || geno_dose(i, j) > ploidy_pr) {
         geno_dose(i, j) = NA_REAL;
       }
     }
   }

   input_data["geno.dose"] = geno_dose;
   return input_data;
 }

