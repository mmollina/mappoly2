/*
 MAPpoly2: a package to construct genetic maps in autopolyploids
 Copyright (C) 2023 Marcelo Mollinari

 This file is part of MAPpoly2.

 MAPpoly2 is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 For a copy of the GNU General Public License, please visit
 <http://www.gnu.org/licenses/>.
 */

/*
 File: two_pts_est.cpp

 Functions Written partially by Marcelo Mollinari.

 Part of this function was adapted from Brent_fmin function,
 which can be found in R/src/library/stats/src/optimize.c
 Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 Copyright (C) 2003-2004  The R Foundation
 Copyright (C) 1998--2014-2018 The R Core Team

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version: December, 2013
 Last update: April, 2023
 */

#include <float.h> /* DBL_EPSILON */
#include <R_ext/Applic.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <Rmath.h>
#include <Rcpp.h>
#include <R_ext/PrtUtil.h>
#include "combinatorics.h"

using namespace std;
using namespace Rcpp;

/* recombination frequency likelihood*/
double twopt_likelihood(double rf,
                        int ploidy_p1,
                        int ploidy_p2,
                        Rcpp::NumericVector dk,
                        Rcpp::NumericVector dk1,
                        Rcpp::NumericVector gen_1,
                        Rcpp::NumericVector gen_2,
                        Rcpp::NumericMatrix count_mat)
{
  int count, count2=0;
  double temp=0.0;
  int offspring_ploidy = (ploidy_p1 + ploidy_p2)/2;
  Rcpp::NumericMatrix Tr(offspring_ploidy+2, offspring_ploidy+2);
  std::fill(Tr.begin(), Tr.end(), 1);
  for(int i = 0; i < dk.size(); i++){
    count=0;
    Tr(dk(i),dk1(i))=0;
    for(int l_p1 = 0; l_p1 <= ploidy_p1/2; l_p1++){
      for(int l_p2 = 0; l_p2 <= ploidy_p2/2; l_p2++) {
        Tr(dk(i),dk1(i))=Tr(dk(i),dk1(i)) +
          count_mat(count2, count) *
          (pow((1-rf), (ploidy_p1/2)-l_p1) * pow(rf, l_p1))/(nChoosek(ploidy_p1/2,l_p1)) *
          (pow((1-rf), (ploidy_p2/2)-l_p2) * pow(rf, l_p2))/(nChoosek(ploidy_p2/2,l_p2));
        count++;
      }
    }
    count2++;
  }
  for(int i = 0; i < gen_1.size(); i++)
    temp += log10(Tr(gen_1(i),gen_2(i)));
  return -temp ;
}
/*recombination frequency estimation*/
RcppExport SEXP pairwise_rf_estimation(SEXP ploidy_p1_R,
                                       SEXP ploidy_p2_R,
                                       SEXP mrk_pairs_R,
                                       SEXP geno_R,
                                       SEXP d_p1_R,
                                       SEXP d_p2_R,
                                       SEXP count_cache_R,
                                       SEXP tol_R,
                                       SEXP swap_parents_R)
{
  Rcpp::NumericMatrix mrk_pairs = Rcpp::as<Rcpp::NumericMatrix>(mrk_pairs_R);
  Rcpp::NumericMatrix geno = Rcpp::as<Rcpp::NumericMatrix>(geno_R);
  Rcpp::NumericVector d_p1 = Rcpp::as<Rcpp::NumericVector>(d_p1_R);
  Rcpp::NumericVector d_p2 = Rcpp::as<Rcpp::NumericVector>(d_p2_R);
  Rcpp::List count_cache = Rcpp::as<Rcpp::List>(count_cache_R);
  bool swap_parents = Rcpp::as<bool>(swap_parents_R);
  Rcpp::NumericVector d_pair(4);
  Rcpp::List out(mrk_pairs.ncol());
  int ploidy_p1 = Rcpp::as<int>(ploidy_p1_R);
  int ploidy_p2 = Rcpp::as<int>(ploidy_p2_R);
  double tol = Rcpp::as<double>(tol_R);
  for(int k=0; k < mrk_pairs.ncol(); k++)
  {
    //Rcpp::Rcout << mrk_pairs(0,k)+1 << " - " << mrk_pairs(1,k)+1 <<  std::endl;
    int id = (ploidy_p1+1) *
        (ploidy_p1+1) *
        (ploidy_p2+1) *
        d_p2[mrk_pairs(1,k)] +
        (ploidy_p1+1) *
        (ploidy_p1+1) *
        d_p2[mrk_pairs(0,k)] +
        (ploidy_p1+1) *
        d_p1[mrk_pairs(1,k)] +
        d_p1[mrk_pairs(0,k)] +
        1;
    //Rcpp::Rcout << "id: " << id <<  std::endl;
    Rcpp::List temp_list = count_cache[(id-1)];
    //Rcpp::List temp_list = count_cache[k];
    if(temp_list.size() > 1)
    {
      NumericVector gen_1 = geno( mrk_pairs(0,k), _);
      NumericVector gen_2 = geno( mrk_pairs(1,k), _);
      Rcpp::NumericMatrix res(3, temp_list.size());
      Rcpp::CharacterVector zn = temp_list.attr( "names" ) ;
      if(swap_parents) { // if ploidy.p1 > ploidy.p2
        Rcpp::CharacterVector zn_rev(zn.size());
        //Rcpp::Rcout << "before: " << zn << "-->";
        for(int i=0; i < zn_rev.size(); i++){
          std::string zn_temp = Rcpp::as<string>(zn(i));
          reverse(zn_temp.begin(), zn_temp.end());
          zn_rev(i) = zn_temp;
        }
       colnames(res) = zn_rev;
       //Rcpp::Rcout << "after: " << zn << "\n";
      } else{
        colnames(res) = zn;
      }
      for(int i=0; i < temp_list.size(); i++)
      {
        Rcpp::NumericMatrix count_mat = temp_list[i] ;
        Rcpp::List dimnames = count_mat.attr( "dimnames" ) ;
        Rcpp::CharacterVector z = dimnames[0];
        Rcpp::NumericVector dk(z.size()), dk1(z.size());
        std::string delimiter = " ";
        for(int j=0; j < z.size(); j++)
        {
          std::string lnames = Rcpp::as<std::string>(z(j));
          dk(j) = std::stoi(lnames.substr(0,lnames.find(delimiter)));
          dk1(j) = std::stoi(lnames.substr(lnames.find(delimiter)+1, lnames.length()));
        }
        //an approximation  x  to the point where  f  attains a minimum  on
        //the interval  (a,b)  is determined.
        // Adapted from Brent_fmin function, which can be found in R/src/library/stats/src/optimize.c
        //   Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
        //   Copyright (C) 2003-2004  The R Foundation
        //   Copyright (C) 1998--2014-2018 The R Core Team
        // This function subprogram is a slightly modified  version  of  the
        // Algol  60 procedure  localmin  given in Richard Brent, Algorithms for
        // Minimization without Derivatives, Prentice-Hall, Inc. (1973).
        // Brent's Minimization Procedure starts
        //  c is the squared inverse of the golden ratio
        const double c = (3. - sqrt(5.)) * .5;
        // Local variables
        double a, b, d, e, p, q, r, u, v, w, x;
        double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
        //  eps is approximately the square root of the relative machine precision.
        eps = DBL_EPSILON;
        tol1 = eps + 1.;// the smallest 1.000... > 1
        eps = sqrt(eps);

        a = 0.0;
        b = 0.5;
        v = a + c * (b - a);
        w = v;
        x = v;

        d = 0.;// -Wall
        e = 0.;
        fx = twopt_likelihood(x,
                              ploidy_p1,
                              ploidy_p2,
                              dk,
                              dk1,
                              gen_1,
                              gen_2,
                              count_mat);
        fv = fx;
        fw = fx;
        tol3 = tol / 3.;

        //  main loop starts here -----------------------------------

        for(;;) {
          xm = (a + b) * .5;
          tol1 = eps * fabs(x) + tol3;
          t2 = tol1 * 2.;

          // check stopping criterion

          if (fabs(x - xm) <= t2 - (b - a) * .5) break;
          p = 0.;
          q = 0.;
          r = 0.;
          if (fabs(e) > tol1) { // fit parabola

            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.;
            if (q > 0.) p = -p; else q = -q;
            r = e;
            e = d;
          }

          if (fabs(p) >= fabs(q * .5 * r) ||
              p <= q * (a - x) || p >= q * (b - x)) { // a golden-section step

            if (x < xm) e = b - x; else e = a - x;
            d = c * e;
          }
          else { // a parabolic-interpolation step

            d = p / q;
            u = x + d;

            // f must not be evaluated too close to ax or bx

            if (u - a < t2 || b - u < t2) {
              d = tol1;
              if (x >= xm) d = -d;
            }
          }

          // f must not be evaluated too close to x

          if (fabs(d) >= tol1)
            u = x + d;
          else if (d > 0.)
            u = x + tol1;
          else
            u = x - tol1;

          fu = twopt_likelihood(u,
                                ploidy_p1,
                                ploidy_p2,
                                dk,
                                dk1,
                                gen_1,
                                gen_2,
                                count_mat);

          //  update  a, b, v, w, and x

          if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w;    w = x;   x = u;
            fv = fw; fw = fx; fx = fu;
          } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
              v = w; fv = fw;
              w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
              v = u; fv = fu;
            }
          }
        }
        // Brent's Minimization Procedure ends
        res(0,i) = x;
        res(1,i) = twopt_likelihood(x,
            ploidy_p1,
            ploidy_p2,
            dk,
            dk1,
            gen_1,
            gen_2,
            count_mat);
        res(2,i) = twopt_likelihood(0.5,
            ploidy_p1,
            ploidy_p2,
            dk,
            dk1,
            gen_1,
            gen_2,
            count_mat);
      }
      out(k)=res;
    }
    else
    {
      if(swap_parents) {
        Rcpp::NumericVector d_out(4);
        d_out(0)=d_p2[mrk_pairs(0,k)];
        d_out(1)=d_p2[mrk_pairs(1,k)];
        d_out(2)=d_p1[mrk_pairs(0,k)];
        d_out(3)=d_p1[mrk_pairs(1,k)];
        out(k)=d_out;
      } else {
        Rcpp::NumericVector d_out(4);
        d_out(0)=d_p1[mrk_pairs(0,k)];
        d_out(1)=d_p1[mrk_pairs(1,k)];
        d_out(2)=d_p2[mrk_pairs(0,k)];
        d_out(3)=d_p2[mrk_pairs(1,k)];
        out(k)=d_out;
      }
    }
  }
  return(out);
}
//end of file two_pts_est.cpp
