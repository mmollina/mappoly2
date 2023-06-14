#ifndef UTILS_H
#define UTILS_H

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

using namespace std;
using namespace Rcpp;

NumericVector segreg_poly(int ploidy_p1, int ploidy_p2, int d_p1, int d_p2);

#endif /* UTILS_H */
