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
#include "utils.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
double calc_loglike_hmm_rcpp(List PH,
                           IntegerMatrix G,
                           NumericMatrix pedigree,
                           NumericVector rf_vec) {
  NumericVector ploidy_p1 = pedigree(_,2)/2 - 1;
  NumericVector ploidy_p2 = pedigree(_,3)/2 - 1;
  int n_ind = pedigree.nrow();
  NumericMatrix temp_phase_mat = PH[0];
  int n_mrk = temp_phase_mat.nrow();
  std::vector<int> pl{2,4,6};
  double s, loglike = 0.0;

  // Getting the maximum ploidy level among founders
  int mpi1 = max(ploidy_p1);
  int mpi2 = max(ploidy_p2);
  int max_ploidy_id;
  if(mpi1 >= mpi2)
    max_ploidy_id = mpi1;
  else
    max_ploidy_id = mpi2;


  // HMM states that should be visited given the phase of
  // the founders, genotype of the offspring and pedigree
  List result = vs_biallelic_Rcpp(PH, G, pedigree);
  List ve = hmm_vectors(result);
  std::vector<std::vector<std::vector<int> > > v = ve["v"];
  std::vector<std::vector<std::vector<double> > > e = ve["e"];

  //Initializing alpha and beta
  std::vector<std::vector<std::vector<double> > > alpha(n_ind);
  for(int ind=0; ind < n_ind; ind++)
  {
    for(int i=0; i < n_mrk; i++)
    {
      std::vector<double> temp3(v[i][ind].size()/2);
      alpha[ind].push_back(temp3);
    }
  }

  std::vector<double> cur_rf(rf_vec.size());
  std::fill(cur_rf.begin(), cur_rf.end(), 0.0);
  for(int i = 0; i < cur_rf.size(); i++){
    cur_rf[i] = rf_vec[i];
  }

  //Initializing transition matrices
  std::vector< std::vector< std::vector< std::vector<double> > > >T;
  for(int j=0; j <= max_ploidy_id; j++)
  {
    std::vector< std::vector< std::vector<double> > > Ttemp;
    for(int i=0; i < n_mrk-1; i++)
    {
      Ttemp.push_back(transition(pl[j], cur_rf[i]));
    }
    T.push_back(Ttemp);
  }
  //Loop over all individuals
  for(int ind=0; ind < n_ind; ind++)
  {
    R_CheckUserInterrupt();
    for(int j=0; (unsigned)j < e[0][ind].size(); j++)
      alpha[ind][0][j] = e[0][ind][j];

    //forward
    for(int k=1; k < n_mrk; k++)
    {
      std::vector<double> temp4(v[k][ind].size()/2);
      temp4 = forward_emit(alpha[ind][k-1],
                           v[k-1][ind],
                                 v[k][ind],
                                     e[k][ind],
                                         T[ploidy_p1[ind]][k-1],
                                                          T[ploidy_p2[ind]][k-1]);

      for(int j=0; (unsigned)j < temp4.size(); j++)
        alpha[ind][k][j]=temp4[j];
    }
  }
  for(int i=0; (unsigned)i < alpha.size(); i++)
  {
    double temp = alpha[i][n_mrk-1][0];
    for(int j=1; (unsigned)j < alpha[i][n_mrk-1].size(); j++)
      temp += alpha[i][n_mrk-1][j];
    if(temp > 1e-50)
      loglike += log(temp);
  }
  return(loglike);
}
