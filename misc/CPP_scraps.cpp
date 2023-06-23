

// if(it ==0 && ind == 0 && k == 1 && s1 < 5)
//   Rcout << "T ["<<ploidy_p1[ind]<<"]["<<k-1<<"]["<<v[k-1][ind][0]<<"]["<< v[k][ind][s1] <<"] : " <<
//     T[ploidy_p1[ind]][k-1][v[k-1][ind][0]][v[k][ind][s1]] <<
//       " + T ["<<ploidy_p1[ind]<<"]["<<k-1<<"]["<<v[k-1][ind][0+ngen_k0]<<"]["<< v[k][ind][s1+ngen_k1] <<"] : " <<
//         T[ploidy_p2[ind]][k-1][v[k-1][ind][0+ngen_k0]][v[k][ind][s1+ngen_k1]] << " = " <<
//           T[ploidy_p1[ind]][k-1][v[k-1][ind][0]][v[k][ind][s1]] +
//           T[ploidy_p2[ind]][k-1][v[k-1][ind][0+ngen_k0]][v[k][ind][s1+ngen_k1]] << "==>" <<
//             "alpha["<<ind<<"]["<<k<<"]["<<s1<<"] : " <<  alpha[ind][k][s1] <<
//               "\n";


// if(it ==0 && ind == 0 && k == 1&& s1 < 5){
//   Rcout << "~~~~~~~~~~~~~~~~~~~~~INIT~~~~~~~~~~~~~~~~~~~~~~~\n";
//   Rcout << "alpha["<<ind<<"]["<<k<<"]["<<s1<<"] : " <<  alpha[ind][k][s1] << "\n";
//   Rcout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
// }


// if(it ==0 && ind == 0 && k == 1 && s1 < 5)
//     Rcout << "T ["<<ploidy_p1[ind]<<"]["<<k-1<<"]["<<v[k-1][ind][s0]<<"]["<< v[k][ind][s1] <<"] : " <<
//       T[ploidy_p1[ind]][k-1][v[k-1][ind][s0]][v[k][ind][s1]] <<
//          " + T ["<<ploidy_p1[ind]<<"]["<<k-1<<"]["<<v[k-1][ind][s0+ngen_k0]<<"]["<< v[k][ind][s1+ngen_k1] <<"] : " <<
//           T[ploidy_p2[ind]][k-1][v[k-1][ind][s0+ngen_k0]][v[k][ind][s1+ngen_k1]] << " = " <<
//             T[ploidy_p1[ind]][k-1][v[k-1][ind][s0]][v[k][ind][s1]] +
//             T[ploidy_p2[ind]][k-1][v[k-1][ind][s0+ngen_k0]][v[k][ind][s1+ngen_k1]] << "==>" <<
//               "alpha["<<ind<<"]["<<k<<"]["<<s1<<"] : " <<  alpha[ind][k][s1] <<
//               "\n";



// if(it ==0 && ind == 0 && k == 1&& s1 < 5){
//   Rcout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
//   Rcout << "alpha["<<ind<<"]["<<k<<"]["<<s1<<"] * emit : " <<  alpha[ind][k][s1] << "\n";
//   Rcout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
// }


//if(it ==0 && ind == 0 && k ==2)
//  Rcout << "T for alpha["<<ind<<"]["<<k<<"]["<<s1<<"] : " <<
//    T[ploidy_p1[ind]][k-1][v[k-1][ind][s1]][v[k][ind][s1]] +
//    T[ploidy_p2[ind]][k-1][v[k-1][ind][s0+ngen_k0]][v[k][ind][s1+ngen_k1]] << "\n";


//if(it ==0)
//  printAlphaBeta(alpha, beta, 0);

//double zeta = 0.0;
//for(int j=0; j < alpha[ind][k].size(); j++)
//zeta += alpha[ind][k][j];
//for(int j=0; j < alpha[ind][k].size(); j++)
//   alpha[ind][k][j] /= zeta;


//double epsilon = 0.0;
// for(int j=0; j < beta[ind][k1].size(); j++)
//  epsilon += beta[ind][k1][j];
// for(int j=0; j < beta[ind][k1].size(); j++)
//   beta[ind][k1][j] /= epsilon;


void printAlphaBeta(vector<vector<vector<double>>>& alpha,
                    vector<vector<vector<double>>>& beta,
                    int ind) {
  int size_alpha_ind = alpha[ind].size();
  int size_beta_ind = beta[ind].size();
  int size_alpha_ind_k, size_beta_ind_k;

  cout << fixed << setprecision(15);

  for(int k = 0; k < size_alpha_ind; k++) {
    size_alpha_ind_k = alpha[ind][k].size();
    for(int j = 0; j < size_alpha_ind_k; j++) {
      cout << "alpha[" << ind << "][" << k << "][" << j << "] = " << alpha[ind][k][j] << endl;
    }
  }

  //for(int k1 = 0; k1 < size_beta_ind; k1++) {
  //  size_beta_ind_k = beta[ind][k1].size();
  //  for(int j = 0; j < size_beta_ind_k; j++) {
  //    cout << "beta[" << ind << "][" << k1 << "][" << j << "] = " << beta[ind][k1][j] << endl;
  //  }
  //}
}

