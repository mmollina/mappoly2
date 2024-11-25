/*
 MAPpoly: a package to construct genetic maps in autopolyploids
 Copyright (C) 2014-2023 Marcelo Mollinari

 This file is part of MAPpoly.

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
 File: twopt_phasing.cpp

 Description:

 Functions Written by Marcelo Mollinari.

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version:       2022
 Last update: May 25, 2023
 */

// [[Rcpp::plugins(cpp11)]]
#include <R.h>
#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include "combinatorics.h"
#include "utils.h"
#include "hmm_map.h"
#include <math.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <unordered_set>
#include <progress.hpp>
#define THRESH 200.0
using namespace Rcpp;
using namespace std;

// Helper function for list concatenation
List cLists(List x, List y) {
  int nsize = x.size();
  int msize = y.size();
  List out(nsize + msize);

  for(int i = 0; i < nsize; i++) {
    out[i] = x[i];
  }
  for(int i = 0; i < msize; i++) {
    out[nsize+i] = y[i];
  }
  return(out);
}

// Helper function for matrix binding
NumericMatrix rbind_cpp(NumericMatrix a, NumericVector b){
  int nra = a.nrow(), nca = a.ncol();
  NumericMatrix out(nra+1, nca);

  // Fill the matrix
  for(int i=0; i < nra; i++)
    out(i, _) = a(i, _);
  out(nra, _) = b;

  return out;
}

// Helper function to initiate phasing matrix
NumericMatrix init_phase_mat(int d, int p){
  NumericMatrix v(1, p);

  if(d == 0)
    return v;
  else {
    for(int i = 0; i < d; i++)
      v(0,i) = 1.0;
    return v;
  }
}

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::export]]
List twopt_phasing_cpp(CharacterVector mrk_id,
                       CharacterVector seg_mrk_id,
                       int ploidy,
                       IntegerVector dose_vec,
                       NumericMatrix S,
                       int max_conf_number,
                       IntegerMatrix G,
                       double tol,
                       double err,
                       int tail,
                       double hmm_thresh,
                       double map_expansion_thresh,
                       bool verbose) {

  // Initialize index vector to keep track of phased markers
  IntegerVector idx = IntegerVector::create(0);

  // Initialize list of haplotype configurations with the first marker
  List H(1);
  H[0] = init_phase_mat(dose_vec[0], ploidy);

  int n = mrk_id.size(); // Total number of markers
  Progress p(n, verbose); // Initialize progress bar

  // Create a mapping from marker names to their indices in mrk_id
  std::map<std::string, int> mrk_id_map;
  for(int i = 0; i < n; ++i) {
    mrk_id_map[as<std::string>(mrk_id[i])] = i;
  }

  // Convert seg_mrk_id to indices into mrk_id
  int m = seg_mrk_id.size(); // Number of segregating markers
  IntegerVector seg_mrk_idx(m);
  for(int i = 0; i < m; ++i) {
    std::string name = as<std::string>(seg_mrk_id[i]);
    seg_mrk_idx[i] = mrk_id_map[name]; // Index into mrk_id
  }

  // Create a mapping from mrk_id indices to positions in G (row indices)
  // G_row_index[i] = position of marker i in seg_mrk_idx (i.e., row in G)
  std::vector<int> G_row_index(n, -1); // Initialize with -1
  for(int i = 0; i < m; ++i) {
    int mrk_idx = seg_mrk_idx[i]; // Index into mrk_id
    G_row_index[mrk_idx] = i;     // Position in G (row index)
  }

  // Initialize a mapping from mrk_id indices to positions in idx (for PH indexing)
  // We'll update this as idx grows
  std::map<int, int> idx_pos_in_idx; // Key: mrk_idx, Value: position in idx
  idx_pos_in_idx[idx[0]] = 0; // idx[0] is 0 (first marker)

  // Initialize previous cumulative distance
  double prev_cumulative_distance = 0.0;

  // Main loop: Iterate over each marker starting from the second one
  for(int i = 1; i < n; ++i) {

    // Check for user interrupt to allow aborting the process
    if (Progress::check_abort())
      return List::create(); // Return an empty list if aborted

    int x = dose_vec[i]; // Dosage of the current marker

    // Case 1: Marker dosage equals ploidy (homozygous presence)
    if(x == ploidy) {
      for(int j = 0; j < H.size(); ++j) {
        // Append a row of ones to each configuration
        H[j] = rbind_cpp(as<NumericMatrix>(H[j]), NumericVector(ploidy, 1.0));
      }
      idx.push_back(i); // Add marker index to the list
      idx_pos_in_idx[i] = idx.size() - 1; // Position in idx

      // Case 2: Marker dosage is zero (homozygous absence)
    } else if(x == 0) {
      for(int j = 0; j < H.size(); ++j) {
        // Append a row of zeros to each configuration
        H[j] = rbind_cpp(as<NumericMatrix>(H[j]), NumericVector(ploidy, 0.0));
      }
      idx.push_back(i); // Add marker index to the list
      idx_pos_in_idx[i] = idx.size() - 1; // Position in idx

      // Case 3: Marker has intermediate dosage (heterozygous)
    } else {
      // Vector to store distances or similarities with previously phased markers
      NumericVector d(idx.size());
      for(int j = 0; j < idx.size(); ++j) {
        d[j] = S(i, idx[j]); // Extract values from S (could be distances or other metrics)
      }

      List vtemp(H.size()); // Temporary list to store valid permutations
      List Hres; // List to accumulate new haplotype configurations

      // Iterate over existing haplotype configurations
      for(int j = 0; j < H.size(); ++j) {
        // Find valid allele permutations for the current configuration
        vtemp[j] = find_valid_permutations(as<NumericMatrix>(H[j]), d, x);
        NumericMatrix z = vtemp[j]; // Matrix of valid permutations

        // If no valid permutations, skip to next configuration
        if(z.nrow() == 0) continue;

        List Htemp(z.nrow()); // Temporary list to store new configurations

        // Append each valid permutation to the current configuration
        for(int k = 0; k < z.nrow(); ++k) {
          Htemp[k] = rbind_cpp(as<NumericMatrix>(H[j]), z(k, _));
        }

        // Filter new configurations and add them to the result list
        Hres = cLists(Hres, filter_matrices(Htemp));
      }

      // Further filter configurations if more than one exists
      if(Hres.size() > 1)
        Hres = filter_matrices(Hres);

      // Skip marker if the number of configurations exceeds the maximum allowed
      if(Hres.size() > max_conf_number || Hres.size() == 0) {
        // Optionally, print a message if verbose is true
        //if(verbose) Rcout << "Skipping marker " << as<std::string>(mrk_id[i]) << "\n";
        continue; // Skip to the next marker
      }

      // Now, for each matrix in Hres, call est_hmm_map_biallelic_single
      // Prepare to collect results
      List est_results(Hres.size()); // To store the results of est_hmm_map_biallelic_single

      // Vectors to store loglikes and cumulative distances
      std::vector<double> loglikes(Hres.size());
      std::vector<double> cumulative_distances(Hres.size());

      // Determine markers that are in both idx (phased markers) and seg_mrk_idx (segregating for the parent)
      std::vector<int> segregating_markers;
      for(int idx_i = 0; idx_i < idx.size(); ++idx_i) {
        int mrk_idx = idx[idx_i]; // Index into mrk_id
        if(G_row_index[mrk_idx] != -1) { // Marker is in seg_mrk_id
          segregating_markers.push_back(mrk_idx);
        }
      }
      // Include current marker if it is segregating for the parent
      if(G_row_index[i] != -1) {
        segregating_markers.push_back(i);
      }

      // Determine the tail length (number of markers to consider)
      int total_seg_markers = segregating_markers.size();
      int tail_length = tail;
      if(total_seg_markers < tail) {
        tail_length = total_seg_markers;
      }

      // Get the indices of the markers in the tail (last 'tail_length' segregating markers)
      std::vector<int> tail_idx(tail_length);
      for(int k = 0; k < tail_length; ++k) {
        tail_idx[k] = segregating_markers[total_seg_markers - tail_length + k];
      }

      // Only proceed if we have at least two markers (required for HMM)
      if(tail_length >= 2) {
        // Set recombination fractions (rf) to 0.01 for all marker pairs in the tail
        NumericVector rf(tail_length - 1, 0.01);

        // For each configuration in Hres
        for(int j = 0; j < Hres.size(); ++j) {
          NumericMatrix PH = as<NumericMatrix>(Hres[j]); // Current phase configuration

          // Construct PH_tail
          NumericMatrix PH_tail(tail_length, ploidy);
          for(int row = 0; row < tail_length; ++row) {
            int mrk_idx = tail_idx[row]; // Index into mrk_id
            int ph_row;

            if (mrk_idx == i) {
              // Current marker is at the last row in PH
              ph_row = PH.nrow() - 1;
            } else if (idx_pos_in_idx.find(mrk_idx) != idx_pos_in_idx.end()) {
              ph_row = idx_pos_in_idx[mrk_idx];
            } else {
              Rcpp::stop("Marker not found in idx_pos_in_idx");
            }

            PH_tail(row, _) = PH(ph_row, _);
          }

          // Extract G_tail as before
          IntegerMatrix G_tail(tail_length, G.ncol());
          for(int row = 0; row < tail_length; ++row) {
            int mrk_idx = tail_idx[row]; // Index into mrk_id
            int g_row = G_row_index[mrk_idx]; // Row index in G
            if(g_row == -1) {
              Rcpp::stop("Marker not found in G");
            }
            G_tail(row, _) = G(g_row, _);
          }

          // Proceed with HMM estimation
          List est_result = est_hmm_map_biallelic_single(PH_tail, G_tail, rf, err, false, false, tol, false);

          // Extract loglike and rf_cur from est_result
          double loglike = as<double>(est_result[0]);
          NumericVector rf_cur = est_result[1];

          // Apply Haldane mapping function to rf_cur to get distances
          NumericVector dist_cur = imf_h(rf_cur);

          // Compute cumulative distance
          double cumulative_distance = std::accumulate(dist_cur.begin(), dist_cur.end(), 0.0);

          // Store loglike and cumulative_distance
          loglikes[j] = loglike;
          cumulative_distances[j] = cumulative_distance;

          // Print distance vector and loglike with rounding to 2 decimal points
          // Rcout << "Marker " << as<std::string>(mrk_id[i]) << ", Configuration " << j + 1 << ":\n";
          // Rcout << "Distance (cM): ";
          // for(int k = 0; k < dist_cur.size(); ++k) {
          //   Rcout << std::fixed << std::setprecision(2) << dist_cur[k] << " ";
          // }
          // Rcout << "\n";
          // Rcout << "loglike: " << std::fixed << std::setprecision(2) << loglike << "\n";

          // Save the result in the est_results list
          est_results[j] = est_result;
        }

        // After collecting est_results, process them according to thresholds

        // Find the maximum loglike
        double max_loglike = *std::max_element(loglikes.begin(), loglikes.end());
        size_t best_config_index = std::distance(loglikes.begin(), std::max_element(loglikes.begin(), loglikes.end()));

        // Filter configurations based on thresholds
        std::vector<int> configs_to_keep;

        for(size_t j = 0; j < Hres.size(); ++j) {
          double loglike_diff = max_loglike - loglikes[j];
          double expansion = cumulative_distances[j] - prev_cumulative_distance;

          if(loglike_diff <= hmm_thresh && expansion <= map_expansion_thresh) {
            configs_to_keep.push_back(j);
          }
        }

        // Check if any configurations remain
        if(configs_to_keep.empty()) {
          //if(verbose) Rcout << "Marker " << as<std::string>(mrk_id[i]) << " removed due to thresholds\n";
          continue; // Skip to the next marker
        } else {
          // Update H to keep only configurations in configs_to_keep
          List H_new(configs_to_keep.size());
          for(size_t k = 0; k < configs_to_keep.size(); ++k) {
            H_new[k] = Hres[configs_to_keep[k]];
          }
          H = H_new;

          // Update prev_cumulative_distance
          prev_cumulative_distance = cumulative_distances[best_config_index];

          // Update idx and idx_pos_in_idx
          idx.push_back(i);
          idx_pos_in_idx[i] = idx.size() - 1;
        }

      } else {
        // Not enough markers to perform HMM estimation
        // Proceed without applying thresholds
        H = Hres;
        idx.push_back(i); // Add marker index to the list
        idx_pos_in_idx[i] = idx.size() - 1; // Position in idx
        // We may choose to update prev_cumulative_distance here if appropriate
      }

    }

    p.increment(); // Update the progress bar
  }

  // Extract the names of the successfully phased markers
  CharacterVector mrk_names(idx.size());
  for(int i = 0; i < idx.size(); i++)
    mrk_names[i] = mrk_id[idx[i]];

  // Return the result as a list containing marker names and phase configurations
  List z = List::create(Named("marker_names") = mrk_names,
                        Named("phase_configs") = H);
  return z;
}



// [[Rcpp::export]]
List phasing_one(CharacterVector mrk_id,
                 IntegerVector dose_vec,
                 NumericMatrix S,
                 NumericMatrix InitPh,
                 bool verbose){

  int ploidy = InitPh.ncol();//ploidy
  int n = mrk_id.size();//markers to be positioned
  List H(n);
  Progress p(n, verbose); //Progress bar
  for(int i = 0; i < n; ++i){
    if (Progress::check_abort())
      return -1.0;
    int x = dose_vec[i];
    if(x == ploidy)
      H[i] = make_mat(1.0, 1, ploidy);
    else if(x == 0)
      H[i] = make_mat(0.0, 1, ploidy);
    else {
      NumericMatrix z = find_valid_permutations(InitPh, S(i,_), x);
      if(z.nrow() == 0)
        H[i] = get_all_combinations(ploidy, x);
      else{
        List Htemp(z.nrow());
        for(int k = 0; k < z.nrow(); ++k)
          Htemp[k] = rbind_cpp(InitPh, z(k, _));
        List Hres = filter_matrices(Htemp);
        NumericMatrix w(Hres.size(), ploidy);
        for(int k = 0; k < Hres.size(); k++){
          NumericMatrix y = Hres[k];
          w(k,_) = y(y.nrow() - 1,_);
        }
        H[i] = w;
      }
    }
  }
  return H;
}




