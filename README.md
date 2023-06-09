# mappoly2

Genetic Linkage Maps in Autopolyploids v 2.0

## Features to be implemented soon

1. [x] DONE ~~Rectify the C++ code to correctly exchange parents with distinct ploidy levels.~~

2. [x] DONE ~~Introduce the functionality to support multilocus analysis for individual parents.~~

3. [x] DONE ~~Prevent underflow in the forward-backward algorithm: probability normalization was employed during both the forward and backward steps. [See this code for more information.](https://github.com/mmollina/mappoly2/commit/ee4d0b8938b0631e377959d4f8f0c6fa27c0c8e7#diff-f405d1ef79df16b745f22994e5c42adddb61716567b5f0d029ce5de6c9b98cadR341)~~

4. [x] DONE ~~Employing a comprehensive approach that combines two-point phasing with a multilocus likelihood-based method.~~

5. [x] DONE ~~Implement the compute_hmm_log_likelihood function employing logarithmic computations within the forward algorithm framework (in C++).~~

6. [x] DONE ~~Integrate emission probabilities to address global genotyping error, within the multilocus approach to improve accuracy.~~ 

7. [x] ~~DONE Develop a procedure to calculate conditional probabilities of genotypes:~~
    - ~~When genotypes from both parents are informative.~~
    - ~~When only the genotype from one parent is informative.~~
    - ~~When genotypes from both parents are informative, with error consideration.~~
    - ~~When only the genotype from one parent is informative, with error consideration.~~
    - ~~Implement the R wrapper functions for C++ calc_genoprob (still to be implemented --> genoprob.R)~~
    
8. [ ] TODO Implement function add_marker, given a phased map. Use pre-computed conditional probabilities of genotypes. 
    - ~~Implement function to list possible phase configurations of a sequence of unmapped 
      markers given a phased map and a pairwise recombination fraction matrix and associated 
      phase statistics~~
    - Implement a function that receives as input a matrix whose rows represent possible phases 
      for a marker, insert it on a precomputed map and return the log-likelihood
    
9. [ ] TODO Use Hidden Markov Models (HMM) to assimilate remaining markers into the phased map.

10. [ ] TODO Enhance the efficiency of two-point calculations by harnessing the power of RcppParallel for its implementation.
