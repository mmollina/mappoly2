# mappoly2

Genetic Linkage Maps in Autopolyploids v 2.0

## Features to be implemented soon

1. [x] DONE ~~Rectify the C++ code to correctly exchange parents with distinct ploidy levels.~~

2. [x] DONE ~~Introduce the functionality to support multilocus analysis for individual parents.~~

3. [x] DONE ~~Prevent underflow in the forward-backward algorithm: probability normalization was employed during both the forward and backward steps. [See this code for more information.](https://github.com/mmollina/mappoly2/commit/ee4d0b8938b0631e377959d4f8f0c6fa27c0c8e7#diff-f405d1ef79df16b745f22994e5c42adddb61716567b5f0d029ce5de6c9b98cadR341)~~

4. [x] DONE ~~Employing a comprehensive approach that combines two-point phasing with a multilocus likelihood-based method.~~

5. [x] DONE ~~Implement the compute_hmm_log_likelihood function employing logarithmic computations within the forward algorithm framework (in C++).~~

6. [x] DONE ~~Integrate emission probabilities to address global genotyping error, within the multilocus approach to improve accuracy.~~ 

7. [x] DONE ~~Develop a procedure to calculate conditional probabilities of genotypes:~~
    - ~~When genotypes from both parents are informative.~~
    - ~~When only the genotype from one parent is informative.~~
    - ~~When genotypes from both parents are informative, with error consideration.~~
    - ~~When only the genotype from one parent is informative, with error consideration.~~
    - ~~Implement the R wrapper functions for C++ calc_genoprob (still to be implemented --> genoprob.R)~~
    
8. [x] DONE ~~Implement function add_marker, given a phased map. Use pre-computed conditional probabilities of genotypes~~. 
    - ~~Implement function to list possible phase configurations of a sequence of unmapped 
      markers given a phased map and a pairwise recombination fraction matrix and associated 
      phase statistics~~
    - ~~Implement the function 'homologprob_to_hmmstates': This subroutine takes in a homolog 
        probability vector (associated with both parent units) and outputs an emission vector. 
        This resultant vector serves as input for the reconstruction of flanking positions, 
        enabling the incorporation of a new marker for subsequent phasing and re-estimation 
        via three-point analysis.~~
    - ~~Implement function to visit states for one unmapped marker flanked by two mapped markers~~
    - ~~Implement a function that receives as input a matrix whose rows represent possible phases 
      for a marker, insert it on a pre-computed map and return the log-likelihood~~
9. [ ] TODO Implement two-point recombination fraction estimation using RcppParallel.    
10. [ ] TODO Implement function to re-phase each marker individually, given the rest of teh map. 
        Implement for a pre-defined segment of the map or the whole map. 
11. [ ] TODO I am updating MAPpoly to streamline its core C++ functions for easier adaptation to 
        multi-population mapping. Once this is done, I will evaluate the adaptability of these 
        updated functions for multi-population scenarios.
12. [ ] TODO Enhance the user-friendliness of the software:
    - A
    - B
    - C
    - etc
13. [ ] TODO Tutorial on how to use MAPpoly2.0 to build maps in diploid families derived from inbred lines (BC, F2, RILs, etc)
