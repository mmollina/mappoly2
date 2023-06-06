---
editor_options: 
  markdown: 
    wrap: 72
---

# mappoly2

Genetic Linkage Maps in Autopolyploids v 2.0

## Features to be implemented soon

1.  [x] DONE ~~Rectify the C++ code to correctly exchange parents with
    distinct ploidy levels.~~

2.  [x] DONE ~~Introduce the functionality to support multilocus
    analysis for individual parents.~~

3.  [x] DONE ~~Prevent underflow in the forward-backward algorithm:
    probability normalization was employed during both the forward and
    backward steps. See [this code for more
    information](%22https://github.com/mmollina/mappoly2/commit/ee4d0b8938b0631e377959d4f8f0c6fa27c0c8e7#diff-f405d1ef79df16b745f22994e5c42adddb61716567b5f0d029ce5de6c9b98cadR341%22)~~

4.  [x] DONE ~~Employing a comprehensive approach that combines
    two-point phasing with a multilocus likelihood-based method.~~

5.  [x] DONE ~~Craft the compute_hmm_log_likelihood function employing
    logarithmic computations within the forward algorithm framework (in
    C++)~~.

6.  [x] DONE ~~Integrate emission probabilities to address global
    genotyping error, within the multilocus approach to improve
    accuracy.~~ 

7.  [ ] TODO Deploy Hidden Markov Models (HMM) to assimilate remaining
    markers into the phased map.

8.  [ ] TODO Enhance the efficiency of two-point calculations by
    harnessing the power of RcppParallel for its implementation.
