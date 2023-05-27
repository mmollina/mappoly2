/*
 MAPpoly: a package to construct genetic maps in autopolyploids
 Copyright (C) 2014-2022 Marcelo Mollinari

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
 File: hmm_elements.h

 Description:

 Functions Written by Marcelo Mollinari.

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version:       2022
 Last update: Feb 21, 2023
 */

double prob_k1_given_k_l_m(int ploidy, int l, double rf);

double log_prob_k1_given_k_l_m(int ploidy, int l, double rf);

std::vector<std::vector<double> > rec_num(int ploidy);

std::vector<std::vector<double> > transition(int ploidy, double rf);

std::vector<std::vector<double> > log_transition(int ploidy, double rf);

std::vector<double> forward_emit(std::vector<double>& fk,
                                 std::vector<int>& ik,
                                 std::vector<int>& ik1,
                                 std::vector<double>& emit,
                                 std::vector<std::vector<double> >& T1,
                                 std::vector<std::vector<double> >& T2);

std::vector<long double> forward_emit_highprec(std::vector<long double>& fk,
                                               std::vector<int>& ik,
                                               std::vector<int>& ik1,
                                               std::vector<double>& emit,
                                               std::vector<std::vector<double> >& T1,
                                               std::vector<std::vector<double> >& T2);

std::vector<double> log_forward_emit(std::vector<double>& fk,
                                     std::vector<int>& ik,
                                     std::vector<int>& ik1,
                                     std::vector<double>& emit,
                                     std::vector<std::vector<double> >& T1,
                                     std::vector<std::vector<double> >& T2);

std::vector<double> backward_emit(std::vector<double>& fk1,
                                  std::vector<int>& ik,
                                  std::vector<int>& ik1,
                                  std::vector<double>& emit,
                                  std::vector<std::vector<double> >& T1,
                                  std::vector<std::vector<double> >& T2);

std::vector<long double> backward_emit_highprec(std::vector<long double>& fk1,
                                                std::vector<int>& ik,
                                                std::vector<int>& ik1,
                                                std::vector<double>& emit,
                                                std::vector<std::vector<double> >& T1,
                                                std::vector<std::vector<double> >& T2);

std::vector<double> log_backward_emit(std::vector<double>& fk1,
                                      std::vector<int>& ik,
                                      std::vector<int>& ik1,
                                      std::vector<double>& emit,
                                      std::vector<std::vector<double> >& T1,
                                      std::vector<std::vector<double> >& T2);

double addlog(double a, double b);
