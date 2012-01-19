This directory contains code for fitting the localized factor model described
in the following paper:

       Deepak Agarwal, Bee-Chung Chen, Bo Long. Localized factor models for 
       multi-context recommendation. KDD 2011.

This model jointly factorizes multiple matrices to provide transfer learning
across multiple contexts or applications.

(0) See Notation.txt for the model first.

(1) To compile, just type make in this directory.
    see Makefile for details.
    
(2) Some C and R functions used are in
    ../R/utils.{c,h,hpp,R}

(3) Please DO NOT check in .o and .so files into the repository

Directories:
  C/           C and C++ code
  R/example    Examples of how to run different functions
  R/model      Code for modeling
  