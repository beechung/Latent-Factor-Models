### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### Compile C code:
### R CMD SHLIB src/R/examples/optim.c -o lib/my_optim.so
###

dyn.load("lib/my_optim.so");

x = c( 0, 0, 0, 0);
w = c(-1, 2, 3, 4);
lower = c(1,1,1,1);
upper = c(3,3,3,3);

.C("min_sum_of_squares",
	as.double(x),
	as.double(w),
	as.double(lower),
	as.double(upper),
	length(w),
	as.integer(1),
	as.integer(1)
);
