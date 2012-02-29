
#	Copyright (c) 2012, Yahoo! Inc.  All rights reserved.	
#   Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
#   Author: Deepak Agarwal

all:c_funcs

srcfiles = arms.c arsspline.c



c_funcs: 
	R CMD SHLIB $(srcfiles) -o c_funcs.so

.PHONY : clean
clean:
	rm *.so *.o