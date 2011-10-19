### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

check_type_size <- function(x, type, size, isNullOK=FALSE, check.NA=TRUE){
	if(is.null(x)){
		if(all(size == 0) || isNullOK) return(TRUE)
		else stop("The input is null");
	}
	if(type == "double"){
		if(!is.double(x)) stop("The input should be double");
	}else if(type == "integer" || type == "int"){
		if(!is.integer(x)) stop("The input should be integer");
	}else stop("Unknown type: ",type,sep="");
	d = dim(x);
	if(is.null(d)) d = length(x);
	if(length(d) != length(size) || any(d != size)) stop("Dimensionality mismatch: (",paste(d,collapse=" x "),") vs (",paste(size,collapse=" x "),")");
	if(check.NA && any(is.na(x))) stop("Some elements are NA");
}

###
### Power iteration: 
### Given A, w, prior, find eigen.vector such that
###   lambda * eigen.vector[i] = sum_j { (w[j] * A[j,i] + (1 - w[j]) * prior[i]) * eigen.vector[j] }
### by repeating this computation until convergence.
### Initially, we set eigen.vector[i] = init[i].
###
### The matrix A is specified in a sparse format:
###     edges=data.frame(from,to,value)
###         A[from,to] = value
###         "from" and "to" are indices starting from 1  << IMPORTANT!!
###
### option == -1: no normalization (run max.iter iterations)
###            0: run max.iter iterations
###            1: stop when nNodes * max_i abs(ev[i] - ev_out[i]) / sum(in) < eps (no more than max.iter)
###
power.iteration <- function(edges, init, prior, w, debug=3, option=1, max.iter=500, eps=1e-4, index.start=1){
    nEdges = as.integer(nrow(edges));
    nNodes = as.integer(length(init));
    out    = rep(double(1), nNodes);
    lambda = rep(as.double(1.0),1);
    num.iter = as.integer(rep(0,1));
	
    if(length(w) == 1) w = rep(w, nNodes);
    if(length(prior) == 1) prior = rep(prior, nNodes);
    
	# set w[j] = 0 if sum_i A[j,i] = 0
	temp = unique(edges$from);
	w[!(1:nNodes %in% temp)] = 0;
	
	if(index.start == 1){
		edges$from = as.integer(edges$from-1);
		edges$to   = as.integer(edges$to-1);
	}else if(index.start != 0) stop("index.start can only be 0 or 1");

    check_type_size(init, "double", nNodes);
    check_type_size(prior, "double", nNodes);
    check_type_size(w, "double", nNodes);
    check_type_size(edges$from, "int", nEdges);
    check_type_size(edges$to, "int", nEdges);
    check_type_size(edges$value, "double", nEdges);
    
	# The following C function is defined in src/C/pagerank.c
	.C("power_iteration",
        out,       # OUTPUT vector: nNodes x 1
        lambda,    # OUTPUT 1x1
        num.iter,  # OUTPUT 1x1
        init,      # INPUT  vector: nNodes x 1
        prior,     # INPUT  prior vector: nNodes x 1
        w,         # INPUT  teleporting probability vector: nNodes x 1
        edges$from,  # INPUT transition probability matrix: nEdges x 1 (index start from 0)
        edges$to,    # INPUT transition probability matrix: nEdges x 1 (index start from 0)
        edges$value, # INPUT transition probability matrix: nEdges x 1
        nNodes, nEdges, as.integer(max.iter), as.double(eps), as.integer(option), 
        as.integer(debug), # 0: no debugging.  1: check w.  2: positivity.  3: sum-up to one.
        DUP=FALSE
    );
    
	if(num.iter >= max.iter && option == 1){
		warning("power.iteration does not converge: maximum number of iterations reached!!");
	}
	
    return(list(eigen.vec=out, lambda=lambda, num.iter=num.iter));
}

###
### Normalize a non-negative vector so that each group sumup to one
###	(set to 0 if all values in a group are zero)
###
normalize.sumToOne.groupBy <- function(x, by){
	len = length(x);
	if(length(by) != len) stop("x and by have different length");
	check_type_size(x,  "double", len);
	check_type_size(by, "int",    len);
	out = rep(double(1), len);
	
	# The following C function is defined in src/C/util.c
	.C("normalize_sumToOne_groupby", 
			out, x, by, len,
			DUP=FALSE
	);
	return(out);
}

###
### Create the edge weights for PageRank
###
###		A_{ji} = fn(b0 + b1 * nPos_{ji} * (1 - (b5 + nPos_{ij}) / (b6 + n_{ij}) )^b3 
###                    - b2 * nNeg_{ji} * (    (b5 + nPos_{ij}) / (b6 + n_{ij}) )^b4
###              )
###		data = data.frame(voter, author, nPos, nNeg, nPos.rev, nNeg.rev);
###		fn is either 'raw', 'log', 'sigmoid' 
###
generate.edges <- function(
	data, b, fn, normalize=TRUE
){
	if(!all(c("voter","author","nPos","nNeg","nPos.rev","nNeg.rev") %in% names(data))) stop("Some columns are missing");
	if(!(length(b) %in% c(3,7))) stop("length(b) = ", length(b)," (should be either 3 or 7)");
	if(!(fn %in% c("raw", "log", "sigmoid"))) stop("fn = '",fn,"'");
	b0 = b[1];  b = b[-1];
	if(length(b) == 2){
		x = b0 + b[1] * data$nPos - b[2] * data$nNeg;
	}else{
		posRate.rev = (b[5] + data$nPos.rev) / (b[6] + data$nPos.rev + data$nNeg.rev);
		x = b0 + b[1] * data$nPos * (1 - posRate.rev)^b[3] - b[2] * data$nNeg * (posRate.rev)^b[4];
	}
	if(fn == "raw"){
		x[x < 0] = 0;
	}else if(fn == "log"){
		x[x < 0] = 0;
		x = log1p(x);
	}else if(fn == "sigmoid"){
		x = 1/(1 + exp(-x));
	}
	edges = data.frame(from=data$voter, to=data$author, value=x);
	edges = edges[edges$value!=0,];
	if(normalize) edges$value = normalize.sumToOne.groupBy(x=edges$value, by=edges$from);
	return(edges);
}

###
### output[i] = sum_j { w[j] * A[j,i] + (1 - w[j]) * prior[i] } x[j]
###
### data = list(edges=data.frame(from,to,value), prior, w);
###     A[from,to] = value
###     "from" and "to" are indices starting from 1
###
compute_one_step_transition <- function(x, data){
    debug = if(is.null(data$debug)) 0 else data$debug;
    ans = power.iteration(edges=data$edges, init=x, prior=data$prior, w=data$w, debug=debug, option=-1, max.iter=1, index.start=1);
    return(ans$eigen.vec);
}

###
### node indices start from 1
###
eigen.arpack <- function(edges, prior, w, debug=0, max.iter=500){
	edges$from = as.integer(edges$from-1);
	edges$to   = as.integer(edges$to-1);
	# set w[j] = 0 if sum_i A[j,i] = 0
	temp = unique(edges$from);
	w[!(1:length(prior) %in% temp)] = 0;
	
	data = list(edges=edges, prior=prior, w=w, debug=debug);
	options = igraph.arpack.default;
	options$n = length(prior);
	options$which = "LM";
	options$nev = 1;
	options$maxiter = max.iter;

	ans = arpack(func=compute_one_step_transition, extra=data, sym=FALSE, options=options);
	
	eigen.vec = as.real(ans$vectors);
	eigen.vec = eigen.vec / sum(eigen.vec);
	
	return(list(eigen.vec=eigen.vec, lambda=as.real(ans$values), num.iter=ans$options$iter));
}

