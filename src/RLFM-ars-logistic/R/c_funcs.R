### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
###
### Author: Liang Zhang

###
### Inversion of a symmetric metrix using Cholesky
###     This function is implemented to make sure the C implementation generates exactly
###     the same output as the R implementation
###
sym_inv.cholesky <- function(x){
    if(!is.double(x)) stop("x is not a double matrix");
    if(!is.matrix(x)) stop("x is not a double matrix");
    if(nrow(x) != ncol(x)) stop("x is not a square matrix");
    n = as.integer(nrow(x));
    ans = .C("sym_inv_byCholesky", x, n, as.integer(1), DUP=TRUE);
    return(ans[[1]]);
}

###
### Eigen decomposition of a symmetric metrix
###     This function is implemented to make sure the C implementation generates exactly
###     the same output as the R implementation
###
sym_eigen <- function(x){
    if(!is.matrix(x)) stop("x is not matrix");
    if(nrow(x) != ncol(x)) stop("nrow(x) != ncol(x)");
    if(!is.double(x)) stop("x is not a double-precision matrix");
    output = list(values  = rep(0.0, nrow(x)),
                  vectors = matrix(0.0, nrow=nrow(x), ncol=ncol(x)));
    .C("sym_eigen", x, nrow(x), output$values, output$vectors, DUP=FALSE);
    return(output);
}

###
### Sum up each row (or column) of a matrix
###     side=1: Sum up each row and return a vector with length nrow
###     side=2: Sum up each column and return a vector with length ncol
###
sum_margin <- function(A, side){
    if(!is.matrix(A)) stop("'A' is not a matrix");
    if(!is.double(A)) stop("'A' is not double");
    if(side == 1) out = double(nrow(A))
    else if(side == 2) out = double(ncol(A))
    else stop("Unknown side=",side," (side=1 for rows, side=2 for columns)");
    ans = .C("sum_margin", out, A, as.integer(nrow(A)), as.integer(ncol(A)), as.integer(side), DUP=FALSE);
    return(out);
}

###
### Draw a random sample from a multivariate normal distribtuion
###     This function is implemented to make sure the C implementation generates exactly
###     the same output as the R implementation
###
my_rmvnorm <- function(n, mu, Sigma, debug=10, tol=1e-8){
    p <- length(mu)
    if (!all(dim(Sigma) == c(p, p))) 
        stop("incompatible arguments")
    eS <- sym_eigen(solve(Sigma))
    ev <- 1/eS$values
    
    if(debug >= 3){
        if(max(abs(eS$vectors %*% diag(ev, p) %*% t(eS$vectors) - Sigma)) > tol * abs(ev[1]))
            stop("sym_eigen(Sigma) seems to have some problems!!");
    }
    
    if (!all(ev >= -tol * abs(ev[1]))) 
        stop("'Sigma' is not positive definite")
    Z <- matrix(rnorm(p * n), n)
    
    # cat("temp:\n");
    # print(eS$vectors %*% diag(sqrt(pmax(ev, 0)), p));
    # cat("rnd: ");
    # print(drop(Z));
    
    X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(Z)
    nm <- names(mu)
    if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) 
        nm <- dn[[1]]
    dimnames(X) <- list(nm, NULL)
    if (n == 1) 
        return(drop(X))
    else return(t(X))
}

###
### The E-Step for the Monte-Carlo EM algorithm (version 1)
###
MC_EStep.C <- function(
    nSamples,
    user, item, y, x, w, z,
    alpha, beta, u, v,
    b, g0, G, d0, D,
    var_y, var_alpha, var_beta, var_u, var_v=1, debug=0
){
    if(debug >= 1) check.input(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_y, var_alpha, var_beta, var_u, var_v);
    nObs     = as.integer(length(y));
    nUsers   = as.integer(length(alpha));
    nItems   = as.integer(length(beta));
    nFactors = as.integer(ncol(u));

    xb  = x %*% b;
    g0w = w %*% g0;
    d0z = z %*% d0;
    Gw  = w %*% G;
    Dz  = z %*% D;

    output = list(
        o.mean       = rep(0.0, nObs),
        o.sumvar     = 0.0,
        alpha.mean   = rep(0.0, nUsers),
        alpha.sumvar = 0.0,
        beta.mean    = rep(0.0, nItems),
        beta.sumvar  = 0.0,
        u.mean       = matrix(0.0, nrow=nUsers, ncol=nFactors),
        u.sumvar     = 0.0,
        v.mean       = matrix(0.0, nrow=nItems, ncol=nFactors),
        v.sumvar     = 0.0,
		pred.y.square = rep(0.0, nObs)
    );
    
    if(!is.integer(user)) stop("!is.integer(user)");
    if(!is.integer(item)) stop("!is.integer(item)");
    if(!is.double(y)) stop("!is.double(y)");
    if(!is.double(xb)) stop("!is.double(xb)");
    if(!is.double(g0w)) stop("!is.double(g0w)");
    if(!is.double(d0z)) stop("!is.double(d0z)");
    if(!is.double(Gw)) stop("!is.double(Gw)");
    if(!is.double(Dz)) stop("!is.double(Dz)");
    if(!is.double(alpha)) stop("!is.double(alpha)");
    if(!is.double(beta)) stop("!is.double(beta)");
    if(!is.double(u)) stop("!is.double(u)");
    if(!is.double(v)) stop("!is.double(v)");
    if(!is.double(var_y)) stop("!is.double(var_y)");
    if(!is.double(var_alpha)) stop("!is.double(var_alpha)");
    if(!is.double(var_beta)) stop("!is.double(var_beta)");
    if(!is.double(var_u)) stop("!is.double(var_u)");
    if(!is.double(var_v)) stop("!is.double(var_v)");

    .C("MCEM_EStep",
        # OUTPUT
        output$o.mean,     output$o.sumvar,
        output$alpha.mean, output$alpha.sumvar,
        output$beta.mean,  output$beta.sumvar,
        output$u.mean,     output$u.sumvar,
        output$v.mean,     output$v.sumvar,
		output$pred.y.square,
        # INPUT
        as.integer(nSamples),
        user, item, y, xb, g0w, d0z, Gw, Dz, 
        alpha, beta, u, v, 
        var_y, var_alpha, var_beta, var_u, var_v,
        nObs, nUsers, nItems, nFactors, as.integer(length(var_y)),
        # OTHER
        as.integer(debug),
        DUP=FALSE
    );
    return(output);
}

###
### The E-Step for the Monte-Carlo EM algorithm (version 1)
###
MC_EStep_logistic.C <- function(
    nSamples, nBurnin,
    user, item, y, x, w, z,
    alpha, beta, u, v,
    b, g0, G, d0, D,
    var_alpha, var_beta, var_u, var_v=1, debug=0
){
    if(debug >= 1) check.input.logistic(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_alpha, var_beta, var_u, var_v);
    nObs     = as.integer(length(y));
    nUsers   = as.integer(length(alpha));
    nItems   = as.integer(length(beta));
    nFactors = as.integer(ncol(u));

    xb  = x %*% b;
    g0w = w %*% g0;
    d0z = z %*% d0;
    Gw  = w %*% G;
    Dz  = z %*% D;

    output = list(
    	o.mean       = rep(0.0, nObs),
        alpha.mean   = rep(0.0, nUsers),
        alpha.sumvar = 0.0,
        beta.mean    = rep(0.0, nItems),
        beta.sumvar  = 0.0,
        u.mean       = matrix(0.0, nrow=nUsers, ncol=nFactors),
        u.sumvar     = 0.0,
        v.mean       = matrix(0.0, nrow=nItems, ncol=nFactors),
        v.sumvar     = 0.0,
	acceptrate.maineff = 0.0,
	acceptrate.fact = 0.0
    );

    if(!is.integer(user)) stop("!is.integer(user)");
    if(!is.integer(item)) stop("!is.integer(item)");
    if(!is.double(y)) stop("!is.double(y)");
    if(!is.double(xb)) stop("!is.double(xb)");
    if(!is.double(g0w)) stop("!is.double(g0w)");
    if(!is.double(d0z)) stop("!is.double(d0z)");
    if(!is.double(Gw)) stop("!is.double(Gw)");
    if(!is.double(Dz)) stop("!is.double(Dz)");
    if(!is.double(alpha)) stop("!is.double(alpha)");
    if(!is.double(beta)) stop("!is.double(beta)");
    if(!is.double(u)) stop("!is.double(u)");
    if(!is.double(v)) stop("!is.double(v)");
    if(!is.double(var_alpha)) stop("!is.double(var_alpha)");
    if(!is.double(var_beta)) stop("!is.double(var_beta)");
    if(!is.double(var_u)) stop("!is.double(var_u)");
    if(!is.double(var_v)) stop("!is.double(var_v)");

    .C("MCEM_EStep_logistic",
        # OUTPUT
	output$o.mean,
        output$alpha.mean, output$alpha.sumvar,
        output$beta.mean,  output$beta.sumvar,
        output$u.mean,     output$u.sumvar,
        output$v.mean,     output$v.sumvar,
	output$acceptrate.maineff, output$acceptrate.fact,
        # INPUT
        as.integer(nSamples), as.integer(nBurnin),
        user, item, y, xb, g0w, d0z, Gw, Dz,
        alpha, beta, u, v,
        var_alpha, var_beta, var_u, var_v,
        nObs, nUsers, nItems, nFactors,
        # OTHER
        as.integer(debug),
        DUP=FALSE
    );
    return(output);
}

###
### The E-Step for the Monte-Carlo EM algorithm (version 2)
###
### Options:
###
###   outputPerUserVar: [TRUE/FALSE] Whether to output per-user variance-covariance matrix
###   outputPerItemVar: [TRUE/FALSE] Whether to output per-item variance-covariance matrix
###
###   isOldUser[i]: [TRUE/FALSE] Whehter the ith user is an old user; default: all FALSE
###                 For old users, the input alpha[i], u[i,] will be used instead of g0w[i] and Gw[i]
###   isOldItem[j]: [TRUE/FALSE] Whehter the jth item is an old item; default: all FALSE
###                 For old items, the input beta[j], v[j,] will be used instead of d0z[j] and Dz[j]
###
###   var_y     (1x1 or nObs x 1) specifies global or per-observation variance
###   var_alpha (1x1 or nUsers x 1) specifies global or per-user variance
###   var_beta  (1x1 or nItems x 1) specifies global or per_item variance
###   var_u     (1x1 or nUsers x nFactors x nFactors) specifies global variance or per-user variance-covariance matrix
###   var_v     (1x1 or nItems x nFactors x nFactors) specifies global variance or per-item variance-covariance matrix
###
### Output:
###
###   output$o.mean[k] = E[ o[k] = alpha[user[k]] + beta[item[k]] + t(u[user[k]]) %*% v[item[k]] ]
###   output$o.sumvar  = sum_k Var(o[k])
###
###   output$alpha.mean[i] = E[ alpha[i] ]
###   output$alpha.sumvar  = sum_i Var( alpha[i] )
###   output$alpha.var[i]  = Var( alpha[i] ), which is available only when outputPerUserVar=TRUE
###
###   output$beta.mean[i] = E[ beta[i] ]
###   output$beta.sumvar  = sum_i Var( beta[i] )
###   output$beta.var[i]  = Var( beta[i] ), which is available only when outputPerItemVar=TRUE
###
###   output$u.mean[i,] = E[ u[i] ];  u.mean is nUsers x nFactors
###   output$u.sumvar   = sum_i Var( u[i] )
###   output$u.var[i,,] = VarCov( u[i] ), which is available only when outputPerUserVar=TRUE
###                       u.var is nUsers x nFactors x nFactors
###
###   output$v.mean[i,] = E[ v[i] ];  v.mean is nItems x nFactors
###   output$v.sumvar   = sum_i Var( v[i] )
###   output$v.var[i,,] = VarCov( v[i] ), which is available only when outputPerItemVar=TRUE
###                       v.var is nItems x nFactors x nFactors
###
MC_EStep2.C <- function(
    nSamples,
    user, item, y, x, w, z,
    alpha, beta, u, v,
    b, g0, G, d0, D,
    var_y, var_alpha, var_beta, var_u, var_v=1,
    outputPerUserVar=FALSE,
    outputPerItemVar=FALSE,
    isOldUser=NULL,
    isOldItem=NULL,
    debug=0
){
    if(debug >= 1) check.input(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_y, var_alpha, var_beta, var_u, var_v, version=2);
    nObs     = as.integer(length(y));
    nUsers   = as.integer(length(alpha));
    nItems   = as.integer(length(beta));
    nFactors = as.integer(ncol(u));

    nVar_y     = as.integer(length(var_y));
    nVar_alpha = as.integer(length(var_alpha));
    nVar_beta  = as.integer(length(var_beta));
    nVar_u     = as.integer(1); if(length(var_u) > 1) nVar_u = as.integer(dim(var_u)[1]);
    nVar_v     = as.integer(1); if(length(var_v) > 1) nVar_v = as.integer(dim(var_v)[1]);
    
    xb  = x %*% b;
    if(is.null(isOldUser)){
        g0w = w %*% g0;
        Gw  = w %*% G;
    }else{
        if(length(isOldUser) != nUsers) stop("length(isOldUser) != nUsers");
        w.new = w[!isOldUser,,drop=FALSE];
        g0w = alpha;
        Gw  = u;
        g0w[!isOldUser] = w.new %*% g0;
        Gw[!isOldUser,] = w.new %*% G;
    }
    
    if(is.null(isOldItem)){
        d0z = z %*% d0;
        Dz  = z %*% D;
    }else{
        if(length(isOldItem) != nItems) stop("length(isOldItem) != nItems");
        z.new = z[!isOldItem,,drop=FALSE];
        d0z = beta;
        Dz  = v;
        d0z[!isOldItem] = z.new %*% d0;
        Dz[!isOldItem,] = z.new %*% D;
    }

    output = list(
        o.mean       = rep(0.0, nObs),
        o.sumvar     = 0.0,
        alpha.mean   = rep(0.0, nUsers),
        alpha.sumvar = 0.0,
        beta.mean    = rep(0.0, nItems),
        beta.sumvar  = 0.0,
        u.mean       = matrix(0.0, nrow=nUsers, ncol=nFactors),
        u.sumvar     = 0.0,
        v.mean       = matrix(0.0, nrow=nItems, ncol=nFactors),
        v.sumvar     = 0.0
    );
    
    if(outputPerUserVar){
        perUserVar = as.integer(1);
        output$alpha.var = rep(0.0, nUsers);
        output$u.var     = array(0.0, dim=c(nUsers, nFactors, nFactors));
    }else{
        perUserVar = as.integer(0);
    }
    
    if(outputPerItemVar){
        perItemVar = as.integer(1);
        output$beta.var = rep(0.0, nItems);
        output$v.var    = array(0.0, dim=c(nItems, nFactors, nFactors));
    }else{
        perItemVar = as.integer(0);
    }
    
    if(!is.integer(user)) stop("!is.integer(user)");
    if(!is.integer(item)) stop("!is.integer(item)");
    if(!is.double(y)) stop("!is.double(y)");
    if(!is.double(xb)) stop("!is.double(xb)");
    if(!is.double(g0w)) stop("!is.double(g0w)");
    if(!is.double(d0z)) stop("!is.double(d0z)");
    if(!is.double(Gw)) stop("!is.double(Gw)");
    if(!is.double(Dz)) stop("!is.double(Dz)");
    if(!is.double(alpha)) stop("!is.double(alpha)");
    if(!is.double(beta)) stop("!is.double(beta)");
    if(!is.double(u)) stop("!is.double(u)");
    if(!is.double(v)) stop("!is.double(v)");
    if(!is.double(var_y)) stop("!is.double(var_y)");
    if(!is.double(var_alpha)) stop("!is.double(var_alpha)");
    if(!is.double(var_beta)) stop("!is.double(var_beta)");
    if(!is.double(var_u)) stop("!is.double(var_u)");
    if(!is.double(var_v)) stop("!is.double(var_v)");

    .C("MCEM_EStep2",
        # OUTPUT
        output$o.mean,     output$o.sumvar,     
        output$alpha.mean, output$alpha.sumvar, output$alpha.var,
        output$beta.mean,  output$beta.sumvar,  output$beta.var,
        output$u.mean,     output$u.sumvar,     output$u.var,
        output$v.mean,     output$v.sumvar,     output$v.var,
        # INPUT
        as.integer(nSamples),
        user, item, y, xb, g0w, d0z, Gw, Dz, 
        alpha, beta, u, v, 
        var_y, var_alpha, var_beta, var_u, var_v,
        nObs, nUsers, nItems, nFactors,
        nVar_y, nVar_alpha, nVar_beta, nVar_u, nVar_v,
        perUserVar, perItemVar,
        # OTHER
        as.integer(debug),
        DUP=FALSE
    );
    return(output);
}


###
### The gradient descent step of the ICM method
###
GradDesc.C <- function(
    nSteps, user, item, y, x, w, z,
    alpha, beta, u, v, b, g0, G, d0, D,
    var_y, var_alpha, var_beta, var_u, var_v, gd.method=3, eta=NULL, tol=1e-8, 
    debug=0, verbose=0
){
    trace = 0;
    if(verbose >= 6) trace = verbose - 5;
    
    if(gd.method != 3) stop("GradDesc.C only supports gd.method=3 (Conjugate Gradient Descend)");
    
    abstol = as.double(0);
    type   = as.integer(1);
    trace  = as.integer(trace);
    
    if(debug >= 1) check.input(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_y, var_alpha, var_beta, var_u, var_v);
    nObs     = as.integer(length(y));
    nUsers   = as.integer(length(alpha));
    nItems   = as.integer(length(beta));
    nFactors = as.integer(ncol(u));

    xb  = x %*% b;
    g0w = w %*% g0;
    d0z = z %*% d0;
    Gw  = w %*% G;
    Dz  = z %*% D;

    output = list(
        alpha = rep(0.0, nUsers),
        beta  = rep(0.0, nItems),
        u     = matrix(0.0, nrow=nUsers, ncol=nFactors),
        v     = matrix(0.0, nrow=nItems, ncol=nFactors)
    );
    
    if(!is.integer(user)) stop("!is.integer(user)");
    if(!is.integer(item)) stop("!is.integer(item)");
    if(!is.double(y)) stop("!is.double(y)");
    if(!is.double(xb)) stop("!is.double(xb)");
    if(!is.double(g0w)) stop("!is.double(g0w)");
    if(!is.double(d0z)) stop("!is.double(d0z)");
    if(!is.double(Gw)) stop("!is.double(Gw)");
    if(!is.double(Dz)) stop("!is.double(Dz)");
    if(!is.double(alpha)) stop("!is.double(alpha)");
    if(!is.double(beta)) stop("!is.double(beta)");
    if(!is.double(u)) stop("!is.double(u)");
    if(!is.double(v)) stop("!is.double(v)");
    if(!is.double(var_y)) stop("!is.double(var_y)");
    if(!is.double(var_alpha)) stop("!is.double(var_alpha)");
    if(!is.double(var_beta)) stop("!is.double(var_beta)");
    if(!is.double(var_u)) stop("!is.double(var_u)");
    if(!is.double(var_v)) stop("!is.double(var_v)");
    
    ans = .C("ICM_CGD",
        # OUTPUT
        output$alpha, output$beta, output$u, output$v,
        # INPUT
        as.integer(nSteps), abstol, as.double(tol), type, trace,
        user, item, y, xb, g0w, d0z, Gw, Dz, 
        alpha, beta, u, v, 
        var_y, var_alpha, var_beta, var_u, var_v,
        nObs, nUsers, nItems, nFactors,
        # OTHER
        as.integer(debug), as.integer(verbose),
        DUP=FALSE
    );
    return(output);
}

###
### The gradient descent step of the ICM method (version 2)
###
GradDesc2.C <- function(
    nSteps, user, item, y, x, w, z,
    alpha, beta, u, v, b, g0, G, d0, D,
    var_y, var_alpha, var_beta, var_u, var_v, gd.method=3, eta=NULL, tol=1e-8, 
    isOldUser=NULL,
    isOldItem=NULL,
    debug=0, verbose=0
){
    trace = 0;
    if(verbose >= 6) trace = verbose - 5;
    
    if(gd.method != 3) stop("GradDesc.C only supports gd.method=3 (Conjugate Gradient Descend)");
    
    abstol = as.double(0);
    type   = as.integer(1);
    trace  = as.integer(trace);

    if(debug >= 1) check.input(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_y, var_alpha, var_beta, var_u, var_v, version=2);
    nObs     = as.integer(length(y));
    nUsers   = as.integer(length(alpha));
    nItems   = as.integer(length(beta));
    nFactors = as.integer(ncol(u));

    nVar_y     = as.integer(length(var_y));
    nVar_alpha = as.integer(length(var_alpha));
    nVar_beta  = as.integer(length(var_beta));
    nVar_u     = as.integer(1); if(length(var_u) > 1) nVar_u = as.integer(dim(var_u)[1]);
    nVar_v     = as.integer(1); if(length(var_v) > 1) nVar_v = as.integer(dim(var_v)[1]);
    
    xb  = x %*% b;
    if(is.null(isOldUser)){
        g0w = w %*% g0;
        Gw  = w %*% G;
    }else{
        if(length(isOldUser) != nUsers) stop("length(isOldUser) != nUsers");
        w.new = w[!isOldUser,,drop=FALSE];
        g0w = alpha;
        Gw  = u;
        g0w[!isOldUser] = w.new %*% g0;
        Gw[!isOldUser,] = w.new %*% G;
    }
    
    if(is.null(isOldItem)){
        d0z = z %*% d0;
        Dz  = z %*% D;
    }else{
        if(length(isOldItem) != nItems) stop("length(isOldItem) != nItems");
        z.new = z[!isOldItem,,drop=FALSE];
        d0z = beta;
        Dz  = v;
        d0z[!isOldItem] = z.new %*% d0;
        Dz[!isOldItem,] = z.new %*% D;
    }

    output = list(
        alpha = rep(0.0, nUsers),
        beta  = rep(0.0, nItems),
        u     = matrix(0.0, nrow=nUsers, ncol=nFactors),
        v     = matrix(0.0, nrow=nItems, ncol=nFactors)
    );
    
    if(!is.integer(user)) stop("!is.integer(user)");
    if(!is.integer(item)) stop("!is.integer(item)");
    if(!is.double(y)) stop("!is.double(y)");
    if(!is.double(xb)) stop("!is.double(xb)");
    if(!is.double(g0w)) stop("!is.double(g0w)");
    if(!is.double(d0z)) stop("!is.double(d0z)");
    if(!is.double(Gw)) stop("!is.double(Gw)");
    if(!is.double(Dz)) stop("!is.double(Dz)");
    if(!is.double(alpha)) stop("!is.double(alpha)");
    if(!is.double(beta)) stop("!is.double(beta)");
    if(!is.double(u)) stop("!is.double(u)");
    if(!is.double(v)) stop("!is.double(v)");
    if(!is.double(var_y)) stop("!is.double(var_y)");
    if(!is.double(var_alpha)) stop("!is.double(var_alpha)");
    if(!is.double(var_beta)) stop("!is.double(var_beta)");
    if(!is.double(var_u)) stop("!is.double(var_u)");
    if(!is.double(var_v)) stop("!is.double(var_v)");
    
    ans = .C("ICM_CGD2",
        # OUTPUT
        output$alpha, output$beta, output$u, output$v,
        # INPUT (version 1)
        as.integer(nSteps), abstol, as.double(tol), type, trace,
        user, item, y, xb, g0w, d0z, Gw, Dz, 
        alpha, beta, u, v, 
        var_y, var_alpha, var_beta, var_u, var_v,
        nObs, nUsers, nItems, nFactors,
        # INPUT (version 2)
        nVar_y, nVar_alpha, nVar_beta, nVar_u, nVar_v,
        # OTHER
        as.integer(debug), as.integer(verbose),
        DUP=FALSE
    );
    return(output);
}

###
### Fit the following model by minimizing the squared error + L2 penalty
###
###   MODEL: y(i,j) = x(i,j)' b  +  w(i)' g0  +  z(j)' d0  +  w(i)' G D' z(j) + error(i,j)
###          y(i,j): Rating of user i about item j
###          x(i,j): Feature vector for (i,j)
###          w(i)  : Feature vector for user i
###          z(j)  : Feature vector for item j
###          b, g0, d0, G, D: Model parameters
###
###   INPUT: y: #Obs x 1,               user: #Obs x 1,            item: #Obs x 1
###          x: #Obs x #JointFeatures,  w: #User x #UserFeatures,  z: #Items x #ItemFeatures
###          (y[k], x[k,], user[k], item[k]) specifies the kth observation.
###          y[k] is the rating on item[k] given by user[k] with feature vector x[k,].
###          Let i = user[k] and j = item[k].
###          The, w[i,] is the user's feature vector and
###               z[j,] is the item's feature vector.
###
###   INITIAL PARAMETER VALUES:
###          b: #JointFeatures x 1,  g0: #UserFeatures x 1,  d0: #ItemFeatures x 1
###          G: #UserFeatures x #Factors,  D: #ItemFeatures x #Factors
###
###   L2 Regularization weights:
###          lambda_b: #JointFeatures x 1,  lambda_g0: #UserFeatures x 1,  lambda_d0: #ItemFeatures x 1
###          lambda_G: #UserFeatures x #Factors,  lambda_D: #ItemFeatures x #Factors
###
###   OUTPUT: b, g0, d0, G, D that minimize
###           sum(error^2) + sum(lambda_b * b^2) + sum(lambda_g0 * g0^2) + sum(lambda_d0 * d0^2)
###                        + sum(lambda_G * G^2) + sum(lambda_D * D^2)
###
leastSquareFit_withL2.C <- function(
    nSteps, user, item, y, x, w, z,
    b, g0, G, d0, D,
    lambda_b=NULL, lambda_g0=NULL, lambda_G=NULL, lambda_d0=NULL, lambda_D=NULL, lambda=0.0,
    tol=1e-8, debug=0, verbose=0
){
    if(is.null(lambda_b))  lambda_b  = rep(lambda, length(b));
    if(is.null(lambda_g0)) lambda_g0 = rep(lambda, length(g0));
    if(is.null(lambda_d0)) lambda_d0 = rep(lambda, length(d0));
    if(is.null(lambda_G))  lambda_G  = matrix(lambda, nrow=nrow(G), ncol=ncol(G));
    if(is.null(lambda_D))  lambda_D  = matrix(lambda, nrow=nrow(D), ncol=ncol(D));

    output = list(
        b  = rep(b,1),
        g0 = rep(g0,1),
        d0 = rep(d0,1),
        G  = matrix(rep(G,1), nrow=nrow(G), ncol=ncol(G)),
        D  = matrix(rep(D,1), nrow=nrow(D), ncol=ncol(D))
    );

    if(debug >= 1) check.input2(user, item, y, x, w, z, output$b, output$g0, output$G, output$d0, output$D, lambda_b, lambda_g0, lambda_G, lambda_d0, lambda_D);

    trace = 0;
    if(verbose >= 6) trace = verbose - 5;
    
    abstol = as.double(0);
    type   = as.integer(1);
    trace  = as.integer(trace);
    
    nObs     = as.integer(length(y));
    nUsers   = as.integer(nrow(w));
    nItems   = as.integer(nrow(z));
    nFactors = as.integer(ncol(G));
    nJointFeatures = as.integer(ncol(x));
    nUserFeatures = as.integer(ncol(w));
    nItemFeatures = as.integer(ncol(z));
    
    if(!is.integer(user)) user = as.integer(user);
    if(!is.integer(item)) item = as.integer(item);
    if(!is.double(y)) stop("!is.double(y)");
    if(!is.double(x)) stop("!is.double(x)");
    if(!is.double(w)) stop("!is.double(w)");
    if(!is.double(z)) stop("!is.double(z)");
    if(!is.double(output$b)) stop("!is.double(b)");
    if(!is.double(output$g0)) stop("!is.double(g0)");
    if(!is.double(output$d0)) stop("!is.double(d0)");
    if(!is.double(output$G)) stop("!is.double(G)");
    if(!is.double(output$D)) stop("!is.double(D)");
    if(!is.double(lambda_b)) stop("!is.double(lambda_b)");
    if(!is.double(lambda_g0)) stop("!is.double(lambda_g0)");
    if(!is.double(lambda_d0)) stop("!is.double(lambda_d0)");
    if(!is.double(lambda_G)) stop("!is.double(lambda_G)");
    if(!is.double(lambda_D)) stop("!is.double(lambda_D)");

    ans = .C("min_sqrLossWithL2",
        # INPUT/OUTPUT
        output$b, output$g0, output$d0, output$G, output$D,
        # INPUT
        as.integer(nSteps), abstol, as.double(tol), type, trace,
        user, item, y, x, w, z,
        lambda_b, lambda_g0, lambda_d0, lambda_G, lambda_D,
        nObs, nUsers, nItems, nFactors, nJointFeatures, nUserFeatures, nItemFeatures,
        # OTHER
        as.integer(debug), as.integer(verbose),
        DUP=FALSE
    );
    return(output);
}

sqrLossWithL2.C <- function(
    user, item, y, x, w, z,
    b, g0, d0, G, D,
    lambda_b, lambda_g0, lambda_d0, lambda_G, lambda_D,
    debug=0
){
    if(debug >= 1) check.input2(user, item, y, x, w, z, b, g0, G, d0, D, lambda_b, lambda_g0, lambda_G, lambda_d0, lambda_D);
    output = 0.0;

    nObs     = as.integer(length(y));
    nUsers   = as.integer(nrow(w));
    nItems   = as.integer(nrow(z));
    nFactors = as.integer(ncol(G));
    nJointFeatures = as.integer(ncol(x));
    nUserFeatures = as.integer(ncol(w));
    nItemFeatures = as.integer(ncol(z));
    
    if(!is.integer(user)) user = as.integer(user);
    if(!is.integer(item)) item = as.integer(item);
    if(!is.double(y)) stop("!is.double(y)");
    if(!is.double(x)) stop("!is.double(x)");
    if(!is.double(w)) stop("!is.double(w)");
    if(!is.double(z)) stop("!is.double(z)");
    if(!is.double(b)) stop("!is.double(b)");
    if(!is.double(g0)) stop("!is.double(g0)");
    if(!is.double(d0)) stop("!is.double(d0)");
    if(!is.double(G)) stop("!is.double(G)");
    if(!is.double(D)) stop("!is.double(D)");
    if(!is.double(lambda_b)) stop("!is.double(lambda_b)");
    if(!is.double(lambda_g0)) stop("!is.double(lambda_g0)");
    if(!is.double(lambda_d0)) stop("!is.double(lambda_d0)");
    if(!is.double(lambda_G)) stop("!is.double(lambda_G)");
    if(!is.double(lambda_D)) stop("!is.double(lambda_D)");

    ans = .C("sqrLossWithL2",
        # OUTPUT
        output,
        # INPUT
        user, item, y, x, w, z,
        b, g0, d0, G, D,
        lambda_b, lambda_g0, lambda_d0, lambda_G, lambda_D,
        nObs, nUsers, nItems, nFactors, nJointFeatures, nUserFeatures, nItemFeatures,
        # OTHER
        as.integer(debug),
        DUP=FALSE    
    );
    return(output);
}

gr_sqrLossWithL2.C <- function(
    user, item, y, x, w, z,
    b, g0, d0, G, D,
    lambda_b, lambda_g0, lambda_d0, lambda_G, lambda_D,
    debug=0
){
    if(debug >= 1) check.input2(user, item, y, x, w, z, b, g0, G, d0, D, lambda_b, lambda_g0, lambda_G, lambda_d0, lambda_D);

    output = list(
        gr_b  = rep(b,1),
        gr_g0 = rep(g0,1),
        gr_d0 = rep(d0,1),
        gr_G  = matrix(rep(G,1), nrow=nrow(G), ncol=ncol(G)),
        gr_D  = matrix(rep(D,1), nrow=nrow(D), ncol=ncol(D))
    );

    nObs     = as.integer(length(y));
    nUsers   = as.integer(nrow(w));
    nItems   = as.integer(nrow(z));
    nFactors = as.integer(ncol(G));
    nJointFeatures = as.integer(ncol(x));
    nUserFeatures = as.integer(ncol(w));
    nItemFeatures = as.integer(ncol(z));
    
    if(!is.integer(user)) user = as.integer(user);
    if(!is.integer(item)) item = as.integer(item);
    if(!is.double(y)) stop("!is.double(y)");
    if(!is.double(x)) stop("!is.double(x)");
    if(!is.double(w)) stop("!is.double(w)");
    if(!is.double(z)) stop("!is.double(z)");
    if(!is.double(b)) stop("!is.double(b)");
    if(!is.double(g0)) stop("!is.double(g0)");
    if(!is.double(d0)) stop("!is.double(d0)");
    if(!is.double(G)) stop("!is.double(G)");
    if(!is.double(D)) stop("!is.double(D)");
    if(!is.double(lambda_b)) stop("!is.double(lambda_b)");
    if(!is.double(lambda_g0)) stop("!is.double(lambda_g0)");
    if(!is.double(lambda_d0)) stop("!is.double(lambda_d0)");
    if(!is.double(lambda_G)) stop("!is.double(lambda_G)");
    if(!is.double(lambda_D)) stop("!is.double(lambda_D)");

    ans = .C("gr_sqrLossWithL2",
        # OUTPUT
        output$gr_b, output$gr_g0, output$gr_d0, output$gr_G, output$gr_D,
        # INPUT
        user, item, y, x, w, z,
        b, g0, d0, G, D,
        lambda_b, lambda_g0, lambda_d0, lambda_G, lambda_D,
        nObs, nUsers, nItems, nFactors, nJointFeatures, nUserFeatures, nItemFeatures,
        # OTHER
        as.integer(debug),
        DUP=FALSE    
    );
    return(output);
}


mainEffect_gradient <- function(
    thisEffIndex, otherEffIndex, y_minus_xb_minus_uv,
    thisEff, fittedEff, otherEff, var_y, var_eff, debug
){
    nObs = as.integer(length(thisEffIndex));
    nThisEff  = as.integer(length(fittedEff));
    nOtherEff = as.integer(length(otherEff));
    
    if(length(thisEff) != nThisEff) stop("length(thisEff) != nThisEff");
    if(length(otherEffIndex) != nObs) stop("length(otherEffIndex) != nObs");
    if(length(y_minus_xb_minus_uv) != nObs) stop("length(y_minus_xb_minus_uv) != nObs");
    if(length(var_y) != 1 || length(var_eff) != 1) stop("length(var_y) != 1 || length(var_eff) != 1");
    
    out = rep(0.0, nThisEff);

    if(!is.double(y_minus_xb_minus_uv)) stop("!is.double(y_minus_xb_minus_uv)");
    if(!is.double(fittedEff)) stop("!is.double(fittedEff)");
    if(!is.double(thisEff))   stop("!is.double(thisEff)");
    if(!is.double(otherEff))  stop("!is.double(otherEff)");
    if(!is.double(var_y))     stop("!is.double(var_y)");
    if(!is.double(var_eff))   stop("!is.double(var_eff)");
    
    ans = .C("mainEffect_gradient",
        # OUTPUT
        out,
        # INPUT
        thisEff, as.integer(thisEffIndex), as.integer(otherEffIndex), y_minus_xb_minus_uv,
        fittedEff, otherEff, var_y, var_eff,
        nObs, nThisEff, nOtherEff, as.integer(debug),
        DUP=FALSE
    );
    return(out);
}

factor_gradient <- function(
    thisEffIndex, otherEffIndex, y_minus_xb_alpha_beta,
    thisEff, fittedEff, otherEff, var_y, var_eff, oi, debug
){

    nObs = as.integer(length(thisEffIndex));
    nThisEff  = as.integer(nrow(fittedEff));
    nOtherEff = as.integer(nrow(otherEff));
    nFactors  = as.integer(ncol(fittedEff));
    
    if(ncol(otherEff) != nFactors) stop("ncol(otherEff) != nFactors");
    if(nrow(thisEff) != nThisEff) stop("nrow(thisEff) != nThisEff");
    if(ncol(thisEff) != nFactors) stop("ncol(thisEff) != nFactors");
    if(length(otherEffIndex) != nObs) stop("length(otherEffIndex) != nObs");
    if(length(y_minus_xb_alpha_beta) != nObs) stop("length(y_minus_xb_alpha_beta) != nObs");
    if(length(var_y) != 1 || length(var_eff) != 1) stop("length(var_y) != 1 || length(var_eff) != 1");
    if(length(oi$obsIndex) != nObs) stop("length(oi$obsIndex) != nObs");
    if(length(oi$start) != nThisEff) stop("length(oi$start) != nThisEff");
    if(length(oi$num) != nThisEff) stop("length(oi$num) != nThisEff");
    
    out = matrix(0.0, nrow=nThisEff, ncol=nFactors);

    if(!is.double(thisEff))   stop("!is.double(thisEff)");
    if(!is.double(y_minus_xb_alpha_beta)) stop("!is.double(y_minus_xb_alpha_beta)");
    if(!is.double(fittedEff)) stop("!is.double(fittedEff)");
    if(!is.double(otherEff)) stop("!is.double(otherEff)");
    if(!is.double(var_y)) stop("!is.double(var_y)");
    if(!is.double(var_eff)) stop("!is.double(var_eff)");
    if(!is.integer(oi$obsIndex)) stop("!is.integer(oi$obsIndex)");
    if(!is.integer(oi$start)) stop("!is.integer(oi$start)");
    if(!is.integer(oi$num)) stop("!is.integer(oi$num)");
    
    ans = .C("factor_gradient",
        # OUTPUT
        out,
        # INPUT
        thisEff, as.integer(thisEffIndex), as.integer(otherEffIndex), y_minus_xb_alpha_beta,
        fittedEff, otherEff, var_y, var_eff,
        nObs, nThisEff, nOtherEff, nFactors,
        oi$obsIndex, oi$start, oi$num,
        as.integer(debug),
        DUP=FALSE
    );

    return(out);
}


mainEffect_condMeanVarSample <- function(
    option, # 1:Sample, 2:Mean&Var, 3:Sample&Mean&Var
    thisEffIndex, otherEffIndex, y_minus_xb_minus_uv,
    fittedEff, otherEff, var_y, var_eff, debug, version=1
){
    nObs = as.integer(length(thisEffIndex));
    nThisEff  = as.integer(length(fittedEff));
    nOtherEff = as.integer(length(otherEff));
    nVar_y    = as.integer(length(var_y));
    nVar_eff  = as.integer(length(var_eff));
    
    if(length(otherEffIndex) != nObs) stop("length(otherEffIndex) != nObs");
    if(length(y_minus_xb_minus_uv) != nObs) stop("length(y_minus_xb_minus_uv) != nObs");
    if(!(nVar_y == 1 || (version >= 2 && nVar_y == nObs))) stop("length(var_y) has problem");
    if(!(nVar_eff == 1 || (version >= 2 && nVar_eff == nThisEff))) stop("length(var_eff) has problem");
    
    out = list(sample=as.double(NULL), mean=as.double(NULL), var=as.double(NULL));
    if(option == 1 || option == 3)  out$sample = rep(0.0, nThisEff);
    if(option == 2 || option == 3){ out$mean   = rep(0.0, nThisEff);   out$var = rep(0.0, nThisEff);}

    if(!is.double(y_minus_xb_minus_uv)) stop("!is.double(y_minus_xb_minus_uv)");
    if(!is.double(fittedEff)) stop("!is.double(fittedEff)");
    if(!is.double(otherEff)) stop("!is.double(otherEff)");
    if(!is.double(var_y)) stop("!is.double(var_y)");
    if(!is.double(var_eff)) stop("!is.double(var_eff)");
    
    if(version == 1){
        ans = .C("mainEffect_condMeanVarSample",
            # OUTPUT
            out$sample, out$mean, out$var,
            # INPUT
            as.integer(option), as.integer(thisEffIndex), as.integer(otherEffIndex), y_minus_xb_minus_uv,
            fittedEff, otherEff, var_y, var_eff, 
            nObs, nThisEff, nOtherEff, as.integer(length(var_y)),
            as.integer(debug), DUP=FALSE
        );
    }else if(version == 2){
        ans = .C("mainEffect_condMeanVarSample2",
            # OUTPUT
            out$sample, out$mean, out$var,
            # INPUT
            as.integer(option), as.integer(thisEffIndex), as.integer(otherEffIndex), y_minus_xb_minus_uv,
            fittedEff, otherEff, var_y, var_eff,
            nObs, nThisEff, nOtherEff, nVar_y, nVar_eff,
            as.integer(debug), DUP=FALSE
        );
    }else stop("Unknown version number: ",version);
    return(out);
}

generateObsIndex <- function(
    effIndex, debug
){
    nEff = as.integer(max(effIndex));
    nObs   = as.integer(length(effIndex));
    
    out = list(
        obsIndex = integer(nObs),
        start    = integer(nEff),
        num      = integer(nEff)
    );
    
    if(!is.integer(effIndex)) stop("!is.integer(effIndex)");
    
    ans = .C("generateObsIndex",
        out$obsIndex, out$start, out$num,
        # INPUT
        effIndex, nObs, nEff,
        # OTHER
        as.integer(debug),
        DUP=FALSE
    );
    return(out);
}

factor_condMeanVarSample <- function(
    option, # 1:Sample, 2:Mean&Var, 3:Sample&Mean&Var
    thisEffIndex, otherEffIndex, y_minus_xb_alpha_beta,
    fittedEff, otherEff, var_y, var_eff, oi, debug, version=1
){

    nObs = as.integer(length(thisEffIndex));
    nThisEff  = as.integer(nrow(fittedEff));
    nOtherEff = as.integer(nrow(otherEff));
    nFactors  = as.integer(ncol(fittedEff));
    
    nVar_y = as.integer(length(var_y));
    nVar_eff = as.integer(length(var_eff));
    if(nVar_eff > 0){
        nVar_eff = as.integer(nrow(var_eff));
    }
    
    if(ncol(otherEff) != nFactors) stop("ncol(otherEff) != nFactors");
    if(length(otherEffIndex) != nObs) stop("length(otherEffIndex) != nObs");
    if(length(y_minus_xb_alpha_beta) != nObs) stop("length(y_minus_xb_alpha_beta) != nObs");
    if(length(oi$obsIndex) != nObs) stop("length(oi$obsIndex) != nObs");
    if(length(oi$start) != nThisEff) stop("length(oi$start) != nThisEff");
    if(length(oi$num) != nThisEff) stop("length(oi$num) != nThisEff");

    if(!(nVar_y == 1 || (version >= 2 && nVar_y == nObs))) stop("length(var_y) has problem");
    if(!(nVar_eff == 1 || (version >= 2 && nVar_eff == nThisEff))) stop("length(var_eff) has problem");
    
    out = list(sample=as.double(NULL), mean=as.double(NULL), var=as.double(NULL));
    if(option == 1 || option == 3)  out$sample = matrix(0.0, nrow=nThisEff, ncol=nFactors);
    if(option == 2 || option == 3){ out$mean   = matrix(0.0, nrow=nThisEff, ncol=nFactors);  out$var = array(0.0, dim=c(nThisEff, nFactors, nFactors));}

    if(!is.double(y_minus_xb_alpha_beta)) stop("!is.double(y_minus_xb_alpha_beta)");
    if(!is.double(fittedEff)) stop("!is.double(fittedEff)");
    if(!is.double(otherEff)) stop("!is.double(otherEff)");
    if(!is.double(var_y)) stop("!is.double(var_y)");
    if(!is.double(var_eff)) stop("!is.double(var_eff)");
    if(!is.integer(oi$obsIndex)) stop("!is.integer(oi$obsIndex)");
    if(!is.integer(oi$start)) stop("!is.integer(oi$start)");
    if(!is.integer(oi$num)) stop("!is.integer(oi$num)");
    
    if(version == 1){
        ans = .C("factor_condMeanVarSample",
            # OUTPUT
            out$sample, out$mean, out$var,
            # INPUT
            as.integer(option), as.integer(thisEffIndex), as.integer(otherEffIndex), y_minus_xb_alpha_beta,
            fittedEff, otherEff, var_y, var_eff,
            nObs, nThisEff, nOtherEff, nFactors,
			oi$obsIndex, oi$start, oi$num, as.integer(length(var_y)),
            as.integer(debug),
            DUP=FALSE
        );
    }else if(version == 2){
        ans = .C("factor_condMeanVarSample2",
            # OUTPUT
            out$sample, out$mean, out$var,
            # INPUT
            as.integer(option), as.integer(thisEffIndex), as.integer(otherEffIndex), y_minus_xb_alpha_beta,
            fittedEff, otherEff, var_y, var_eff,
            nObs, nThisEff, nOtherEff, nFactors, nVar_y, nVar_eff,
            oi$obsIndex, oi$start, oi$num,
            as.integer(debug),
            DUP=FALSE
        );
    }else stop("Unknown version number: ", version);
    
    return(out);
}

###
### The E-Step for the Monte-Carlo EM algorithm (version 1)
###
MC_EStep_logistic_adaptive_sampling.C <- function(
    nSamples, nBurnin,
    user, item, y, x, w, z,
    alpha, beta, u, v,
    ars_XI_alpha, ars_XI_beta,
    ars_XI_u, ars_XI_v,
    b, g0, G, d0, D,
    var_alpha, var_beta, var_u, var_v=1, 
    ars_ninit, ars_qcent, ars_xl, ars_xu, ars_alpha,
    debug=0
){
    if(debug >= 1) check.input.logistic(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_alpha, var_beta, var_u, var_v);
    nObs     = as.integer(length(y));
    nUsers   = as.integer(length(alpha));
    nItems   = as.integer(length(beta));
    nFactors = as.integer(ncol(u));

    xb  = x %*% b;
    g0w = w %*% g0;
    d0z = z %*% d0;
    Gw  = w %*% G;
    Dz  = z %*% D;

    output = list(
    	o.mean       = rep(0.0, nObs),
        alpha.mean   = rep(0.0, nUsers),
        alpha.sumvar = 0.0,
        beta.mean    = rep(0.0, nItems),
        beta.sumvar  = 0.0,
        u.mean       = matrix(0.0, nrow=nUsers, ncol=nFactors),
        u.sumvar     = 0.0,
        v.mean       = matrix(0.0, nrow=nItems, ncol=nFactors),
        v.sumvar     = 0.0,
	ars_XI_alpha = as.double(ars_XI_alpha),
	ars_XI_beta = as.double(ars_XI_beta),
	ars_XI_u = as.double(ars_XI_u),
	ars_XI_v = as.double(ars_XI_v)
    );

    if(!is.integer(user)) stop("!is.integer(user)");
    if(!is.integer(item)) stop("!is.integer(item)");
    if(!is.double(y)) stop("!is.double(y)");
    if(!is.double(xb)) stop("!is.double(xb)");
    if(!is.double(g0w)) stop("!is.double(g0w)");
    if(!is.double(d0z)) stop("!is.double(d0z)");
    if(!is.double(Gw)) stop("!is.double(Gw)");
    if(!is.double(Dz)) stop("!is.double(Dz)");
    if(!is.double(alpha)) stop("!is.double(alpha)");
    if(!is.double(beta)) stop("!is.double(beta)");
    if(!is.double(u)) stop("!is.double(u)");
    if(!is.double(v)) stop("!is.double(v)");
    if(!is.double(var_alpha)) stop("!is.double(var_alpha)");
    if(!is.double(var_beta)) stop("!is.double(var_beta)");
    if(!is.double(var_u)) stop("!is.double(var_u)");
    if(!is.double(var_v)) stop("!is.double(var_v)");

    .C("MCEM_EStep_logistic_adaptive_sampling",
        # OUTPUT
	output$o.mean,
        output$alpha.mean, output$alpha.sumvar,
        output$beta.mean,  output$beta.sumvar,
        output$u.mean,     output$u.sumvar,
        output$v.mean,     output$v.sumvar,
	output$ars_XI_alpha, output$ars_XI_beta,
	output$ars_XI_u, output$ars_XI_v,
        # INPUT
        as.integer(nSamples), as.integer(nBurnin),
        user, item, y, xb, g0w, d0z, Gw, Dz,
        alpha, beta, u, v,
        var_alpha, var_beta, var_u, var_v,
        nObs, nUsers, nItems, nFactors,
	as.integer(ars_ninit), as.double(ars_qcent), 
	as.double(ars_xl), as.double(ars_xu), as.double(ars_alpha),
        # OTHER
        as.integer(debug),
        DUP=FALSE
    );
    return(output);
}
MC_EStep_logistic_arsc.C <- function(
    nSamples, nBurnin,
    user, item, y, x, w, z,
    alpha, beta, u, v,
    ars_XI_alpha, ars_XI_beta,
    ars_XI_u, ars_XI_v,
    b, g0, G, d0, D,
    var_alpha, var_beta, var_u, var_v=1, 
    ars_ninit, ars_qcent, ars_xl, ars_xu, ars_alpha, 
    debug=0,
    beta.int = F, center=F,  main.effects = F
){
    if(debug >= 1) check.input.logistic(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_alpha, var_beta, var_u, var_v);
    nObs     = as.integer(length(y));
    nUsers   = as.integer(length(alpha));
    nItems   = as.integer(length(beta));
    nFactors = as.integer(ncol(u));

    xb  = as.double(x %*% b);
    g0w = as.double(w %*% g0);
    if(beta.int==F) d0z = as.double(z %*% d0) else d0z = as.double(cbind(1,z) %*% d0);
    Gw  = as.double(w %*% G);
    Dz  = as.double(z %*% D);

    output = list(
    	o.mean       = rep(0.0, nObs),
        alpha.mean   = rep(0.0, nUsers),
        alpha.sumvar = 0.0,
        beta.mean    = rep(0.0, nItems),
        beta.sumvar  = 0.0,
        u.mean       = matrix(0.0, nrow=nUsers, ncol=nFactors),
        u.sumvar     = 0.0,
        v.mean       = matrix(0.0, nrow=nItems, ncol=nFactors),
        v.sumvar     = 0.0,
	ars_XI_alpha = as.double(ars_XI_alpha),
	ars_XI_beta = as.double(ars_XI_beta),
	ars_XI_u = as.double(ars_XI_u),
	ars_XI_v = as.double(ars_XI_v)
    );

    if(!is.integer(user)) stop("!is.integer(user)");
    if(!is.integer(item)) stop("!is.integer(item)");
    if(!is.double(y)) stop("!is.double(y)");
    if(!is.double(xb)) stop("!is.double(xb)");
    if(!is.double(g0w)) stop("!is.double(g0w)");
    if(!is.double(d0z)) { cat(d0z,"\n"); stop("!is.double(d0z)"); }
    if(!is.double(Gw)) stop("!is.double(Gw)");
    if(!is.double(Dz)) stop("!is.double(Dz)");
    if(!is.double(alpha)) stop("!is.double(alpha)");
    if(!is.double(beta)) stop("!is.double(beta)");
    if(!is.double(u)) stop("!is.double(u)");
    if(!is.double(v)) stop("!is.double(v)");
    if(!is.double(var_alpha)) stop("!is.double(var_alpha)");
    if(!is.double(var_beta)) stop("!is.double(var_beta)");
    if(!is.double(var_u)) stop("!is.double(var_u)");
    if(!is.double(var_v)) stop("!is.double(var_v)");
    if(!is.logical(main.effects)) stop("!is.logical(main.effects)");
    if(!is.logical(center)) stop("!is.logical(center)");
    if(!is.logical(beta.int)) stop("is.logical(beta.int)");
#browser()
    .C("MCEM_EStep_arsc",
        # OUTPUT
	output$o.mean,
        output$alpha.mean, output$alpha.sumvar,
        output$beta.mean,  output$beta.sumvar,
        output$u.mean,     output$u.sumvar,
        output$v.mean,     output$v.sumvar,
	output$ars_XI_alpha, output$ars_XI_beta,
	output$ars_XI_u, output$ars_XI_v,
        # INPUT
        as.integer(nSamples), as.integer(nBurnin),
        user, item, y, xb, g0w, d0z, Gw, Dz,
        alpha, beta, u, v,
        var_alpha, var_beta, var_u, var_v,
        nObs, nUsers, nItems, nFactors,
	as.integer(ars_ninit), as.double(ars_qcent), 
	as.double(ars_xl), as.double(ars_xu), as.double(ars_alpha),
        # OTHER
        as.integer(debug),
        as.integer(main.effects),
        as.integer(beta.int),
        as.integer(center),
        DUP=FALSE
    );
#browser()
    return(output);
}

MC_EStep_logistic_arscid.C <- function(
    nSamples, nBurnin,
    user, item, y, x, w, z,
    alpha, beta, u, v,
    ars_XI_alpha, ars_XI_beta,
    ars_XI_u, ars_XI_v,
    b, g0, G, d0, D,
    var_alpha, var_beta, var_u, var_v=rep(1,nFactors),
    ars_ninit, ars_qcent, ars_xl, ars_xu, ars_alpha,
    debug=0,
    beta.int = F, center=F,  main.effects = F
){
    if(debug >= 1) check.input.logistic(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_alpha, var_beta, var_u, var_v);
    nObs     = as.integer(length(y));
    nUsers   = as.integer(length(alpha));
    nItems   = as.integer(length(beta));
    nFactors = as.integer(ncol(u));

    xb  = as.double(x %*% b);
    g0w = as.double(w %*% g0);
    if(beta.int==F) d0z = as.double(z %*% d0) else d0z = as.double(cbind(1,z) %*% d0);
    Gw  = as.double(w %*% G);
    Dz  = as.double(z %*% D);

    output = list(
        o.mean       = rep(0.0, nObs),
        alpha.mean   = rep(0.0, nUsers),
        alpha.sumvar = 0.0,
        beta.mean    = rep(0.0, nItems),
        beta.sumvar  = 0.0,
        u.mean       = matrix(0.0, nrow=nUsers, ncol=nFactors),
        u.sumvar     = rep(0.0,nFactors),
        v.mean       = matrix(0.0, nrow=nItems, ncol=nFactors),
        v.sumvar     = rep(0.0,nFactors),
        ars_XI_alpha = as.double(ars_XI_alpha),
        ars_XI_beta = as.double(ars_XI_beta),
        ars_XI_u = as.double(ars_XI_u),
        ars_XI_v = as.double(ars_XI_v)
    );
    if(!is.integer(user)) stop("!is.integer(user)");
    if(!is.integer(item)) stop("!is.integer(item)");
    if(!is.double(y)) stop("!is.double(y)");
    if(!is.double(xb)) stop("!is.double(xb)");
    if(!is.double(g0w)) stop("!is.double(g0w)");
    if(!is.double(d0z)) { cat(d0z,"\n"); stop("!is.double(d0z)"); }
    if(!is.double(Gw)) stop("!is.double(Gw)");
    if(!is.double(Dz)) stop("!is.double(Dz)");
    if(!is.double(alpha)) stop("!is.double(alpha)");
    if(!is.double(beta)) stop("!is.double(beta)");
    if(!is.double(u)) stop("!is.double(u)");
    if(!is.double(v)) stop("!is.double(v)");
    if(!is.double(var_alpha)) stop("!is.double(var_alpha)");
    if(!is.double(var_beta)) stop("!is.double(var_beta)");
    if(!is.double(var_u)) stop("!is.double(var_u)");
    if(!is.double(var_v)) stop("!is.double(var_v)");
    if(!is.logical(main.effects)) stop("!is.logical(main.effects)");
    if(!is.logical(center)) stop("!is.logical(center)");
    if(!is.logical(beta.int)) stop("is.logical(beta.int)");
    .C("MCEM_EStep_arscid",
        # OUTPUT
        output$o.mean,
        output$alpha.mean, output$alpha.sumvar,
        output$beta.mean,  output$beta.sumvar,
        output$u.mean,     output$u.sumvar,
        output$v.mean,     output$v.sumvar,
        output$ars_XI_alpha, output$ars_XI_beta,
        output$ars_XI_u, output$ars_XI_v,
        # INPUT
        as.integer(nSamples), as.integer(nBurnin),
        user, item, y, xb, g0w, d0z, Gw, Dz,
        alpha, beta, u, v,
        var_alpha, var_beta, var_u, var_v,
        nObs, nUsers, nItems, nFactors,
        as.integer(ars_ninit), as.double(ars_qcent),
        as.double(ars_xl), as.double(ars_xu), as.double(ars_alpha),
        # OTHER
        as.integer(debug),
        as.integer(main.effects),
        as.integer(beta.int),
        as.integer(center),
        DUP=FALSE
    );
#browser()
    return(output);
}