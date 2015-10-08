### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
###
### Author: Liang Zhang

# compute probability given eta and alpha. eta here could be a vector of values
get.splinep <- function(knotval,eta){
    neg <- eta < 0
	    A <- knotval
		  ans <- rep(NA,length(eta))
		    ans[neg] <- 2*A/(1+exp(-2.0*(1-A)*eta[neg]))
			  ans[!neg] <- 2*A - 1 + 2*(1-A)/(1+exp(-2*A*eta[!neg]))
			    ans
}

splinefn <- function(knotval,etapos,etaneg){
    -(sum(log(get.splinep(knotval,etapos))) + sum(log(1-get.splinep(knotval,etaneg))))
}

estalpha <- function(Y,mu,o){
    pos <- Y>0
	    etapos <- o[pos]+mu; etaneg <- o[!pos]+mu
		  optimize(f=splinefn,lower=.001,upper=.8,etapos=etapos,etaneg=etaneg)[[1]]
}

predict.y.from.factors <- function(user, item, x, alpha, beta, u, v, b, use.C=FALSE){
    if(use.C) return(x %*% b + alpha[user] + beta[item] + sum_margin(u[user,,drop=FALSE] * v[item,,drop=FALSE], 1))
    else      return(x %*% b + alpha[user] + beta[item] + apply(u[user,,drop=FALSE] * v[item,,drop=FALSE], 1, sum));
}

predict.from.factors <- function(user, item, x, alpha, beta, u, v, b, is.logistic, use.C=FALSE){
    pred.y = predict.y.from.factors(user, item, x, alpha, beta, u, v, b, use.C);
	if(is.logistic){
		pred.y = 1/(1+exp(-pred.y));
	}
	return(pred.y);
}

check.input.logistic <- function(
    user, item, y, x, w, z,
    alpha, beta, u, v,
    b, g0, G, d0, D,
    var_alpha, var_beta, var_u, var_v=1,
    version=1,
    check.NA=FALSE
){
    if(!is.vector(b))  stop("b should be a vector");
    if(!is.vector(g0) && !is.matrix(g0)) stop("g0 should be a vector");
    if(!is.vector(d0) && !is.matrix(d0)) stop("d0 should be a vector");
    if(!is.matrix(G))  stop("G should be a matrix");
    if(!is.matrix(D))  stop("D should be a matrix");
    if(!is.vector(y))  stop("y should be a vector");
    if(!is.vector(user))   stop("user should be a vector");
    if(!is.vector(item))   stop("item should be a vector");

    nObs   = length(y);
    nUsers = length(alpha);
    nItems = length(beta);
    nJointFeatures = length(b);
    nUserFeatures  = length(g0);
    nItemFeatures  = length(d0);
    nFactors       = ncol(G);

    check.individual("feature$x_obs", x, c("double", "dgCMatrix"), c(nObs, nJointFeatures), isNullOK=FALSE, stopIfAnyNull=list("param$b"=b), check.NA=check.NA);
    check.individual("feature$x_src", w, c("double", "dgCMatrix"), c(nUsers, nUserFeatures), isNullOK=FALSE, stopIfAnyNull=list("param$g0"=g0), check.NA=check.NA);
    check.individual("feature$x_dst", z, c("double", "dgCMatrix"), c(nItems, nItemFeatures), isNullOK=FALSE, stopIfAnyNull=list("param$d0"=d0), check.NA=check.NA);

    check.individual(
            "factor$alpha", alpha, "double", list(c(nUsers, 1), nUsers), isNullOK=FALSE,
            stopIfAnyNull=list("obs$y"=y,"obs$src.id"=user,"feature$x_src"=w,"param$g0"=g0,"param$var_alpha"=var_alpha),
            check.NA=check.NA
    );
    check.individual(
            "factor$beta", beta, "double", list(c(nItems, 1), nItems), isNullOK=FALSE,
            stopIfAnyNull=list("obs$y"=y,"obs$dst.id"=item,"feature$x_dst"=z,"param$d0"=d0,"param$var_beta"=var_beta),
            check.NA=check.NA
    );
    check.individual(
            "factor$u", u, "double", c(nUsers, nFactors), isNullOK=(out$nFactors==0),
            stopIfAnyNull=list("obs$y"=y,"obs$src.id"=user,"feature$x_src"=w,"param$G"=G,"param$var_u"=var_u),
            check.NA=check.NA
    );
    check.individual(
            "factor$v", v, "double", c(nItems, nFactors), isNullOK=(out$nFactors==0),
            stopIfAnyNull=list("obs$y"=y,"obs$dst.id"=item,"feature$x_dst"=z,"param$D"=D,"param$var_v"=var_v),
            check.NA=check.NA
    );
    
    if(version == 1){
        if(!length(var_alpha) == 1) stop("var_alpha should have length 1");
        if(!length(var_beta) == 1)  stop("var_beta should have length 1");
	if(!length(var_u) == 1 && !length(var_u)==nFactors)     stop("var_u should have length 1 or nFactors");
        if(!length(var_v) == 1 && !length(var_v)==nFactors)     stop("var_v should have length 1 or nFactors");

        if(var_alpha < 0) stop("var_alpha < 0");
        if(var_beta < 0)  stop("var_beta < 0");
        if(any(var_u < 0))     stop("var_u < 0");
        if(any(var_v < 0))     stop("var_v < 0");
    }else if(version == 2){
        if(!(length(var_alpha) == 1 || length(var_alpha) == length(alpha))) stop("var_alpha should have length 1 or length(alpha)");
        if(!(length(var_beta) == 1  || length(var_beta)  == length(beta)))  stop("var_beta should have length 1 or length(beta)");
        if(!(length(var_u) == 1 || all(dim(var_u) == c(nUsers, nFactors, nFactors)))) stop("var_u should have length 1 or nUsers x nFactors x nFactors");
        if(!(length(var_v) == 1 || all(dim(var_v) == c(nItems, nFactors, nFactors)))) stop("var_v should have length 1 or nItems x nFactors x nFactors");

        if(any(var_alpha < 0)) stop("var_alpha < 0");
        if(any(var_beta < 0))  stop("var_beta < 0");
        if(length(var_u) == 1){
            if(var_u < 0) stop("var_u < 0");
        }else{
            for(f in 1:nFactors){ if(any(var_u[,f,f] < 0)) stop("var_u < 0");}
        }
        if(length(var_v) == 1){
            if(var_v < 0) stop("var_v < 0");
        }else{
            for(f in 1:nFactors){ if(any(var_v[,f,f] < 0)) stop("var_v < 0");}
        }
    }else stop("Unkown version number: version = ",version);

    if(ncol(D) != nFactors) stop("ncol(D) != nFactors");
    if(nrow(G) != nUserFeatures) stop("nrow(G) != nUserFeatures");
    if(nrow(D) != nItemFeatures) stop("nrow(D) != nItemFeatures");
    if(nObs < nUsers || nObs < nItems) stop("nObs < nUsers || nObs < nItems");
    if(length(user) != nObs) stop("length(user) != nObs");
    if(length(item) != nObs) stop("length(item) != nObs");
}

logLikelihood.logistic.old <- function(
    user, item, y, x, w, z,
    alpha, beta, u, v,
    b, g0, G, d0, D,
    var_alpha, var_beta, var_u, var_v=1, debug=0, use.C=FALSE
){
    if(debug >= 1) check.input.logistic(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_alpha, var_beta, var_u, var_v);
    nObs     = length(y);
    nUsers   = length(alpha);
    nItems   = length(beta);
    nFactors = ncol(u);

    ans = 0;
    o = predict.y.from.factors(user, item, x, alpha, beta, u, v, b, use.C);
    ans = ans + sum(y*o) - sum(log(1+exp(o)));
    if (length(var_u)==1) {
       err = u - w %*% G;
       ans = ans - (sum(err^2) / var_u + nUsers * nFactors * log(var_u))/2;
       err = v - z %*% D;
       ans = ans - (sum(err^2) / var_v + nItems * nFactors * log(var_v))/2;
    } else {
       err = u - w %*% G;
       for (k in 1:nFactors) {
       	   ans = ans - (sum(err[,k]^2)/var_u[k] + nUsers * log(var_u[k]))/2;
       }
       err = v - z %*% D;
       for (k in 1:nFactors) {
       	   ans = ans - (sum(err[,k]^2)/var_v[k] + nItems * log(var_v[k]))/2;
       }
    }
    err = alpha - w %*% g0;
    ans = ans - (sum(err^2) / var_alpha + nUsers * log(var_alpha))/2;
    err = beta - z %*% d0;
    ans = ans - (sum(err^2) / var_beta  + nItems * log(var_beta))/2;
    return(ans);
}

logLikelihood.logistic <- function(
    user, item, y, x, w, z,
    alpha, beta, u, v,
    b, g0, G, d0, D,
    var_alpha, var_beta, var_u, var_v=1, ars_alpha=0.5,
    beta.int = F, debug=0, use.C=FALSE
){
    if(debug >= 1) check.input.logistic(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_alpha, var_beta, var_u, var_v);
    nObs     = length(y);
    nUsers   = length(alpha);
    nItems   = length(beta);
    nFactors = ncol(u);

    ans = 0;
    o = predict.y.from.factors(user, item, x, alpha, beta, u, v, b, use.C);
    p = 2*ars_alpha/(1+exp(-2*(1-ars_alpha)*o));
    p[o>=0] = 2*ars_alpha - 1 + 2*(1-ars_alpha)/(1+exp(-2*ars_alpha*o[o>=0]));
    ans = ans + sum(y*log(p)+(1-y)*log(1-p));
    #ans = ans + sum(y*o) - sum(log(1+exp(o)));
    if (length(var_u)==1) {
       err = u - w %*% G;
       ans = ans - (sum(err^2) / var_u + nUsers * nFactors * log(var_u))/2;
       err = v - z %*% D;
       ans = ans - (sum(err^2) / var_v + nItems * nFactors * log(var_v))/2;
    } else {
       err = u - w %*% G;
       for (k in 1:nFactors) {
           ans = ans - (sum(err[,k]^2)/var_u[k] + nUsers * log(var_u[k]))/2;
       }
       err = v - z %*% D;
       for (k in 1:nFactors) {
           ans = ans - (sum(err[,k]^2)/var_v[k] + nItems * log(var_v[k]))/2;
       }
    }
    err = alpha - w %*% g0;
    ans = ans - (sum(err^2) / var_alpha + nUsers * log(var_alpha))/2;
    if (beta.int==F) err = beta - z %*% d0 else err = beta - cbind(1,z) %*% d0;
    ans = ans - (sum(err^2) / var_beta  + nItems * log(var_beta))/2;
    return(ans);
}

