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

check.input <- function(
    user, item, y, x, w, z,
    alpha, beta, u, v,
    b, g0, G, d0, D,
    var_y, var_alpha, var_beta, var_u, var_v=1,
    version=1
){
    if(!is.vector(b))  stop("b should be a vector");
    if(!is.vector(g0)) stop("g0 should be a vector");
    if(!is.vector(d0)) stop("d0 should be a vector");
    if(!is.matrix(G))  stop("G should be a matrix");
    if(!is.matrix(D))  stop("D should be a matrix");
    if(!is.vector(y))  stop("y should be a vector");
    if(!is.matrix(x))  stop("x should be a matrix");
#    if(!is.matrix(w))  stop("w should be a matrix");
#    if(!is.matrix(z))  stop("z should be a matrix");
    if(!is.matrix(u))  stop("u should be a matrix");
    if(!is.matrix(v))  stop("v should be a matrix");
    if(!is.vector(alpha))  stop("alpha should be a vector");
    if(!is.vector(beta))   stop("beta should be a vector");
    if(!is.vector(user))   stop("user should be a vector");
    if(!is.vector(item))   stop("item should be a vector");
    
    nObs   = length(y);
    nUsers = length(alpha);
    nItems = length(beta);
    nJointFeatures = length(b);
    nUserFeatures  = length(g0);
    nItemFeatures  = length(d0);
    nFactors       = ncol(G);

    if(version == 1){
        if(!(length(var_y) == 1 || length(var_y) == nObs)) stop("var_y should have length 1 or #obs");
        if(!length(var_alpha) == 1) stop("var_alpha should have length 1");
        if(!length(var_beta) == 1)  stop("var_beta should have length 1");
	if(!length(var_u) == 1 && !length(var_u)==nFactors)     stop("var_u should have length 1 or nFactors");
        if(!length(var_v) == 1 && !length(var_v)==nFactors)     stop("var_v should have length 1 or nFactors");
        
        if(any(var_y < 0)) stop("var_y < 0");
        if(var_alpha < 0) stop("var_alpha < 0");
        if(var_beta < 0)  stop("var_beta < 0");
        if(any(var_u < 0))     stop("var_u < 0");
        if(any(var_v < 0))     stop("var_v < 0");
        
    }else if(version == 2){
        if(!(length(var_y) == 1 || length(var_y) == length(y))) stop("var_y should have length 1 or length(y)");
        if(!(length(var_alpha) == 1 || length(var_alpha) == length(alpha))) stop("var_alpha should have length 1 or length(alpha)");
        if(!(length(var_beta) == 1  || length(var_beta)  == length(beta)))  stop("var_beta should have length 1 or length(beta)");
        if(!(length(var_u) == 1 || all(dim(var_u) == c(nUsers, nFactors, nFactors)))) stop("var_u should have length 1 or nUsers x nFactors x nFactors");
        if(!(length(var_v) == 1 || all(dim(var_v) == c(nItems, nFactors, nFactors)))) stop("var_v should have length 1 or nItems x nFactors x nFactors");
        
        if(any(var_y < 0))     stop("var_y < 0");
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
    if(nrow(u) != nUsers)   stop("nrow(u) != nUsers");
    if(ncol(u) != nFactors) stop("ncol(u) != nFactors");
    if(nrow(v) != nItems)   stop("nrow(u) != nItems");
    if(ncol(v) != nFactors) stop("ncol(u) != nFactors");
    if(length(user) != nObs) stop("length(user) != nObs");
    if(length(item) != nObs) stop("length(item) != nObs");
    if(nrow(x) != nObs) stop("nrow(x) != nObs");
    if(ncol(x) != nJointFeatures) stop("ncol(x) != nJointFeatures");
    if(nrow(w) != nUsers)        stop("nrow(w) != nUsers");
    if(ncol(w) != nUserFeatures) stop("ncol(w) != nUserFeatures");
    if(nrow(z) != nItems)        stop("nrow(z) != nItems");
    if(ncol(z) != nItemFeatures) stop("ncol(z) != nItemFeatures");
}

check.input2 <- function(
    user, item, y, x, w, z,
    b, g0, G, d0, D,
    lambda_b, lambda_g0, lambda_G, lambda_d0, lambda_D
){
    if(!is.vector(b))  stop("b should be a vector");
    if(!is.vector(g0)) stop("g0 should be a vector");
    if(!is.vector(d0)) stop("d0 should be a vector");
    if(!is.matrix(G))  stop("G should be a matrix");
    if(!is.matrix(D))  stop("D should be a matrix");
    if(!is.vector(y))  stop("y should be a vector");
    if(!is.matrix(x))  stop("x should be a matrix");
#if(!is.matrix(w))  stop("w should be a matrix");
#    if(!is.matrix(z))  stop("z should be a matrix");
    if(!is.vector(lambda_b))  stop("lambda_b should be a vector");
    if(!is.vector(lambda_g0)) stop("lambda_g0 should be a vector");
    if(!is.vector(lambda_d0)) stop("lambda_d0 should be a vector");
    if(!is.matrix(lambda_G))  stop("lambda_G should be a matrix");
    if(!is.matrix(lambda_D))  stop("lambda_D should be a matrix");
    if(!is.vector(user))   stop("user should be a vector");
    if(!is.vector(item))   stop("item should be a vector");

    nObs   = length(y);
    nUsers = nrow(w);
    nItems = nrow(z);
    nJointFeatures = length(b);
    nUserFeatures  = length(g0);
    nItemFeatures  = length(d0);
    nFactors       = ncol(G);
    
    if(ncol(D) != nFactors) stop("ncol(D) != nFactors");
    if(nrow(G) != nUserFeatures) stop("nrow(G) != nUserFeatures");
    if(nrow(D) != nItemFeatures) stop("nrow(D) != nItemFeatures");
    if(nObs < nUsers || nObs < nItems) stop("nObs < nUsers || nObs < nItems");
    if(length(user) != nObs) stop("length(user) != nObs");
    if(length(item) != nObs) stop("length(item) != nObs");
    if(nrow(x) != nObs) stop("nrow(x) != nObs");
    if(ncol(x) != nJointFeatures) stop("ncol(x) != nJointFeatures");
    if(ncol(w) != nUserFeatures) stop("ncol(w) != nUserFeatures");
    if(ncol(z) != nItemFeatures) stop("ncol(z) != nItemFeatures");
    
    if(length(b)  != length(lambda_b))  stop("length(b)  != length(lambda_b)");
    if(length(g0) != length(lambda_g0)) stop("length(g0) != length(lambda_g0)");
    if(length(d0) != length(lambda_d0)) stop("length(d0) != length(lambda_d0)");
    if(any(dim(G) != dim(lambda_G)))    stop("any(dim(G) != dim(lambda_G))");
    if(any(dim(D) != dim(lambda_D)))    stop("any(dim(D) != dim(lambda_D))");
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
    version=1
){
    if(!is.vector(b))  stop("b should be a vector");
    if(!is.vector(g0)) stop("g0 should be a vector");
    if(!is.vector(d0)) stop("d0 should be a vector");
    if(!is.matrix(G))  stop("G should be a matrix");
    if(!is.matrix(D))  stop("D should be a matrix");
    if(!is.vector(y))  stop("y should be a vector");
    if(!is.matrix(x))  stop("x should be a matrix");
#if(!is.matrix(w))  stop("w should be a matrix");
#    if(!is.matrix(z))  stop("z should be a matrix");
    if(!is.matrix(u))  stop("u should be a matrix");
    if(!is.matrix(v))  stop("v should be a matrix");
    if(!is.vector(alpha))  stop("alpha should be a vector");
    if(!is.vector(beta))   stop("beta should be a vector");
    if(!is.vector(user))   stop("user should be a vector");
    if(!is.vector(item))   stop("item should be a vector");

    nObs   = length(y);
    nUsers = length(alpha);
    nItems = length(beta);
    nJointFeatures = length(b);
    nUserFeatures  = length(g0);
    nItemFeatures  = length(d0);
    nFactors       = ncol(G);

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
    if(nrow(u) != nUsers)   stop("nrow(u) != nUsers");
    if(ncol(u) != nFactors) stop(paste("ncol(u)=",ncol(u)," != ",nFactors,"=nFactors",sep=""));
    if(nrow(v) != nItems)   stop("nrow(v) != nItems");
    if(ncol(v) != nFactors) stop("ncol(u) != nFactors");
    if(length(user) != nObs) stop("length(user) != nObs");
    if(length(item) != nObs) stop("length(item) != nObs");
    if(nrow(x) != nObs) stop("nrow(x) != nObs");
    if(ncol(x) != nJointFeatures) stop("ncol(x) != nJointFeatures");
    if(nrow(w) != nUsers)        stop("nrow(w) != nUsers");
    if(ncol(w) != nUserFeatures) stop("ncol(w) != nUserFeatures");
    if(nrow(z) != nItems)        stop("nrow(z) != nItems");
    if(ncol(z) != nItemFeatures) stop("ncol(z) != nItemFeatures");
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


# log likelihood with the constant term removed
logLikelihood <- function(
    user, item, y, x, w, z,
    alpha, beta, u, v,
    b, g0, G, d0, D,
    var_y, var_alpha, var_beta, var_u, var_v=1, is.logistic=FALSE, debug=0, use.C=FALSE
){
    if(debug >= 1) check.input(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_y, var_alpha, var_beta, var_u, var_v);
    nObs     = length(y);
    nUsers   = length(alpha);
    nItems   = length(beta);
    nFactors = ncol(u);
    
    pred.y = predict.y.from.factors(user, item, x, alpha, beta, u, v, b, use.C);
	ans = obsLoglik.from.gaussian(pred.y, y, var_y, is.logistic);
	if(var_u>0){
		ans = ans + loglik.gaussian(w %*% G, u, var_u);
	}
	if(var_v >0){
		ans = ans + loglik.gaussian(z %*% D, v, var_v);
	}
	if(var_alpha > 0){
		ans = ans + loglik.gaussian(w %*% g0, alpha, var_alpha);
	}
	if(var_beta>0){
		ans = ans + loglik.gaussian(z %*% d0, beta,  var_beta);
	}
	
    return(ans);
}

loglik.gaussian <- function(pred.x, x, var_x){
	if(length(var_x) != 1) stop("length(var_x) != 1");
	loglik = -(1/2) * ( sum((x - pred.x)^2 / var_x) + length(x) * log(var_x) );
	return(loglik);
}

###
### For non-Gaussian response
### 
check.obs <- function(y, is.logistic){
	if(is.logistic){
		labels = unique(y);
		if(length(labels) != 2) stop("The response is not binary: ",paste(labels[1:min(10,length(labels))],collapse=", "));
		labels = sort(labels);
		if(any(labels != c(0,1)) && any(labels != c(-1,1))) stop("Binary response must be {0,1} or {-1,1}");
	}
}

# Generate Gaussian response for Laplace approximation
generate.response.laplace <- function(y, response, var_y, eta, is.logistic, verbose=0){
  if(is.logistic){
    # laplace approximation
    if(length(eta) != length(y)) stop("length(eta) != nObs");
    if(length(response) != length(y)) stop("length(response) != nObs");
    p = 1.0/(1+exp(-eta));
    one_minus_p = 1.0/(1+exp(eta));
    y = eta + (response-p)/(p*one_minus_p);
    var_y = 1/p/one_minus_p;
  }
  return(list(y=y,var_y=var_y));
}


# Generate Gaussain response
generate.response <- function(y, response, var_y, xi, is.logistic, verbose=0){
	if(is.logistic){
		# variational approximation
		if(verbose >= 1) cat("generate gaussian response for logistic\n");
		if(length(xi) != length(y)) stop("length(xi) != nObs");
		if(length(response) != length(y)) stop("length(response) != nObs");
		if(all(response %in% c(0,1))) response[response == 0] = -1;
		if(!all(response %in% c(-1,1))) stop("Binary response must be {-1,1} at this point");
		var_y = 1/(2 * logistic.lambda(xi));
		y = response * var_y / 2;
	}
	return(list(y=y, var_y=var_y));
}
# Update the parameters
update.xi <- function(mc_e, xi, is.logistic){
	if(is.logistic){
		if(length(mc_e$pred.y.square) != length(xi)) { cat("length(mc_e$pred.y.square) =", length(mc_e$pred.y.square), " length(xi)= ",length(xi),"\n",sep="");  stop("length(mc_e$pred.y.square) != length(xi)"); }
		xi = sqrt(mc_e$pred.y.square);
	}
	return(xi);
}
# input pred.y is based on the Gaussian model
obsLoglik.from.gaussian <- function(pred.y, y, var_y, is.logistic){
	if(length(y) != length(pred.y)) stop("length(y) != length(pred.y)");
	if(is.logistic){
		if(all(y %in% c(0,1))) y[y == 0] = -1;
		if(!all(y %in% c(-1,1))) stop("Binary response must be {-1,1} at this point");
		loglik = sum( -log1p(exp(-y * pred.y) ) );
		attr(loglik, "loss") = -loglik / length(y);
	}else{
		loglik = loglik.gaussian(pred.x=pred.y, x=y, var_x=var_y);
		attr(loglik, "loss") = sqrt(mean( (y - pred.y)^2 ));
	}
	return(loglik);
}
# input pred.y is based on the Gaussian model
# output$pred.y is the input for the Gaussian model
#               is the predicted probability for the Logistic model
predict.response.from.gaussian <- function(pred.y, y, var_y, is.logistic){
	if(!is.null(y)) loglik = obsLoglik.from.gaussian(pred.y=pred.y, y=y, var_y=var_y, is.logistic=is.logistic)
	else            loglik = NULL;
	if(is.logistic){
		pred.y = 1/(1+exp(-pred.y));
		if(!is.null(y)){
			if(min(y) == -1) temp.y = 2*pred.y - 1
			else             temp.y = pred.y;
		}
	}else{
		temp.y = pred.y;
	}
	if(!is.null(y)){
		rmse = sqrt(mean( (y - temp.y)^2 ));
		mae  = mean( abs(y - temp.y) );
	}else{
		rmse = NULL;
		mae  = NULL;
	}
	return(list(pred.y=pred.y, true.y=y, rmse=rmse, mae=mae, loglik=loglik, test.loss=attr(loglik, "loss")));
}
logistic.lambda <- function(xi){
	return(tanh(xi/2) / (4*xi));
}

# Conditional mean and variance for alpha
# See output$mean and output$var
alpha.condMeanVar.R <- function(user, item, y, xb, g0w, beta, u, v, var_y, var_alpha){
    output = list();
    nUsers = nrow(u);

    if(length(var_y) == 1) inv_var = rep(1/var_y, length(y))
    else                   inv_var = 1/var_y;
    
    # temp[i,2] is the result for user i has
    temp = aggregate(inv_var, list(user), sum);
    if(any(temp[,1] != 1:nUsers)) stop("any(temp[,1] != 1:nUsers)");
    
    output$var = 1/(temp[,2] + 1/var_alpha);
	
	#cat("sum 1/var_y:\n"); print(temp);
	#cat("var_alpha: ",var_alpha,"\n");
	
    o = y - beta[item] - xb - apply(u[user,,drop=FALSE] * v[item,,drop=FALSE], 1, sum);
    temp = aggregate(o/var_y, list(user), sum);
    if(any(temp[,1] != 1:nUsers)) stop("any(temp[,1] != 1:nUsers)");
    
    output$mean = drop(output$var * (temp[,2] + g0w/var_alpha));
	
    return(output);
}

# Conditional mean and variance for beta
# See output$mean and output$var
beta.condMeanVar.R <- function(user, item, y, xb, d0z, alpha, u, v, var_y, var_beta){
    output = list();
    nItems = nrow(v);

    if(length(var_y) == 1) inv_var = rep(1/var_y, length(y))
    else                   inv_var = 1/var_y;

    # temp[j,2] is the result for item j
    temp = aggregate(inv_var, list(item), sum);
    if(any(temp[,1] != 1:nItems)) stop("any(temp[,1] != 1:nItems)");
    
    output$var = 1/(temp[,2] + 1/var_beta);

    o = y - alpha[user] - xb - apply(u[user,,drop=FALSE] * v[item,,drop=FALSE], 1, sum);
    temp = aggregate(o/var_y, list(item), sum);
    if(any(temp[,1] != 1:nItems)) stop("any(temp[,1] != 1:nItems)");
    
    output$mean = drop(output$var * (temp[,2] + d0z/var_beta));

    return(output);
}

# Conditional mean and variance for u
# See output$mean and output$var 
#   (output$var is a 3D array, var[i,,] is the variance-covariance matrix for user i)
u.condMeanVar.R <- function(user, item, y, xb, Gw, alpha, beta, v, var_y, var_u){
    output = list();
    nUsers   = length(alpha);
    nFactors = ncol(v);

    obsIndex.forUser = tapply(1:length(user), list(user), c, simplify=F);
    if(any(names(obsIndex.forUser) != 1:nUsers)) stop("any(names(obsIndex.forUser) != 1:nUsers)");
    
    output$mean = matrix(NA, nrow=nUsers, ncol=nFactors);
    output$var  = array(NA, dim=c(nUsers, nFactors, nFactors));
    
    for(i in 1:nUsers){
        selectedObs   = obsIndex.forUser[[i]];
        selectedItems = item[selectedObs];
        v.selected = v[selectedItems,,drop=FALSE];
        o.selected = y[selectedObs] - (alpha[i] + beta[selectedItems] + xb[selectedObs]);
        
        if(length(var_y) == 1){
            sum.vv = (t(v.selected) %*% v.selected) / var_y;
            sum.ov = (t(v.selected) %*% o.selected) / var_y;
        }else{
            var_y.selected = var_y[selectedObs];
            sum.vv = t(v.selected) %*% diag(1/var_y.selected, nrow=length(var_y.selected)) %*% v.selected;
            sum.ov = t(v.selected) %*% diag(1/var_y.selected, nrow=length(var_y.selected)) %*% o.selected;
        }
        if(length(var_u) == 1) inv_var_ui = diag(1/var_u, nrow=nFactors)
        else                   inv_var_ui = solve(var_u[i,,]);
        
        output$var[i,,] = solve(sum.vv + inv_var_ui);
                
        # print(drop(sum.ov));
        # print(drop(Gw[i,]));
        
        output$mean[i,] = output$var[i,,] %*% (sum.ov + inv_var_ui %*% Gw[i,]);
    }
    
    return(output);
}

# Conditional mean and variance for v
# See output$mean and output$var 
#   (output$var is a 3D array, var[j,,] is the variance-covariance matrix for item j)
v.condMeanVar.R <- function(user, item, y, xb, Dz, alpha, beta, u, var_y, var_v=1){
    output = list();
    nItems   = length(beta);
    nFactors = ncol(u);

    obsIndex.forItem = tapply(1:length(item), list(item), c, simplify=F);
    if(any(names(obsIndex.forItem) != 1:nItems)) stop("any(names(obsIndex.forItem) != 1:nItems)");
    
    output$mean = matrix(NA, nrow=nItems, ncol=nFactors);
    output$var  = array(NA, dim=c(nItems, nFactors, nFactors));
    
    for(j in 1:nItems){
        selectedObs   = obsIndex.forItem[[j]];
        selectedUsers = user[selectedObs];
        u.selected = u[selectedUsers,,drop=FALSE];
        o.selected = y[selectedObs] - (alpha[selectedUsers] + beta[j] + xb[selectedObs]);

        if(length(var_y) == 1){
            sum.uu = (t(u.selected) %*% u.selected) / var_y;
            sum.ou = (t(u.selected) %*% o.selected) / var_y;
        }else{
            var_y.selected = var_y[selectedObs];
            sum.uu = t(u.selected) %*% diag(1/var_y.selected, nrow=length(var_y.selected)) %*% u.selected;
            sum.ou = t(u.selected) %*% diag(1/var_y.selected, nrow=length(var_y.selected)) %*% o.selected;
        }
        if(length(var_v) == 1) inv_var_vj = diag(1/var_v, nrow=nFactors)
        else                   inv_var_vj = solve(var_v[j,,]);
        
        output$var[j,,] = solve(sum.uu + inv_var_vj);

        # print(drop(sum.ou));
        # print(drop(Dz[j,]));

        output$mean[j,] = output$var[j,,] %*% (sum.ou + inv_var_vj %*% Dz[j,]);
    }
    
    return(output);
}

# Get multivariate normal sample
#   Need to do: library(MASS)
#   mean[k,] is the mean vector for the kth sample point
#   var[k,,] is the variance-covariance matrix for the kth sample point
# output[k,] is the kth sample point
getMVNSample <- function(mean, var, FUN=mvrnorm){
    
    nPoints = nrow(mean);
    nDim    = ncol(mean);
    temp = dim(var);
    if(temp[1] != nPoints || temp[2] != nDim || temp[3] != nDim) stop("size mismatch");
    
    if(nDim == 1) return(matrix(rnorm(nPoints, mean, sqrt(var)), nrow=nPoints, ncol=1));
    
    output = matrix(NA, nrow=nPoints, ncol=nDim);
    for(k in 1:nPoints){
        output[k,] = FUN(1, mu=mean[k,], Sigma=var[k,,]);
    }
    return(output);
}


# Monte-Carlo Expectation to find mean and variance of o, alpha, beta, u and v
#   o = alpha[user] + beta[item] + apply(u[user,,drop=FALSE] * v[item,,drop=FALSE], 1, sum)
#   See output$o.mean, output$o.sumvar, ...
MC_EStep.R <- function(
    nSamples,
    user, item, y, x, w, z,
    alpha, beta, u, v,
    b, g0, G, d0, D,
    var_y, var_alpha, var_beta, var_u, var_v=1, debug=0, 
    func.rmvnorm=mvrnorm, print.path=FALSE,
    outputPerUserVar=FALSE, outputPerItemVar=FALSE,
    isOldUser=NULL, isOldItem=NULL
){
    if(debug >= 1) check.input(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_y, var_alpha, var_beta, var_u, var_v, version=2);
    nObs     = length(y);
    nUsers   = length(alpha);
    nItems   = length(beta);
    nFactors = ncol(u);
    
    # sos means sum of squares; sop means sums of products of pairs of factors
    o.sum     = rep(0, nObs);       o.sos     = rep(0, nObs);
    alpha.sum = rep(0, nUsers);     alpha.sos = rep(0, nUsers);
    beta.sum  = rep(0, nItems);     beta.sos  = rep(0, nItems);

    u.sum = matrix(0, nrow=nUsers, ncol=nFactors);
    v.sum = matrix(0, nrow=nItems, ncol=nFactors);
    
    u.sos = matrix(0, nrow=nUsers, ncol=nFactors);
    v.sos = matrix(0, nrow=nItems, ncol=nFactors);

    if(outputPerUserVar) u.sop = array(0, dim=c(nUsers, nFactors, nFactors))    
    if(outputPerItemVar) v.sop = array(0, dim=c(nItems, nFactors, nFactors))
    
    if(print.path)
        u11 <- v11 <- u21 <- v21 <- rep(NA,length(nSamples));

    xb  = x %*% b;

    g0w = w %*% g0;
    Gw  = w %*% G;
    if(!is.null(isOldUser)){
        g0w[isOldUser] = alpha[isOldUser];
        Gw[isOldUser,] = u[isOldUser,]
    }

    d0z = z %*% d0;
    Dz  = z %*% D;
    if(!is.null(isOldItem)){
        d0z[isOldItem] = beta[isOldItem];
        Dz[isOldItem,] = v[isOldItem,]
    }
    
    for(s in 0:nSamples){
        
        ans = alpha.condMeanVar.R(user, item, y, xb, g0w, beta, u, v, var_y, var_alpha);
        if(length(ans$mean) != nUsers || length(ans$var) != nUsers) stop("alpha.condMeanVar.R output error");
        alpha = rnorm(nUsers, mean=ans$mean, sd=sqrt(ans$var));
        
        ans = beta.condMeanVar.R(user, item, y, xb, d0z, alpha, u, v, var_y, var_beta);
        if(length(ans$mean) != nItems || length(ans$var) != nItems) stop("beta.condMeanVar.R output error");
        beta = rnorm(nItems, mean=ans$mean, sd=sqrt(ans$var));
        
        ans = u.condMeanVar.R(user, item, y, xb, Gw, alpha, beta, v, var_y, var_u);
        if(nrow(ans$mean) != nUsers || nrow(ans$var) != nUsers) stop("u.condMeanVar.R output error");
        u = getMVNSample(ans$mean, ans$var, func.rmvnorm);
        
        ans = v.condMeanVar.R(user, item, y, xb, Dz, alpha, beta, u, var_y, var_v);
        if(nrow(ans$mean) != nItems || nrow(ans$var) != nItems) stop("v.condMeanVar.R output error");
        v = getMVNSample(ans$mean, ans$var, func.rmvnorm);

        if(s > 0){
            o = alpha[user] + beta[item] + apply(u[user,,drop=FALSE] * v[item,,drop=FALSE], 1, sum);
            o.sum = o.sum + o;
            o.sos = o.sos + o^2;
            alpha.sum = alpha.sum + alpha;
            alpha.sos = alpha.sos + alpha^2;
            beta.sum  = beta.sum + beta;
            beta.sos  = beta.sos + beta^2;

            u.sum  = u.sum + u;
            v.sum  = v.sum + v;
            u.sos = u.sos + u^2;
            v.sos = v.sos + v^2;

            if(outputPerUserVar){
                for(i in 1:nUsers){ u.sop[i,,] = u.sop[i,,] + u[i,] %*% t(u[i,]);}
            }
            
            if(outputPerItemVar){
                for(j in 1:nItems){ v.sop[j,,] = v.sop[j,,] + v[j,] %*% t(v[j,]);}
            }
            
            if(print.path){u11[s] <- u[1,1];u21[s] <- u[2,1];v11[s] <- v[1,1];v21[s] <- v[2,1]}
        }
        if(debug >= 1) check.input(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_y, var_alpha, var_beta, var_u, var_v, version=2);
    }
    
    output = list();

    output$o.mean     = o.sum / nSamples;
    output$alpha.mean = alpha.sum / nSamples;
    output$beta.mean  = beta.sum / nSamples;
    output$u.mean     = u.sum / nSamples;
    output$v.mean     = v.sum / nSamples;

    alpha.var = alpha.sos / (nSamples-1) - (nSamples/(nSamples-1)) * output$alpha.mean^2;
    beta.var  = beta.sos  / (nSamples-1) - (nSamples/(nSamples-1)) * output$beta.mean^2;
    
    if(outputPerUserVar){
        u.var = array(NA, dim=c(nUsers, nFactors, nFactors));
        for(i in 1:nUsers){ u.var[i,,] = u.sop[i,,] / (nSamples-1) - (nSamples/(nSamples-1)) * (output$u.mean[i,] %*% t(output$u.mean[i,])); }
        output$alpha.var = alpha.var;
        output$u.var     = u.var;
    }
    if(outputPerItemVar){
        v.var = array(NA, dim=c(nItems, nFactors, nFactors));
        for(i in 1:nItems){ v.var[i,,] = v.sop[i,,] / (nSamples-1) - (nSamples/(nSamples-1)) * (output$v.mean[i,] %*% t(output$v.mean[i,])); }
        output$beta.var = beta.var;
        output$v.var    = v.var;
    }
    
    output$o.sumvar     = sum(o.sos / (nSamples-1) - (nSamples/(nSamples-1)) * output$o.mean^2);
    output$alpha.sumvar = sum(alpha.var);
    output$beta.sumvar  = sum(beta.var);
    output$u.sumvar     = sum(u.sos / (nSamples-1) - (nSamples/(nSamples-1)) * output$u.mean^2);
    output$v.sumvar     = sum(v.sos / (nSamples-1) - (nSamples/(nSamples-1)) * output$v.mean^2);
    
    if(print.path){output$u11 <- u11;output$u21 <- u21;output$v11 <- v11;output$v21 <- v21}
    
    return(output);
}

###
### The following are functions for the ICM method
###

# Compute the derivatives w.r.t. alpha
# See output$d1 and output$d2
alpha.derivatives.R <- function(user, item, y, xb, g0w, alpha, beta, u, v, var_y, var_alpha){
    output = list();
    nUsers = length(alpha);

    # temp[i,2] is the number of observations that user i has
    temp = aggregate(rep(1,length(user)), list(user), sum);
    if(any(temp[,1] != 1:nUsers)) stop("any(temp[,1] != 1:nUsers)");
    
    output$d2 = 2*(temp[,2]/var_y + 1/var_alpha);

    o = y - beta[item] - xb - apply(u[user,,drop=FALSE] * v[item,,drop=FALSE], 1, sum);
    temp = aggregate(o, list(user), sum);
    if(any(temp[,1] != 1:nUsers)) stop("any(temp[,1] != 1:nUsers)");
    
    output$d1 = drop(output$d2 * alpha - 2*(temp[,2]/var_y + g0w/var_alpha));

    return(output);
}

# Compute the derivatives w.r.t. beta
# See output$d1 and output$d2
beta.derivatives.R <- function(user, item, y, xb, d0z, alpha, beta, u, v, var_y, var_beta){
    output = list();
    nItems = length(beta);

    # temp[j,2] is the number of observations that item j has
    temp = aggregate(rep(1,length(item)), list(item), sum);
    if(any(temp[,1] != 1:nItems)) stop("any(temp[,1] != 1:nItems)");
    
    output$d2 = 2*(temp[,2]/var_y + 1/var_beta);

    o = y - alpha[user] - xb - apply(u[user,,drop=FALSE] * v[item,,drop=FALSE], 1, sum);
    temp = aggregate(o, list(item), sum);
    if(any(temp[,1] != 1:nItems)) stop("any(temp[,1] != 1:nItems)");
    
    output$d1 = drop(output$d2 * beta - 2*(temp[,2]/var_y + d0z/var_beta));

    return(output);
}

# Compute the derivatives w.r.t. u
# See output$d1 and output$d2
#   (output$d2 is a 3D array, d2[i,,] is the Hessian for user i)
u.derivatives.R <- function(user, item, y, xb, Gw, alpha, beta, u, v, var_y, var_u){
    output = list();
    nUsers   = length(alpha);
    nFactors = ncol(v);

    obsIndex.forUser = tapply(1:length(user), list(user), c, simplify=F);
    if(any(names(obsIndex.forUser) != 1:nUsers)) stop("any(names(obsIndex.forUser) != 1:nUsers)");
    
    output$d1 = matrix(NA, nrow=nUsers, ncol=nFactors);
    output$d2 = array(NA, dim=c(nUsers, nFactors, nFactors));
    
    for(i in 1:nUsers){
        selectedObs   = obsIndex.forUser[[i]];
        selectedItems = item[selectedObs];
        v.selected = v[selectedItems,,drop=FALSE];
        sum.vv = t(v.selected) %*% v.selected;
        output$d2[i,,] = 2*(sum.vv / var_y + diag(1/var_u, nrow=nFactors));
        
        o.selected = y[selectedObs] - (alpha[i] + beta[selectedItems] + xb[selectedObs]);
        sum.ov = t(v.selected) %*% o.selected;
                
        output$d1[i,] = output$d2[i,,] %*% u[i,] - 2*(sum.ov / var_y + Gw[i,] / var_u);
    }
    return(output);
}

# Compute the derivatives w.r.t. v
# See output$d1 and output$d2
#   (output$d2 is a 3D array, d2[j,,] is the Hessian for item j)
v.derivatives.R <- function(user, item, y, xb, Dz, alpha, beta, u, v, var_y, var_v=1){
    output = list();
    nItems   = length(beta);
    nFactors = ncol(u);

    obsIndex.forItem = tapply(1:length(item), list(item), c, simplify=F);
    if(any(names(obsIndex.forItem) != 1:nItems)) stop("any(names(obsIndex.forItem) != 1:nItems)");
    
    output$d1 = matrix(NA, nrow=nItems, ncol=nFactors);
    output$d2 = array(NA, dim=c(nItems, nFactors, nFactors));
    
    for(j in 1:nItems){
        selectedObs   = obsIndex.forItem[[j]];
        selectedUsers = user[selectedObs];
        u.selected = u[selectedUsers,,drop=FALSE];
        sum.uu = t(u.selected) %*% u.selected;
        output$d2[j,,] = 2*(sum.uu / var_y + diag(1/var_v, nrow=nFactors));
        
        o.selected = y[selectedObs] - (alpha[selectedUsers] + beta[j] + xb[selectedObs]);
        sum.ou = t(u.selected) %*% o.selected;

        output$d1[j,] = output$d2[j,,] %*% v[j,] - 2*(sum.ou / var_y + Dz[j,] / var_v);
    }
    return(output);
}

obj.lineSearch <- function(
    eta,
    h.alpha, h.beta, h.u, h.v,
    user, item, y, x, w, z,
    alpha, beta, u, v,
    b, g0, G, d0, D,
    var_y, var_alpha, var_beta, var_u, var_v, debug=0, use.C=FALSE
){
    return(logLikelihood(
        user, item, y, x, w, z,
        alpha+eta*h.alpha, beta+eta*h.beta, u+eta*h.u, v+eta*h.v,
        b, g0, G, d0, D,
        var_y, var_alpha, var_beta, var_u, var_v, debug, use.C
    ));
}

restore.factors <- function(factors, nUsers, nItems, nFactors){

    if(length(factors) != (nUsers + nItems + nUsers*nFactors + nItems*nFactors)) stop("'factors' has length mismatch");
    
    return(list(
        alpha = factors[1:nUsers],
        beta  = factors[nUsers + (1:nItems)],
        u     = matrix(factors[nUsers + nItems + (1:(nUsers*nFactors))], nrow=nUsers, ncol=nFactors),
        v     = matrix(factors[(nUsers + nItems + (nUsers*nFactors) + 1):length(factors)], nrow=nItems, ncol=nFactors)
    ));
}

negLikelihood.gr <- function(
    factors, # c(alpha, beta, u, v);
    user, item, y, x, w, z,
    b, g0, G, d0, D,
    var_y, var_alpha, var_beta, var_u, var_v, debug=0
){

    f = restore.factors(factors, nrow(w), nrow(z), ncol(G));
    
    if(debug >= 1) check.input(user, item, y, x, w, z, f$alpha, f$beta, f$u, f$v, b, g0, G, d0, D, var_y, var_alpha, var_beta, var_u, var_v);
    
    xb  = x %*% b;
    g0w = w %*% g0;
    d0z = z %*% d0;
    Gw  = w %*% G;
    Dz  = z %*% D;
    
    d.alpha = alpha.derivatives.R(user, item, y, xb, g0w, f$alpha, f$beta, f$u, f$v, var_y, var_alpha);
    d.beta  = beta.derivatives.R(user, item, y, xb, d0z, f$alpha, f$beta, f$u, f$v, var_y, var_beta);
    d.u     = u.derivatives.R(user, item, y, xb, Gw, f$alpha, f$beta, f$u, f$v, var_y, var_u);
    d.v     = v.derivatives.R(user, item, y, xb, Dz, f$alpha, f$beta, f$u, f$v, var_y, var_v);
    
    return(c(d.alpha$d1, d.beta$d1, d.u$d1, d.v$d1));
}

negLikelihood <- function(
    factors, # c(alpha, beta, u, v);
    user, item, y, x, w, z,
    b, g0, G, d0, D,
    var_y, var_alpha, var_beta, var_u, var_v, debug=0
){
    f = restore.factors(factors, nrow(w), nrow(z), ncol(G));
    
    if(debug >= 1) check.input(user, item, y, x, w, z, f$alpha, f$beta, f$u, f$v, b, g0, G, d0, D, var_y, var_alpha, var_beta, var_u, var_v);

    return(-2*logLikelihood(user, item, y, x, w, z, f$alpha, f$beta, f$u, f$v, b, g0, G, d0, D, var_y, var_alpha, var_beta, var_u, var_v, debug));
}

# Gradient descent for optimizing alpha, beta, u and v
#   See output$alpha, output$beta, ...
GradDesc.R <- function(
    nSteps,
    user, item, y, x, w, z,
    alpha, beta, u, v,
    b, g0, G, d0, D,
    var_y, var_alpha, var_beta, var_u, var_v=1, method=3, eta=1, tol=1e-8, debug=0, verbose=0
){
    if(debug >= 1) check.input(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_y, var_alpha, var_beta, var_u, var_v);
    nObs     = length(y);
    nUsers   = length(alpha);
    nItems   = length(beta);
    nFactors = ncol(u);
    
    xb  = x %*% b;
    g0w = w %*% g0;
    d0z = z %*% d0;
    Gw  = w %*% G;
    Dz  = z %*% D;
    
    if(verbose >= 3){
        ll = logLikelihood(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_y, var_alpha, var_beta, var_u, var_v, debug);
        cat("    beginning logLikelihood = ",ll," + constant\n",sep="");
    }

    if(method == 3){
        trace = 0;
        if(verbose >= 6) trace = verbose - 5;
        par = c(alpha, beta, u, v);
        ans = optim(par=par, fn=negLikelihood, gr=negLikelihood.gr, method="CG", control=list(reltol=tol, maxit=nSteps, trace=trace),
                    user=user, item=item, y=y, x=x, w=w, z=z, b=b, g0=g0, G=G, d0=d0, D=D,
                    var_y=var_y, var_alpha=var_alpha, var_beta=var_beta, var_u=var_u, var_v=var_v, debug=debug);
        output = restore.factors(ans$par, nrow(w), nrow(z), ncol(G));
        if(debug >= 1) check.input(user, item, y, x, w, z, output$alpha, output$beta, output$u, output$v, b, g0, G, d0, D, var_y, var_alpha, var_beta, var_u, var_v);
        if(verbose >= 5){
            cat("    Optim output:\n");
            print(ans);
        }
        return(output);
    }
    
    if(method == 1){
        prev.h.alpha = 0;
        prev.h.beta  = 0;
        prev.h.u     = 0;
        prev.h.v     = 0;
        prev.alpha.d1 = 1;
        prev.beta.d1  = 1;
        prev.u.d1     = 1;
        prev.v.d1     = 1;
        prev.ll = -Inf;
    }
            
    for(s in 1:nSteps){
        
        d.alpha = alpha.derivatives.R(user, item, y, xb, g0w, alpha, beta, u, v, var_y, var_alpha);
        d.beta  = beta.derivatives.R(user, item, y, xb, d0z, alpha, beta, u, v, var_y, var_beta);
        d.u     = u.derivatives.R(user, item, y, xb, Gw, alpha, beta, u, v, var_y, var_u);
        d.v     = v.derivatives.R(user, item, y, xb, Dz, alpha, beta, u, v, var_y, var_v);
        
        if(method == 0){
            alpha = alpha - eta * d.alpha$d1;
            beta  = beta  - eta * d.beta$d1;
            u     = u     - eta * d.u$d1;
            v     = v     - eta * d.v$d1;
        }else if(method == 1){
            ##
            ## Conjugate gradient descent/ascent
            ##
            gamma = (sum((d.alpha$d1 - prev.alpha.d1)*d.alpha$d1) + 
                     sum((d.beta$d1  - prev.beta.d1) *d.beta$d1)  +
                     sum((d.u$d1 - prev.u.d1) * d.u$d1) +
                     sum((d.v$d1 - prev.v.d1) * d.v$d1)) /
                    (sum(prev.alpha.d1^2) + sum(prev.beta.d1^2) + sum(prev.u.d1^2) + sum(prev.v.d1^2));

            h.alpha = -d.alpha$d1 + gamma * prev.h.alpha;
            h.beta  = -d.beta$d1  + gamma * prev.h.beta;
            h.u = -d.u$d1 + gamma * prev.h.u;
            h.v = -d.v$d1 + gamma * prev.h.v;
            
            ans = optimize(f=obj.lineSearch, lower=0, upper=1, maximum=TRUE, tol=tol,
                h.alpha=h.alpha, h.beta=h.beta, h.u=h.u, h.v=h.v,
                user=user, item=item, y=y, x=x, w=w, z=z,
                alpha=alpha, beta=beta, u=u, v=v,
                b=b, g0=g0, G=G, d0=d0, D=D,
                var_y=var_y, var_alpha=var_alpha, var_beta=var_beta, var_u=var_u, var_v=var_v
            );
            
            eta = ans$maximum;
            if(ans$objective <= prev.ll){
                if(verbose >= 2) cat("    Early stopping: eta=",eta,"   new logLikelihood=",ans$objective,"\n",sep="");
                break;
            }
            if(verbose >= 4){
                cat("      eta=",eta,"\n",
                    "      objective=",ans$objective,"\n",sep="");
                if(verbose >= 10){
                    cat("      h.alpha=\n"); print(h.alpha);
                    cat("      h.beta=\n");  print(h.beta);
                    cat("      h.u=\n");     print(h.u);
                    cat("      h.v=\n");     print(h.v);
                }
            }
            
            alpha = alpha + eta*h.alpha;
            beta  = beta  + eta*h.beta;
            u     = u     + eta*h.u;
            v     = v     + eta*h.v;

            prev.h.alpha = h.alpha;
            prev.h.beta  = h.beta;
            prev.h.u     = h.u;
            prev.h.v     = h.v;
            prev.alpha.d1 = d.alpha$d1;
            prev.beta.d1  = d.beta$d1;
            prev.u.d1     = d.u$d1;
            prev.v.d1     = d.v$d1;
            prev.ll     = ans$objective;
            
        }else if(method == 2){
            alpha = alpha - eta * (1/d.alpha$d2) * d.alpha$d1;
            beta  = beta  - eta * (1/d.beta$d2)  * d.beta$d1;
            for(i in 1:nUsers) u[i,] = drop(u[i,] - eta * solve(d.u$d2[i,,]) %*% d.u$d1[i,]);
            for(j in 1:nItems) v[j,] = drop(v[j,] - eta * solve(d.v$d2[j,,]) %*% d.v$d1[j,]);
        }
        if(debug >= 1) check.input(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_y, var_alpha, var_beta, var_u, var_v);
        if(verbose >= 3){
            ll = logLikelihood(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_y, var_alpha, var_beta, var_u, var_v, debug);
            cat("    step ",s,": logLikelihood = ",ll," + constant\n",sep="");
        }
    }
    
    output = list(alpha=alpha, beta=beta, u=u, v=v);
    return(output);
}

###
### Other utility functions
###
compare.two.lists <- function(list1, list2){
    name1 = names(list1);
    name2 = names(list2);
    
    name.all = union(name1, name2);
    num = length(name.all);
    output = data.frame(row.names=name.all, MaxAbsErr=rep(NA,num), MeanAbsErr=rep(NA,num), List1_MinAbs=rep(NA,num), List1_MaxAbs=rep(NA,num), List2_MinAbs=rep(NA,num), List2_MaxAbs=rep(NA,num));
    
    for(c in name.all){
        if(c %in% name1){
            output[c,'List1_MinAbs'] = min(abs(list1[[c]]));
            output[c,'List1_MaxAbs'] = max(abs(list1[[c]]));
        }
        if(c %in% name2){
            output[c,'List2_MinAbs'] = min(abs(list2[[c]]));
            output[c,'List2_MaxAbs'] = max(abs(list2[[c]]));
        }
        if(c %in% name1 && c %in% name2){
            output[c,'MaxAbsErr'] = max(abs(list1[[c]] - list2[[c]]));
            output[c,'MeanAbsErr'] = mean(abs(list1[[c]] - list2[[c]]));
        }
    }
    return(output);
}

