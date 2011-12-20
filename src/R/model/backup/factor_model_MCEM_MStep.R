### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen


# Use bayesglm to fit the coefficients (g0, d0 or c0) for the main effects (alpha, beta, gamma)
#   E.g., output = fit.forMainEffect.bayesglm(alpha, x_user);
#   See output$coef (coefficient vector) and output$rss (residual sum of squares)
fit.forMainEffect.bayesglm <- function(
    target, feature, set.to.zero=F, lm=F,...
){
    if(length(target) != nrow(feature)) stop("length(target) != nrow(feature)");
    if(set.to.zero || all(feature == 0)){
        fit = list(coefficients=rep(0,ncol(feature)),fitted.values=0);
	}else if(ncol(feature)==1 && all(feature == 1)){
		mean.target = mean(target);
		fit = list(coefficients=mean.target,fitted.values=mean.target);
	}else{
        fit = if(!lm) bayesglm(target ~ feature -1, model=FALSE, x=FALSE, y=FALSE, ...) 
              else    lm(target ~ feature -1, model=FALSE, x=FALSE, y=FALSE, ...);
    }
    if(length(fit$coef) != ncol(feature)) stop("length(fit$coef) != ncol(feature)");
    output = list();
    output$coef = fit$coefficients;
    output$rss  = sum((fit$fitted.values - target)^2);
    return(output);
}

# Use bayesglm to fit the coefficients (G, D, H) for the factors (u, v, s)
#   E.g., output = fit.forFactors.bayesglm(u, x_user);
#   See output$coef (coefficient matrix) and output$rss (sum of squared errors)
#   Note: target is a matrix
fit.forFactors.bayesglm <- function(
    target, feature, set.to.zero=FALSE, lm=FALSE,...
){
    if(nrow(target) != nrow(feature)) stop("nrow(target) != nrow(feature)");
    nFeatures = ncol(feature);
    nFactors  = ncol(target);
    output = list();
    output$coef = matrix(NA, nrow=nFeatures, ncol=nFactors);
    rss = 0;

    for(f in 1:nFactors){
        if(set.to.zero || all(feature == 0)){
            fit = list(coefficients=rep(0,ncol(feature)),fitted.values=0);
		}else if(ncol(feature)==1 && all(feature == 1)){
			mean.target = mean(target[,f]);
			fit = list(coefficients=mean.target,fitted.values=mean.target);
		}else{
            fit = if(!lm) bayesglm(target[,f] ~ feature -1, model=FALSE, x=FALSE, y=FALSE, ...) 
                  else    lm(target[,f] ~ feature -1, model=FALSE, x=FALSE, y=FALSE, ...);
        }
        if(length(fit$coef) != ncol(feature)) stop("length(fit$coef) != ncol(feature)");
        output$coef[,f] = fit$coefficients;
        rss = rss + sum((fit$fitted.values - target[,f])^2);
    }
    
    output$rss = rss;
    return(output);
}

#
# Monte-Carlo EM (M-step)
# OUTPUT: 
#   list(b, g0, d0, G, var_y, var_alpha, var_beta, var_v);
#
# NOTE: Do regression for FACTOR only if param$var_FACTOR is not NULL.
#
MCEM_MStep <- function(
    factor.mean, factor.sumvar, obs, feature, param, fit.var_v, is.G.zero, is.d0.zero,
	lm=F, debug=0, verbose=0, ...
){
    size = syncheck.factorModel.spec(factor=factor.mean, obs=obs, feature=feature, param=param);
    
    nObs     = size$nObs;
    nUsers   = size$nUsers;
    nFactors = size$nFactors;

    # determine b and var_y
    b.time = proc.time();
    target = factor.mean$fErr;
    
    if(all(feature$x_dyad == 0)){
        fit = list(fitted.values=0, coefficients=rep(0,ncol(feature$x_dyad)));
    }else{
        fit = lm(target ~ feature$x_dyad - 1, model=FALSE);
    }
    if(length(fit$coef) != ncol(feature$x_dyad)) stop("length(fit$coef) != ncol(feature$x_dyad)");
    rss = sum((fit$fitted.values - target)^2);
    
    param$b = fit$coefficients;
    param$var_y = (factor.sumvar$fErr + rss) / nObs;

    time.used = proc.time() - b.time;
    if(verbose > 0){
        LL = logLikelihood(obs, factor.mean, feature, param, is.logistic=FALSE);
        cat("after fitting b : complete data logLikelihood = ",LL,",  used ",time.used[3]," sec\n", sep="");
    }
    
    # determin g0 and var_alpha
    if(!is.null(param$var_alpha)){
        b.time = proc.time();
        
		fit = fit.forMainEffect.bayesglm(factor.mean$alpha, feature$x_user, lm=lm, ...);
		param$g0 = fit$coef;
        param$var_alpha = (fit$rss + factor.sumvar$alpha) / nUsers;
		
        time.used = proc.time() - b.time;
        if(verbose > 0){
            LL = logLikelihood(obs, factor.mean, feature, param, is.logistic=FALSE);
            cat("after fitting g0: complete data logLikelihood = ",LL,",  used ",time.used[3]," sec\n", sep="");
        }
    }
    # determin d0 and var_beta
    if(!is.null(param$var_beta)){
        b.time = proc.time();
        fit = fit.forMainEffect.bayesglm(factor.mean$beta, feature$x_user, set.to.zero=is.d0.zero, lm=lm,...);
        param$d0 = fit$coef;
        param$var_beta = (fit$rss + factor.sumvar$beta) / nUsers;
        time.used = proc.time() - b.time;
        if(verbose > 0){
            LL = logLikelihood(obs, factor.mean, feature, param, is.logistic=FALSE);
            cat("after fitting d0: complete data logLikelihood = ",LL,",  used ",time.used[3]," sec\n", sep="");
        }
    }

    if(nFactors > 0){
        # determin G and var_v
        if(!is.null(param$var_v)){
            b.time = proc.time();

			fit = fit.forFactors.bayesglm(factor.mean$v, feature$x_user, set.to.zero=is.G.zero, lm=lm,...);
	        param$G = fit$coef;

			if(fit.var_v){
            	param$var_v = (fit$rss + factor.sumvar$v) / (nUsers * nFactors);
			}
            time.used = proc.time() - b.time;
            if(verbose > 0){
                LL = logLikelihood(obs, factor.mean, feature, param, is.logistic=FALSE);
                cat("after fitting G : complete data logLikelihood = ",LL,",  used ",time.used[3]," sec\n", sep="");
            }
        }
    }
    
    return(param);
}
