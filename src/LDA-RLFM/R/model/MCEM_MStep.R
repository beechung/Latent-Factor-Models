### Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

# Use bayesglm to fit the coefficients (g0, d0 or c0) for the main effects (alpha, beta, gamma)
#   E.g., output = fit.forMainEffect.bayesglm(alpha, x_user);
#   See output$coef (coefficient vector) and output$rss (residual sum of squares)
fit.forMainEffect.bayesglm <- function(
    target, feature, lm=F,...
){
    if(length(target) != nrow(feature)) stop("length(target) != nrow(feature)");
    if(all(feature == 0)){
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
    target, feature, lm=F,...
){
    if(nrow(target) != nrow(feature)) stop("nrow(target) != nrow(feature)");
    nFeatures = ncol(feature);
    nFactors  = ncol(target);
    output = list();
    output$coef = matrix(NA, nrow=nFeatures, ncol=nFactors);
    rss = 0;

    for(f in 1:nFactors){
        if(all(feature == 0)){
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
#   list(b, g0, d0, c0, G, D, H, var_y, var_alpha, var_beta, var_gamma, var_u, var_v, var_s, lambda, eta);
#
# NOTE: Do regression for FACTOR only if param$var_FACTOR is not NULL.
#
MCEM_MStep <- function(
    factor.mean, factor.sumvar, obs, corpus, feature, param, try, 
    lda.objval, lm=F,
    debug=0, verbose=0, ...
){
    size = syncheck.LDA_RLFM.spec(factor=factor.mean, obs=obs, corpus=corpus, feature=feature, param=param,
                                  is.corpus_topic.matrix=if(!is.null(factor.mean$corpus_topic)) is.matrix(factor.mean$corpus_topic) else FALSE);
    
    nObs     = size$nObs;
    nUsers   = size$nUsers;
    nItems   = size$nItems;
    nFactors = size$nFactors;
    nTopics  = size$nTopics;

    # determine b and var_y
    b.time = proc.time();
    weight = factor.mean$gamma2[obs$user];
    target = factor.mean$o_gamma / weight;
    
    if(all(feature$x_dyad == 0)){
        fit = list(fitted.values=0, coefficients=rep(0,ncol(feature$x_dyad)));
    }else{
        fit = lm(target ~ feature$x_dyad - 1, weights=weight, model=FALSE);
    }
    if(length(fit$coef) != ncol(feature$x_dyad)) stop("length(fit$coef) != ncol(feature$x_dyad)");
    rss = sum((fit$fitted.values - target)^2);
    
    param$b = fit$coefficients;
    param$var_y = (factor.sumvar$o_adj + weight * rss) / nObs;

    time.used = proc.time() - b.time;
    if(verbose > 0){
        LL = logLikelihood(obs, factor.mean, feature, param, corpus);
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
            LL = logLikelihood(obs, factor.mean, feature, param, corpus);
            cat("after fitting g0: complete data logLikelihood = ",LL,",  used ",time.used[3]," sec\n", sep="");
        }
    }
    # determin d0 and var_beta
    if(!is.null(param$var_beta)){
        b.time = proc.time();
        fit = fit.forMainEffect.bayesglm(factor.mean$beta, feature$x_item, lm=lm,...);
        param$d0 = fit$coef;
        param$var_beta = (fit$rss + factor.sumvar$beta) / nItems;
        time.used = proc.time() - b.time;
        if(verbose > 0){
            LL = logLikelihood(obs, factor.mean, feature, param, corpus);
            cat("after fitting d0: complete data logLikelihood = ",LL,",  used ",time.used[3]," sec\n", sep="");
        }
    }

    # determin c0 and var_gamma
    if(!is.null(param$var_gamma)){
        b.time = proc.time();
        fit = fit.forMainEffect.bayesglm(factor.mean$gamma, feature$x_user, lm=lm, ...);
        param$c0 = fit$coef;
        param$var_gamma = (fit$rss + factor.sumvar$gamma) / nUsers;
        time.used = proc.time() - b.time;
        if(verbose > 0){
            LL = logLikelihood(obs, factor.mean, feature, param, corpus);
            cat("after fitting c0: complete data logLikelihood = ",LL,",  used ",time.used[3]," sec\n", sep="");
        }
    }
    
    if(nFactors > 0){
        # determin G and var_u
        if(!is.null(param$var_u)){
            b.time = proc.time();
            fit = fit.forFactors.bayesglm(factor.mean$u, feature$x_user, lm=lm,...);
            param$G = fit$coef;
            param$var_u = (fit$rss + factor.sumvar$u) / (nUsers * nFactors);
            time.used = proc.time() - b.time;
            if(verbose > 0){
                LL = logLikelihood(obs, factor.mean, feature, param, corpus);
                cat("after fitting G : complete data logLikelihood = ",LL,",  used ",time.used[3]," sec\n", sep="");
            }
        }
        # determin D and var_v
        if(!is.null(param$var_v)){
            b.time = proc.time();
            fit = fit.forFactors.bayesglm(factor.mean$v, feature$x_item, lm=lm, ...);
            param$D = fit$coef;
            param$var_v = (fit$rss + factor.sumvar$v) / (nItems * nFactors);
            time.used = proc.time() - b.time;
            if(verbose > 0){
                LL = logLikelihood(obs, factor.mean, feature, param, corpus);
                cat("after fitting D : complete data logLikelihood = ",LL,",  used ",time.used[3]," sec\n", sep="");
            }
        }
    }
    
    if(nTopics > 0){
        # determin H and var_s
        if(!is.null(param$var_s)){
            b.time = proc.time();
            fit = fit.forFactors.bayesglm(factor.mean$s, feature$x_user, lm=lm,...);
            param$H = fit$coef;
            param$var_s = (fit$rss + factor.sumvar$s) / (nUsers * nTopics);
            time.used = proc.time() - b.time;
            if(verbose > 0){
                LL = logLikelihood(obs, factor.mean, feature, param, corpus);
                cat("after fitting H : complete data logLikelihood = ",LL,",  used ",time.used[3]," sec\n", sep="");
            }
        }
        
        if(verbose > 0){
            b.time = proc.time();
            count = getTopicCounts(corpus, factor.mean$corpus_topic, nItems, nTopics, size$nTerms);
            time.used = proc.time() - b.time;
            cat("    call getTopicCounts: used ",time.used[3]," sec\n");
        }
        
        # determine eta
        if(length(try$eta) > 0){
            index = which.min(lda.objval$eta);
            param$eta = try$eta[index];
            if(verbose > 0){
                print(data.frame(
                    eta=try$eta, objval=lda.objval$eta, 
                    CD.negLL=compute_LDA_negLL(try$eta, nTopics, size$nTerms, count$cnt_topic, count$cnt_topic_term)
                ));
                LL = logLikelihood(obs, factor.mean, feature, param, corpus);
                cat("after fitting    eta: complete data logLikelihood = ",LL," (eta=",param$eta,")\n", sep="");
            }
        }
        
        # determine lambda
        if(length(try$lambda) > 0){
            index = which.min(lda.objval$lambda);
            param$lambda = try$lambda[index];
            if(verbose > 0){
                print(data.frame(
                    lambda=try$lambda, objval=lda.objval$lambda, 
                    CD.negLL=compute_LDA_negLL(try$lambda, nItems, nTopics, count$cnt_item, count$cnt_item_topic)
                ));
                LL = logLikelihood(obs, factor.mean, feature, param, corpus);
                cat("after fitting lambda: complete data logLikelihood = ",LL," (lambda=",param$lambda,")\n", sep="");
            }
        }
    }
    
    return(param);
}
