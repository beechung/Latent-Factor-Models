### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
###
### Author: Liang Zhang

# Use bayesglm to fit the coefficients (g0 or d0) for the main effects (alpha or beta)
#   E.g., output = fit.forMainEffect.bayesglm(alpha, w);
#   See output$coef (coefficient vector) and output$rss (residual sum of squares)
fit.forMainEffect.bayesglm <- function(
    target, feature, lm=F,...
){
    if(length(target) != nrow(feature)) stop("length(target) != nrow(feature)");
        flagSingleCol=0;
    if(all(feature == 0)){
        fit = list(coefficients=rep(0,ncol(feature)),fitted.values=0);
    } else{
           if( ncol(feature)==1 ) {
                        flagSingleCol=1;
                        feature=cbind(feature,matrix(0,nrow=length(feature),ncol=1)); #workaround bayesglm bug
           }
        fit = if(!lm)bayesglm(target ~ feature -1, model=FALSE, x=FALSE, y=FALSE, ...) else lm(target ~ feature -1, model=FALSE, x=FALSE, y=FALSE, ...);;
    }

    if(length(fit$coef) != ncol(feature)) stop("length(fit$coef) != ncol(feature)");
        output = list();
        if(flagSingleCol==1){
                output$coef = fit$coefficients[1];
        } else{
        output$coef = fit$coefficients;
        }
        output$rss  = sum((fit$fitted.values - target)^2);
    return(output);
}

# Use bayesglm to fit the coefficients (G or D) for the factors (u or v)
#   E.g., output = fit.forFactors.bayesglm(u, w);
#   See output$coef (coefficient matrix) and output$rss (standard error)
#   Note: target is a matrix
fit.forFactors.bayesglm <- function(
    target, feature, lm=F,...
){
    if(nrow(target) != nrow(feature)) stop("nrow(target) != nrow(feature)");
    nFeatures = ncol(feature);
    nFactors  = ncol(target);
        flagSingleCol=0;
    output = list();
    output$coef = matrix(NA, nrow=nFeatures, ncol=nFactors);
    rss = rep(0,nFactors);

        if( ncol(feature)==1 ) {
                flagSingleCol=1;
                feature=cbind(feature,matrix(0,nrow=length(feature),ncol=1)); #workaround bayesglm bug
        }

        for(f in 1:nFactors){
        if(all(feature == 0)){
            fit = list(coefficients=rep(0,ncol(feature)),fitted.values=0);
        }else{

            fit = if(!lm)bayesglm(target[,f] ~ feature -1, model=FALSE, x=FALSE, y=FALSE, ...) else lm(target[,f] ~ feature -1, model=FALSE, x=FALSE, y=FALSE, ...);
        }
        if(length(fit$coef) != ncol(feature)) stop("length(fit$coef) != ncol(feature)");
                if(flagSingleCol==1){
                  output$coef[,f] = fit$coefficients[1];
                } else{
                output$coef[,f] = fit$coefficients;
                }
        rss[f] = rss[f] + sum((fit$fitted.values - target[,f])^2);
    }

    output$rss = rss;
    return(output);
}

fit.forMainEffect.glmnet<- function(
    target, feature, ...
)
{
                cat("MainEffect for ",dim(target), "and feature dim = ", dim(feature) ,"\n");
            if(length(target) != nrow(feature)) stop("length(target) != nrow(feature)");
#           if(!is.matrix(feature)) stop("feature has to be a matrix");
            #cat("now run glmnet\n");
            if (sum(abs(feature[,1]-1))==0)
            {
                        if (ncol(feature)==1 ) {
                                output=fit.forMainEffect.bayesglm(target,as.matrix(feature), ...);
                                return(output);
                        } else {
                                x = Matrix(feature[,2:ncol(feature)],sparse=T);
                        }
                } else {
                        stop("feature should have a constant term for intercept");
            }

            fit0 = cv.glmnet(x,target,nfolds=3,type="mse");
            #cat("lambda=",fit0$lambda.1se,"\n");
            #fit = glmnet(x,target,family="gaussian",lambda=fit0$lambda.min);
            lambdaind = which(fit0$lambda==fit0$lambda.min);
            a0 = fit0$glmnet.fit$a0[lambdaind];
            beta = fit0$glmnet.fit$beta[,lambdaind];
            output = list();
            output$coef = as.vector(c(a0,beta));
            fitted.values = feature%*%output$coef;
            output$rss = sum((fitted.values - target)^2);

                return(output);
}
fit.forFactors.glmnet <- function(
    target, feature, ...
){
    cat("FactorEffect for ",dim(target), "and feature dim = ", dim(feature) ,"\n");
    if(nrow(target) != nrow(feature)) stop("nrow(target) != nrow(feature)");
#if(!is.matrix(feature)) stop("feature has to be a matrix");
    if (sum(abs(feature[,1]-1))==0)
    {
                if (ncol(feature)==1 ) {
                        output=fit.forFactors.bayesglm(target=target,feature=as.matrix(feature));
                        return(output);
                }else{
                        x = Matrix(feature[,2:ncol(feature)],sparse=T);
                }
        } else {
                stop("feature should has a constant term for intercept");
    }

           nFeatures = ncol(feature);
           nFactors  = ncol(target);
           output = list();
           output$coef = matrix(NA, nrow=nFeatures, ncol=nFactors);
           rss = rep(0,nFactors);

           for(f in 1:nFactors){
                fit0 = cv.glmnet(x,target[,f],nfolds=3,type="mse");
           #fit = glmnet(x,target[,f],family="gaussian",lambda=fit0$lambda.min);
           lambdaind = which(fit0$lambda==fit0$lambda.min);
               a0 = fit0$glmnet.fit$a0[lambdaind];
               beta = fit0$glmnet.fit$beta[,lambdaind];
               output$coef[,f] = as.vector(c(a0,beta));
           fitted.values = feature%*%output$coef[,f];
               rss[f] = rss[f] + sum((fitted.values - target[,f])^2);
           }
           output$rss = rss;
           return(output);
}


# Monte-Carlo EM (M-step)
#   See output$b, output$var_y, ...
MC_MStep <- function(
    user, item, y, x, w, z,
    o, o.sumvar, alpha, alpha.sumvar, beta, beta.sumvar, u, u.sumvar, v, v.sumvar,
        bfixed=FALSE,
        b_old=NULL, # not used unless bfixed=T
    debug=0, lm=F, use.glmnet=F, ...
){
    nObs     = length(y);
    nUsers   = length(alpha);
    nItems   = length(beta);
    nFactors = ncol(u);

    output = list();

    # determine b and var_y
    target = y - o;
        rss=0;
        if(bfixed == FALSE){
                 if(all(x == 0)){
                     fit = list(fitted.values=0, coefficients=rep(0,ncol(x)));
                 }else{
                     fit = lm(target ~ x -1, model=FALSE);
                 }
        if(length(fit$coef) != ncol(x)) stop("length(fit$coef) != ncol(x)");
        rss = sum((fit$fitted.values - target)^2);

        output$b = fit$coefficients;
        }       else{
                if(is.null(b_old)) { error("No b_old provided in MC_MStep when bfixed=T\n");}
                output$b=b_old;
        rss = sum(( (x %*% b_old) - target)^2);
        }

        output$var_y = (rss + o.sumvar) / nObs;

    cat("use.glmnet=",use.glmnet,"\n");
    # determin g0 and var_alpha
    if (use.glmnet==F )
    {
        fit = fit.forMainEffect.bayesglm(alpha, w, lm=lm, ...);
    } else
    {
        fit = fit.forMainEffect.glmnet(alpha, w, ...);
    }
    output$g0 = fit$coef;
    output$var_alpha = (sum(fit$rss) + alpha.sumvar) / nUsers;

    # determin d0 and var_beta
    if (use.glmnet==F)
    {
        fit = fit.forMainEffect.bayesglm(beta, z, lm=lm,...);
    } else
    {
        fit = fit.forMainEffect.glmnet(beta, z, ...);
    }
    output$d0 = fit$coef;
    output$var_beta = (sum(fit$rss) + beta.sumvar) / nItems;

    # determin G and var_u
    if (use.glmnet==F)
    {
        fit = fit.forFactors.bayesglm(u, w, lm=lm,...);
    } else
    {
        fit = fit.forFactors.glmnet(u, w, ...);
    }
    output$G = fit$coef;
        if(var_u_fixed==FALSE){
            output$var_u = (sum(fit$rss) + sum(u.sumvar)) / (nUsers * nFactors);
        } else{
                if(use_nFactors_based_fixed_u_var==TRUE){
                        output$var_u = C * 6/sqrt(nFactors);
                } else{
            output$var_u = forever_fixed_u_var;
        }
            cat("var_u is fixed, = ",output$var_u,"\n",sep="");
        }

    # determin D and var_v
    if (use.glmnet==F)
    {
        fit = fit.forFactors.bayesglm(v, z, lm=lm, ...);
    } else
    {
        fit = fit.forFactors.glmnet(v, z, ...);
    }
    output$D = fit$coef;

        if(var_v_fixed==FALSE){
            output$var_v = (sum(fit$rss) + sum(v.sumvar)) / (nItems * nFactors);
        } else{
                if(use_nFactors_based_fixed_v_var==TRUE){
                  output$var_v = C * 6/sqrt(nFactors);
                } else{
                        output$var_v = forever_fixed_v_var;
                }
                cat("var_v is fixed, = ",output$var_u,"\n",sep="");
        }

    return(output);
}

MC_MStep_logistic <- function(
    user, item, y, x, w, z,
    o, alpha, alpha.sumvar, beta, beta.sumvar, u, u.sumvar, v, v.sumvar,
    debug=0, lm=F, use.glmnet=F, ...
){
    nObs     = length(y);
    nUsers   = length(alpha);
    nItems   = length(beta);
    nFactors = ncol(u);

    output = list();

    # determine b and var_y
    if(all(x == 0)){
        fit = list(fitted.values=0, coefficients=rep(0,ncol(x)));
    }else{
        fit = glm(y ~ x -1, family=binomial(link = "logit"),offset=o, model=F);
    }
    if(length(fit$coef) != ncol(x)) stop("length(fit$coef) != ncol(x)");

    output$b = fit$coefficients;

    cat("use.glmnet=",use.glmnet,"\n");
    # determin g0 and var_alpha
    if (use.glmnet==F)
    {
        fit = fit.forMainEffect.bayesglm(alpha, w, lm=lm, ...);
    } else
    {
        fit = fit.forMainEffect.glmnet(alpha, w, ...);
    }
    output$g0 = fit$coef;
    output$var_alpha = (sum(fit$rss) + alpha.sumvar) / nUsers;

    # determin d0 and var_beta
    if (use.glmnet==F)
    {
        fit = fit.forMainEffect.bayesglm(beta, z, lm=lm,...);
    } else
    {
        fit = fit.forMainEffect.glmnet(beta, z, ...);
    }
    output$d0 = fit$coef;
    output$var_beta = (sum(fit$rss) + beta.sumvar) / nItems;

    # determin G and var_u
    if (use.glmnet==F)
    {
        fit = fit.forFactors.bayesglm(u, w, lm=lm,...);
    } else
    {
        fit = fit.forFactors.glmnet(u, w, ...);
    }
    output$G = fit$coef;
    output$var_u = (sum(fit$rss) + sum(u.sumvar)) / (nUsers * nFactors);

    # determin D and var_v
    if (use.glmnet==F)
    {
        fit = fit.forFactors.bayesglm(v, z, lm=lm, ...);
    } else
    {
        fit = fit.forFactors.glmnet(v, z, ...);
    }
    output$D = fit$coef;
    output$var_v = (sum(fit$rss) + sum(v.sumvar)) / (nItems * nFactors);
        cat("fit.rss=",sum(fit$rss)," v.sumvar=",sum(v.sumvar)," var_v=",output$var_v,"\n");
    return(output);
}


MC_MStep_factoronly <- function(
    user, item, y, x, w, z,
    o, o.sumvar, alpha, alpha.sumvar, beta, beta.sumvar, u, u.sumvar, v, v.sumvar,
    debug=0, lm=F, use.glmnet=F, ...
){
    nObs     = length(y);
    nUsers   = length(alpha);
    nItems   = length(beta);
    nFactors = ncol(u);

    output = list();

    # determine b and var_y
    target = y - o;
    if(all(x == 0)){
        fit = list(fitted.values=0, coefficients=rep(0,ncol(x)));
    }else{
        fit = lm(target ~ x -1, model=FALSE);
    }
    if(length(fit$coef) != ncol(x)) stop("length(fit$coef) != ncol(x)");
    rss = sum((fit$fitted.values - target)^2);

    output$b = fit$coefficients;
    output$var_y = (rss + o.sumvar) / nObs;

    # determin G and var_u
    if (use.glmnet==F)
    {
        fit = fit.forFactors.bayesglm(u, w, lm=lm,...);
    } else
    {
        fit = fit.forFactors.glmnet(u, w, ...);
    }
    output$G = fit$coef;
    output$var_u = (sum(fit$rss) + sum(u.sumvar)) / (nUsers * nFactors);

    # determin D and var_v
    if (use.glmnet==F)
    {
        fit = fit.forFactors.bayesglm(v, z, lm=lm, ...);
    } else
    {
        fit = fit.forFactors.glmnet(v, z, ...);
    }
    output$D = fit$coef;
    output$var_v = (sum(fit$rss) + sum(v.sumvar)) / (nItems * nFactors);

    return(output);
}

# ICM (regression)
#   See output$b, output$var_y, ...
ICM_Regression <- function(
    user, item, y, x, w, z,
    alpha, beta, u, v,
    debug=0, lm=F,...
){
    nObs     = length(y);
    nUsers   = length(alpha);
    nItems   = length(beta);
    nFactors = ncol(u);

    output = list();

    # determine b and var_y
    target = y - alpha[user] - beta[item] - apply(u[user,,drop=FALSE] * v[item,,drop=FALSE], 1, sum);
    if(all(x == 0)){
        fit = list(fitted.values=0, coefficients=rep(0,ncol(x)));
    }else{
        fit = lm(target ~ x -1, model=FALSE);
    }
    if(length(fit$coef) != ncol(x)) stop("length(fit$coef) != ncol(x)");
    rss = sum((fit$fitted.values - target)^2);

    output$b = fit$coefficients;
    output$var_y = rss / nObs;

    # determin g0 and var_alpha
    fit = fit.forMainEffect.bayesglm(alpha, w,lm=lm, ...);
    output$g0 = fit$coef;
    output$var_alpha = sum(fit$rss) / nUsers;

    # determin d0 and var_beta
    fit = fit.forMainEffect.bayesglm(beta, z,lm=lm, ...);
    output$d0 = fit$coef;
    output$var_beta = sum(fit$rss) / nItems;

    # determin G and var_u
    fit = fit.forFactors.bayesglm(u, w, lm=lm,...);
    output$G = fit$coef;
    output$var_u = sum(fit$rss) / (nUsers * nFactors);

    # determin D and var_v
    fit = fit.forFactors.bayesglm(v, z, lm=lm,...);
    output$D = fit$coef;
    output$var_v = sum(fit$rss) / (nItems * nFactors);

    return(output);
}

MC_MStep_logistic_arscid <- function(
    user, item, y, x, b, w, z,
    o, alpha, alpha.sumvar, beta, beta.sumvar, u, u.sumvar, v, v.sumvar,
    debug=0, lm=F, use.glmnet=F, fit.ars.alpha=F, fit.regression=T,
    beta.int=T, main.effects=F,...
){
    nObs     = length(y);
    nUsers   = length(alpha);
    nItems   = length(beta);
    nFactors = ncol(u);

    output = list();

    # find ars alpha ...
    if (fit.ars.alpha )
      {
        ars_alpha = estalpha(y,b,o);
        output$ars_alpha = ars_alpha;
      }
    else
      {
        output$ars_alpha = 0.5;
        ars_alpha = 0.5;
      }

    # determine b and var_y
    if(beta.int == F){
      if(all(x == 0)){
        fit = list(fitted.values=0, coefficients=rep(0,ncol(x)));
      }else{
        #fit = glm(y ~ x -1, family=binomial(link = "logit"),offset=o, model=F);
        # fit as covariate
        nobs = length(x)
        x = cbind(matrix(x,length(x),1), rep(0,length(x)))

        fit = bayesglm(y ~ x - 1, family=binomial(link="logit"), offset=o ,model=F, prior.scale = 5);

        fit$coef = fit$coef[1]; fit$coefficients=fit$coefficients[1]; x = matrix(x[,1],nobs,1)
      }
      #if(length(fit$coef) != ncol(x)) stop("length(fit$coef) != ncol(x)");
      output$b = fit$coefficients;
    } else {
      #fit in random effects heirarchy inside of beta ... but still center for computation ...
      output$b = b;
    }

    if (fit.ars.alpha && ncol(as.matrix(x))>1)
    {
        stop("Currently fit.ars.alpha=T only works for ncol(x)==1");
    }
    #cat("use.glmnet=",use.glmnet,"\n");
    cat("fit regression=",fit.regression,"\n")
    cat("intercept in beta prior =",beta.int,"\n")

    #If fit.regression=T do the normal thing, if not just find vars
    if(fit.regression){
    # determin g0 and var_alpha
      if (use.glmnet==F  )
        {
          fit = fit.forMainEffect.bayesglm(alpha, w, lm=lm,...);
        } else
      {
        fit = fit.forMainEffect.glmnet(alpha, w, ...);
      }
      output$g0 = fit$coef;
      output$var_alpha = (sum(fit$rss) + alpha.sumvar) / nUsers;

    # determin d0 and var_beta ( and b if in heirarchy )
      if (beta.int) z2 = cbind(1,z) else z2=z

      if (use.glmnet==F    )
        {
          fit = fit.forMainEffect.bayesglm(beta, z2, lm=lm,...);
        } else
      {
        fit = fit.forMainEffect.glmnet(beta, z2, ...);
      }

      output$d0 = fit$coef;
      output$var_beta = (sum(fit$rss) + beta.sumvar) / nItems;
      #output$var_beta = 1;

      if(!main.effects){
                                        # determin G and var_u
        if (use.glmnet==F)
          {
            fit = fit.forFactors.bayesglm(u, w, lm=lm,...);
          } else
        {
          fit = fit.forFactors.glmnet(u, w, ...);
        }
        output$G = fit$coef;
	
	#cat("var_u.rss=",fit$rss,"\n");
	#cat("u.sumvar=",u.sumvar,"\n");
	#output$var_u = (fit$rss + u.sumvar)/nUsers;
        #output$var_u = (fit$rss + u.sumvar) / (nUsers * nFactors);
        #output$var_u = 1;
	output$var_u = rep(1,nFactors);

        # determin D and var_v
        if (use.glmnet==F)
          {
            fit = fit.forFactors.bayesglm(v, z, lm=lm, ...);
          } else
        {
          fit = fit.forFactors.glmnet(v, z, ...);
        }
        output$D = fit$coef;
	cat("var_v.rss=",fit$rss,"\n");
        cat("v.sumvar=",v.sumvar,"\n");
	output$var_v = (fit$rss + v.sumvar)/nItems;
        #output$var_v = (fit$rss + v.sumvar) / (nItems * nFactors);
        #output$var_v = rep(1,nFactors);

      }
    } else
    {
      output$var_alpha = (sum(alpha^2) + alpha.sumvar) / nUsers;
      #output$var_beta = 1;
      output$var_beta = (sum(beta^2) + beta.sumvar) / nItems;
      #if(!main.effects) output$var_u = (sum(u^2) + u.sumvar) / (nUsers * nFactors);
      #output$var_v = rep(1,nFactors);
      #if (!main.effects) output$var_u = (apply(u^2,2,sum)+u.sumvar)/nItems;
      output$var_u = rep(1,nFactors);
      if (!main.effects) output$var_v = (apply(v^2,2,sum)+v.sumvar)/nItems;
      #if(!main.effects) output$var_v = (sum(v^2) + v.sumvar) / (nItems * nFactors);
    }
    if( beta.int == T && fit.regression == F)
      {
        #z = cbind(rep(1,length(beta)),rep(0,length(beta)))
        #fit = fit.forMainEffect.bayesglm(beta, z, lm=lm,...);
        output$d0 = rep(0,dim(z)[2]); output$d0[1] = mean(beta);
        output$var_alpha = (sum(alpha^2) + alpha.sumvar) / nUsers;
        #output$var_beta = 1;
        output$var_beta = (sum((beta - mean(beta))^2) + beta.sumvar) / nItems;
        #if(!main.effects) output$var_u = (sum(u^2) + u.sumvar) / (nUsers * nFactors);
        #output$var_u = 1;
        #if(!main.effects) output$var_v = (sum(v^2) + v.sumvar) / (nItems * nFactors);
	#if (!main.effects) output$var_u = (apply(u^2,2,sum)+u.sumvar)/nItems;
        #output$var_v = rep(1,nFactors);
	output$var_u = rep(1,nFactors);
	if (!main.effects) output$var_v = (apply(v^2,2,sum)+v.sumvar)/nItems;
      }

    return(output);
}
