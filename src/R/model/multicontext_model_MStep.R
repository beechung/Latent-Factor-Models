### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### Set x = 0 to disable b
###
fit.b <- function(
	y, x, fScore.mean, fScore.var, lambda, weight=NULL, algo=NULL, control=NULL, max.nObs=NULL
){
	if(!any(x != 0)){
		out = list(
			b = rep(0, ncol(x)),
			var_y = sum( (y - fScore.mean)^2 + fScore.var ) / length(y)
		);
		return(out);
	}
	
	if(is.null(algo) || ncol(x) == 1){
		ans = fit.lm.random.effect(
			response.mean = y - fScore.mean, 
			 feature.mean = x,
			response.var  = fScore.var, 
			weight=weight, lambda=lambda
		);
	}else{
		ans = reg.train(
			x = x, y = y - fScore.mean, algo=algo, 
			y.var=fScore.var, weight=weight, control=control, max.nObs=max.nObs
		);
	}
	return(list(b=ans$coeff, var_y=ans$var));
}

###
### x[subset[[k]],] selects the nodes for context k
###
fit.mainEffectPrior <- function(
	x, local.mean, global.mean, local.var, global.var, local.cov,
	lambda, subset=NULL, zero.mean=FALSE, algo=NULL, control=NULL
){
	nContexts = if(is.vector(local.mean)) 1 else ncol(local.mean);
	nFeatures = ncol(x);
	if(is.null(zero.mean) || is.na(zero.mean)) zero.mean=FALSE;
	if(nContexts == 1){
		if(!is.null(subset)){
			if(is.null(subset[[1]])) stop("is.null(subset[[1]])");
			x = x[subset[[1]],,drop=FALSE];
			local.mean = local.mean[subset[[1]]];
			local.var  = local.var[subset[[1]]];
		}
		if(!any(x != 0) || zero.mean){
			return(list(
				coeff = matrix(0, nrow=ncol(x), ncol=1),
				  var = sum(local.mean^2 + local.var) / length(local.mean)
			));
		}
		if(is.null(algo) || ncol(x) == 1){
			ans = fit.lm.random.effect(
				response.mean = local.mean, 
				 feature.mean = x,
				response.var  = local.var, 
				  lambda      = lambda
			);
			return(list(coeff=matrix(ans$coeff,ncol=1), var=ans$var));
		}else{
			ans = reg.train(
				x = x, y = local.mean, algo=algo, 
				y.var=local.var, control=control
			);
			return(list(coeff=list(ans$coeff), var=ans$var));
		}
	}
	out = list(coeff=matrix(NA, nrow=nFeatures, ncol=nContexts), slope=rep(NA, nContexts), var=rep(NA,nContexts));
	for(k in 1:nContexts){
		if(!is.null(subset)){
			if(is.null(subset[[k]])) stop("is.null(subset[[",k,"]])");
			x_k = x[subset[[k]],,drop=FALSE];
			local.mean_k = local.mean[subset[[k]],k];
			local.var_k  = local.var[subset[[k]],k];
			local.cov_k  = local.cov[subset[[k]],k];
			global.mean_k = global.mean[subset[[k]]];
			global.var_k  = global.var[subset[[k]]];
		}else{
			x_k = x;
			local.mean_k = local.mean[,k];
			local.var_k  = local.var[,k];
			local.cov_k  = local.cov[,k];
			global.mean_k = global.mean;
			global.var_k  = global.var;
		}
		if(!any(x != 0) || zero.mean){
			out$coeff[,k] = 0;
			z = matrix(global.mean_k, ncol=1);
			Delta = matrix(global.var_k, ncol=1);
			c = matrix(local.cov_k, ncol=1);
		}else{
			z = cbind(global.mean_k, as.matrix(x_k));
			Delta = cbind(global.var_k, matrix(0, nrow=nrow(z), ncol=nFeatures));
			c = cbind(local.cov_k, matrix(0, nrow=nrow(z), ncol=nFeatures));
		}
		ans = fit.lm.random.effect(
			response.mean = local.mean_k, 
			response.var  = local.var_k,
			 feature.mean = z,
			 feature.var  = Delta, 
			 feature.cov  = c,
			  lambda      = lambda
		);
		if(!any(x != 0) || zero.mean){
			out$slope[k] = ans$coeff;
		}else{
			out$slope[k] = ans$coeff[1];
			out$coeff[,k] = ans$coeff[2:length(ans$coeff)];
		}
		out$var[k] = ans$var;
	}
	return(out);
}

### If there are context-specific factors,
### 	x[subset[[k]],] selects the nodes for context k.
### Otherwise,
###		x[subset,] selects the nodes
fit.interactionPrior <- function(
	x, factor.mean, factor.var, lambda, zero.mean=FALSE, subset=NULL, algo=NULL, control=NULL
){
	nFeatures = ncol(x);
	nFactors  = ncol(factor.mean);
	nNodes    = nrow(x);
	if(nrow(factor.mean) != nNodes) stop("nrow(x) != nrow(factor.mean)");
	if(ncol(factor.var) != nFactors) stop("ncol(factor.var) != nFactors");
	
	out = list();
	loss = rep(NA, nFactors);
	num  = rep(NA, nFactors);
	
	if(is.null(zero.mean) || is.na(zero.mean)) zero.mean=FALSE;
	
	if(is.null(algo) || !any(x != 0) || zero.mean || ncol(x) == 1) out$coeff = matrix(NA, nrow=nFeatures, ncol=nFactors)
	else                                                           out$coeff = list();
	
	for(f in 1:nFactors){
		if(is.null(subset)){
			f.mean = factor.mean[,f];
			f.var  = factor.var[ ,f];
			x.f    = x;
		}else if(is.list(subset)){
			nContexts = length(subset);
			k = get.context.index(f, nContexts, nFactors);
			if(is.null(subset[[k]])) stop("is.null(subset[[",k,"]])");
			f.mean = factor.mean[subset[[k]],f];
			f.var  = factor.var[ subset[[k]],f];
			x.f    = x[subset[[k]],,drop=FALSE];
		}else if(is.vector(subset)){
			f.mean = factor.mean[subset,f];
			f.var  = factor.var[ subset,f];
			x.f    = x[subset,,drop=FALSE];
		}else stop("Unknown type for select");
		
		if(!any(x.f != 0) || zero.mean){
			out$coeff[,f] = 0;
			loss[f] = sum(f.mean^2 + f.var);
			num[f]  = length(f.mean)
		}else{
			if(is.null(algo) || ncol(x) == 1){
				ans = fit.lm.random.effect(
					response.mean = f.mean, 
					response.var  = f.var, 
					feature.mean = x.f,
					lambda      = lambda
				);
				out$coeff[,f] = ans$coeff;
			}else{
				ans = reg.train(
					x = x.f, y = f.mean, algo=algo, 
					y.var=f.var, control=control
				);
				out$coeff[[f]] = ans$coeff;
			}
			loss[f] = ans$loss;  num[f] = ans$num;
		}
	}
	if(is.list(subset)){
		nContexts = length(subset);
		out$var = rep(NA, nFactors);
		for(k in 1:nContexts){
			f = select.factor.indices(k, nContexts, nFactors);
			out$var[f] = sum(loss[f]) / sum(num[f]);
		}
	}else{
		out$var = sum(loss) / sum(num);
	}
	return(out);
}


###
### Monte-Carlo EM (M-step)
###
###   zero.mean["alpha"] = TRUE  ->  g = 0.
###   zero.mean["v"]     = TRUE  ->  D = 0.
###   fix.var["FACTOR"]  = a     ->  var_FACTOR = a
###
MCEM_MStep.multicontext <- function(
    factor.mean, factor.var, factor.cov, obs, feature, param, subset.info,
	ridge.lambda=rep(0,0), zero.mean=rep(0,0), fix.var=list(), is.logistic, max.nObs.for.b=NULL,
	debug=0, verbose=0, ...
){
    size = syncheck.multicontext.spec(factor=factor.mean, obs=obs, feature=feature, param=param);
	
	if(!all(names(ridge.lambda) %in% c("b","g0","d0","h0","G","D","H")))  stop("names(ridge.lambda)");
	if(!all(names(zero.mean) %in% c("alpha","beta","gamma","u","v","w"))) stop("names(zero.mean)");
	if(!all(names(fix.var) %in% c("alpha","beta","gamma","u","v","w")))   stop("names(fix.var)");
	
    nObs      = size$nObs;
    nSrcNodes = size$nSrcNodes;
    nDstNodes = size$nDstNodes;
	nFactors  = size$nFactors;
	
    # determine b and var_y
    b.time = proc.time();
	if(param$approx.interaction){
		# predict E[uv] as E[u]E[v].
		fScore = predict.y.from.factors(obs, factor.mean, feature, param, ignore.xb=TRUE);
	}else{
		fScore = factor.mean$fScore;
	}
	weight = NULL;
	if(is.logistic) weight = 1/param$var_y; # for logistic
	ans = fit.b(
		y=obs$y, x=feature$x_obs, lambda=ridge.lambda["b"],
		fScore.mean=fScore, fScore.var=factor.var$fScore, weight=weight,
		algo=param$reg.algo, control=param$reg.control, max.nObs=max.nObs.for.b
	);
	param$b = ans$b;
	if(is.null(weight)) param$var_y = ans$var_y
	else                param$var_y = ans$var_y / weight;
    time.used = proc.time() - b.time;
	if(verbose >= 2 || debug > 0) loglik = get.logLikelihood(
		prefix="  after   obs", suffix=sprintf("  fitting time %.2f sec",time.used[3]), prev.loglik=NULL,
		obs=obs, factor=factor.mean, feature=feature, param=param, factor.var=factor.var, factor.cov=factor.cov, subset.info=subset.info, verbose=verbose-1, is.logistic=is.logistic
	);

    # determin g0 and var_alpha
    if(!is.null(param$var_alpha)){
        b.time = proc.time(); 
		
		ans = fit.mainEffectPrior(
				x=feature$x_src, lambda=ridge.lambda["g0"], zero.mean=zero.mean["alpha"],
				local.mean=factor.mean$alpha, global.mean=factor.mean$alpha_global, 
				local.var =factor.var$alpha,  global.var =factor.var$alpha_global,
				local.cov =factor.cov$alpha,  subset=subset.info$src.context,
				algo=param$reg.algo, control=param$reg.control
		);
		param$g0 = ans$coeff;
		if(!is.null(ans$slope)) param[["q"]] = ans$slope;
		if("alpha" %in% names(fix.var)) param$var_alpha[] = fix.var[["alpha"]]
        else                            param$var_alpha   = ans$var;
		
        time.used = proc.time() - b.time;
		if(verbose >= 2 || debug > 0) loglik = get.logLikelihood(
			prefix="  after alpha", suffix=sprintf("  fitting time %.2f sec",time.used[3]), prev.loglik=loglik,
			obs=obs, factor=factor.mean, feature=feature, param=param, factor.var=factor.var, factor.cov=factor.cov, subset.info=subset.info, verbose=verbose-1, is.logistic=is.logistic
		);
	}
    # determin d0 and var_beta
    if(!is.null(param$var_beta)){
		b.time = proc.time(); 
		
		ans = fit.mainEffectPrior(
				x=feature$x_dst, lambda=ridge.lambda["d0"], zero.mean=zero.mean["beta"],
				local.mean=factor.mean$beta, global.mean=factor.mean$beta_global, 
				local.var =factor.var$beta,  global.var =factor.var$beta_global,
				local.cov =factor.cov$beta,  subset=subset.info$dst.context,
				algo=param$reg.algo, control=param$reg.control
		);
		param$d0 = ans$coeff;
		if(!is.null(ans$slope)) param[["r"]] = ans$slope;
		if("beta" %in% names(fix.var))  param$var_beta[] = fix.var[["beta"]]
		else                            param$var_beta   = ans$var;
		
		time.used = proc.time() - b.time;
		if(verbose >= 2 || debug > 0) loglik = get.logLikelihood(
			prefix="  after  beta", suffix=sprintf("  fitting time %.2f sec",time.used[3]), prev.loglik=loglik,
			obs=obs, factor=factor.mean, feature=feature, param=param, factor.var=factor.var, factor.cov=factor.cov, subset.info=subset.info, verbose=verbose-1, is.logistic=is.logistic
		);
	}
	# determin h0 and var_gamma
	if(!is.null(param$var_gamma)){
		b.time = proc.time(); 
		
		ans = fit.mainEffectPrior(
				x=feature$x_ctx, lambda=ridge.lambda["h0"], zero.mean=zero.mean["gamma"],
				local.mean=factor.mean$gamma, global.mean=NULL, 
				local.var =factor.var$gamma,  global.var =NULL,
				local.cov =NULL,              subset=NULL
		);
		param$h0 = if(is.list(ans$coeff)) ans$coeff[[1]] else ans$coeff[,1];
		if("gamma" %in% names(fix.var))  param$var_gamma[] = fix.var[["gamma"]]
		else                             param$var_gamma   = ans$var;
		
		time.used = proc.time() - b.time;
		if(verbose >= 2 || debug > 0) loglik = get.logLikelihood(
			prefix="  after gamma", suffix=sprintf("  fitting time %.2f sec",time.used[3]), prev.loglik=loglik,
			obs=obs, factor=factor.mean, feature=feature, param=param, factor.var=factor.var, factor.cov=factor.cov, subset.info=subset.info, verbose=verbose-1, is.logistic=is.logistic
		);
	}
	# determin G and var_u
	if(nFactors > 0 && !is.null(param$var_u)){
		b.time = proc.time(); 
		if(is.null(param$nLocalFactors)) subset = subset.info$src.id
		else                             subset = subset.info$edge.context;
		
		ans = fit.interactionPrior(
				x=feature$x_src, lambda=ridge.lambda["G"], zero.mean=zero.mean["u"],
				factor.mean=factor.mean$u, factor.var=factor.var$u, subset=subset,
				algo=param$reg.algo, control=param$reg.control
		);
		param$G = ans$coeff;
		if("u" %in% names(fix.var)) param$var_u[] = fix.var[["u"]]
		else                        param$var_u   = ans$var;
		
		time.used = proc.time() - b.time;
		if(verbose >= 2 || debug > 0) loglik = get.logLikelihood(
			prefix="  after     u", suffix=sprintf("  fitting time %.2f sec",time.used[3]), prev.loglik=loglik,
			obs=obs, factor=factor.mean, feature=feature, param=param, factor.var=factor.var, factor.cov=factor.cov, subset.info=subset.info, verbose=verbose-1, is.logistic=is.logistic
		);
	}
	# determin D and var_v
	if(nFactors > 0 && !is.null(param$var_v)){
		b.time = proc.time(); 
		if(!is.null(factor.mean$u)){
			if(is.null(param$nLocalFactors)) subset = subset.info$dst.id
			else                             subset = subset.info$edge.context;
		}else{
			if(is.null(param$nLocalFactors)) subset = subset.info$any.id
			else                             subset = subset.info$edge.context;
		}
		
		ans = fit.interactionPrior(
				x=feature$x_dst, lambda=ridge.lambda["D"], zero.mean=zero.mean["v"],
				factor.mean=factor.mean$v, factor.var=factor.var$v, subset=subset,
				algo=param$reg.algo, control=param$reg.control
		);
		param$D = ans$coeff;
		if("v" %in% names(fix.var)) param$var_v[] = fix.var[["v"]]
		else                        param$var_v   = ans$var;
		
		time.used = proc.time() - b.time;
		if(verbose >= 2 || debug > 0) loglik = get.logLikelihood(
			prefix="  after     v", suffix=sprintf("  fitting time %.2f sec",time.used[3]), prev.loglik=loglik,
			obs=obs, factor=factor.mean, feature=feature, param=param, factor.var=factor.var, factor.cov=factor.cov, subset.info=subset.info, verbose=verbose-1, is.logistic=is.logistic
		);
	}
	# determin H and var_w
	if(nFactors > 0 && !is.null(param$var_w) && is.null(param$nLocalFactors)){
		b.time = proc.time(); 
		
		ans = fit.interactionPrior(
				x=feature$x_ctx, lambda=ridge.lambda["H"], zero.mean=zero.mean["w"],
				factor.mean=factor.mean$w, factor.var=factor.var$w, subset=NULL,
				algo=param$reg.algo, control=param$reg.control
		);
		param$H = ans$coeff;
		if("w" %in% names(fix.var)) param$var_w[] = fix.var[["w"]]
		else                        param$var_w   = ans$var;
		
		time.used = proc.time() - b.time;
		if(verbose >= 2 || debug > 0) loglik = get.logLikelihood(
			prefix="  after     w", suffix=sprintf("  fitting time %.2f sec",time.used[3]), prev.loglik=loglik,
			obs=obs, factor=factor.mean, feature=feature, param=param, factor.var=factor.var, factor.cov=factor.cov, subset.info=subset.info, verbose=verbose-1, is.logistic=is.logistic
		);
	}
	if(!is.null(param$nLocalFactors) && param$var_w != 0) stop("!is.null(param$nLocalFactors) && param$var_w != 0");
	
    return(param);
}
