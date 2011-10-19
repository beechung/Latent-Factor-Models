### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### Hierarchical smoothing: 1-dimensional, 2 levels
###
###		data = data.frame(item, category, obs, var);
###
### 		(item[n], category[n], obs[n], var[n])
###               i            k
###      This tuple specifies that we have observation obs[n] for item i in category k
###      with observation variance var[n].
###      An item can be in multiple categories.
### 
###	IMPORTANT NOTE: Item and category indices start from 1 (NOT 0)
###
### Model:
### 	 obs[n] ~ N(mean=b[i,k],    var=var_obs[n])
###      b[i,k] ~ N(mean=q[k]*a[i], var=var_b[k])
###      a[i]   ~ N(mean=0,         var=var_a)
###
### Output:
###		a$mean[i] =   E[a[i] | obs]
###		a$var[i]  = Var[a[i] | obs]
###	  b$mean[i,k] =   E[b[i,k] | obs]
###	  b$var[i,k]  = Var[b[i,k] | obs]
###   b$cov[i,k]  = Cov[b[i,k], a[i] | obs]
###
hierarchical_smoothing_1D_2Levels <- function(
		data, q, var_b, var_a, nItems, verbose=0, debug=0
){
	nCategories = as.integer(length(q));
	nObs = as.integer(nrow(data));
	nItems = as.integer(nItems);
	if(max(data$item) > nItems) stop("max(data$item) > nItems");
	if(max(data$category) > nCategories) stop("max(data$category) > nCategories");
	a = list(
			mean=rep(double(1), nItems),
			var =rep(double(1), nItems)
	);
	b = list(
			mean=matrix(rep(double(1), nItems*nCategories), nrow=nItems, ncol=nCategories),
			var =matrix(rep(double(1), nItems*nCategories), nrow=nItems, ncol=nCategories),
			cov =matrix(rep(double(1), nItems*nCategories), nrow=nItems, ncol=nCategories)
	);
	check_type_size(data$item, "integer", nObs);
	check_type_size(data$category, "integer", nObs);
	check_type_size(data$obs, "double", nObs);
	check_type_size(data$var, "double", nObs);
	check_type_size(q, "double", nCategories);
	check_type_size(var_b, "double", nCategories);
	check_type_size(var_a, "double", 1);
	
	.C("hierarchical_smoothing_1D_2Levels",
			a$mean, a$var,        # nItems x 1
			b$mean, b$var, b$cov, # nItems x nCategories
			# Input
			data$item, data$category, data$obs, data$var, # nObs x 1
			q, var_b, # nCategories x 1
			var_a, nObs, nItems, nCategories,
			# Control
			as.integer(verbose), as.integer(debug),
			DUP=FALSE
	);
	out = list(a=a, b=b);
	return(out);
}


###
### EM hierarchical smoothing: 1-dimensional, 2 levels
###
###		data = data.frame(item, category, obs, var);
###
### 		(item[n], category[n], obs[n], var[n])
###               i            k
###      This tuple specifies that we have observation obs[n] for item i in category k
###      with observation variance var[n].
###      An item can be in multiple categories.
### 
###	Observations can be pre-aggregated. The following two are equivalent:
###		          obs[n] ~ Normal(mu, sigma^2), for n=1,...,N
###     (sum_n obs[n])/N ~ Normal(mu, sigma^2/N)
###
###	IMPORTANT NOTE: Item and category indices start from 1 (NOT 0)
###
### Model:
### 	 obs[n] ~ N(mean = b[i,k] + sum(w_obs * x_obs[n,]) + offset[n],  var=var[n]*var_obs_adj)
###      b[i,k] ~ N(mean = q[k]*a[i] + sum(w_b[k,] * x_b[i,]),           var=var_b[k])
###      a[i]   ~ N(mean = 0,                                            var=var_a)
###
###   Features:  feature = list(x_obs, x_b);
###    Factors:  factor  = list(a=list(mean,var), b=list(mean,var,cov));
###	Parameters:  param   = list(var_obs_adj, var_b, var_a, q, w_obs, w_b);
###
###	if fix.param == TRUE, only one E-step will take place and no M-step.
###
fit.hierModel.1Dim2Levels <- function(
	data, lambda, nEM.iter, feature=NULL, factor=NULL, param=NULL, q.fixed=NULL, offset=NULL,
	init.var=1, init.q=1, init.factor.a.sd=1, init.factor.b.sd=1, fix.param=FALSE, var_a.fixed=NULL,
	data.test=NULL, feature.test=NULL, offset.test=NULL, smooth.test=FALSE, ignore.var.w_b=FALSE, verbose=1, debug=0
){
	if(!is.integer(data$item))     data$item = as.integer(data$item);
	if(!is.integer(data$category)) data$category = as.integer(data$category);
	if(is.null(param)) init.w_obs = TRUE else init.w_obs = FALSE;
	x_obs = feature$x_obs;
	x_b   = feature$x_b;
	nObs = nrow(data);
	nItems = if(!is.null(x_b)) nrow(x_b) else max(data$item);
	nCategories = max(data$category);
	
	if(smooth.test) warning("You probably should NOT smooth the test data!!");
	
	if(smooth.test && !is.null(data.test)){
		test.items = unique(data.test$item);
		select = data$item %in% test.items;
		data.test = rbind(data.test, data[select,]);
		if(!is.null(feature.test)) feature.test$x_obs = rbind(feature.test$x_obs, feature$x_obs[select,,drop=FALSE]);
		if(!is.null(offset.test))  offset.test = c(offset.test, offset[select]);
	}
	
	###
	###	Initialization
	###
	if(is.null(factor) && !fix.param){
		factor = list(
			a=list(mean=rnorm(nItems,mean=0,sd=init.factor.a.sd),
				    var=rep(1,nItems)),
		    b=list(mean=matrix(rnorm(nItems*nCategories,mean=0,sd=init.factor.b.sd), nrow=nItems),
				    var=matrix(1, nrow=nItems, ncol=nCategories),
					cov=matrix(1, nrow=nItems, ncol=nCategories))
		 );
	}
	if(is.null(param)){
		param = list(var_obs_adj=init.var, q=rep(init.q, nCategories), var_b=rep(init.var, nCategories), var_a=init.var)
		if(!is.null(x_obs)) param$w_obs = rep(0, ncol(x_obs));
		if(!is.null(x_b))   param$w_b   = matrix(0, nrow=nCategories, ncol=ncol(x_b));
	}
	if(is.null(data$var)) data$var = 1.0;
	
	if(!is.null(q.fixed)){
		if(length(q.fixed) == 1) q.fixed=rep(q.fixed, nCategories);
		param$q = q.fixed;
	}
	if(!is.null(var_a.fixed)) param$var_a = var_a.fixed;
	
	###
	### Sanity check
	###
	if(!fix.param){
		check_size(factor$a$mean, size=nItems); check_size(factor$a$var, size=nItems);
		check_size(factor$b$mean, size=c(nItems, nCategories));
		check_size(factor$b$var,  size=c(nItems, nCategories));
		check_size(factor$b$cov,  size=c(nItems, nCategories));
	}
	check_size(param$var_obs_adj, 1); check_size(param$var_a, 1);
	check_size(param$q, nCategories); check_size(param$var_b, nCategories);
	if(any(data$item > nItems)) stop("data$item > nItems");
	if(!is.null(x_obs)) check_size(param$w_obs, ncol(x_obs));
	if(!is.null(x_b))   check_size(param$w_b,   c(nCategories, ncol(x_b)));
	if(!is.null(x_b)   && nrow(x_b) != nItems) stop("nrow(x_b) != nItems");
	if(!is.null(x_obs) && nrow(x_obs) != nObs) stop("nrow(x_obs) != nObs");
	if(!is.null(data.test)){
		if(is.null(feature) != is.null(feature.test)) stop("is.null(feature) != is.null(feature.test)");
		if(is.null(offset)  != is.null(offset.test)) stop("is.null(offset) != is.null(offset.test)");
	}
	
	###
	### Fit w_obs first
	###
	if(!fix.param && !is.null(x_obs) && init.w_obs){
		y = data$obs;
		if(!is.null(offset)) y = y - offset;
		model = lm(y ~ x_obs - 1, weights=1/data$var, model=FALSE, qr=FALSE);
		param$w_obs = model$coefficients;
		param$w_obs[is.na(param$w_obs)] = 0;
	}
	if(verbose > 0 && !fix.param){
		Eloglik = Eloglik.hierModel.1Dim2Levels(data=data, a=factor$a, b=factor$b, param=param, lambda=lambda, x_obs=x_obs, x_b=x_b, offset=offset);
		cat("Start fitting Hierarchical Model: 1-dimensional, 2 levels\n",sep="");
		if(is.null(feature)) cat("    ** NO FEATURE IS SPECIFIED **\n");
		if(is.null(offset))  cat("    ** NO OFFSET  IS SPECIFIED **\n");
		cat("    Initial E[loglik] = ",Eloglik," + const\n",sep="");
		if(!is.null(data.test)){
			if(smooth.test) factor.test = smooth.hierModel.1Dim2Levels(data=data.test, param=param, feature=feature.test, offset=offset.test)
			else            factor.test = factor;
			rmse.weighted = predict.hierModel.1Dim2Levels(
					data=data.test, factor=factor.test, param=param, feature=feature.test, offset=offset.test
			)$rmse.weighted;
			cat("              RMSE.wt = ",rmse.weighted,"\n",sep="");
		}
	}
	
	data.EStep = data;
	ik = cbind(data$item, data$category);
	
	###
	### EM iterations
	###
	for(iter in 1:nEM.iter){
		if(verbose > 0) cat("  EM iteration: ",iter,"\n",sep="");
		
		data.EStep$obs = data$obs;
		if(!is.null(offset)) data.EStep$obs = data.EStep$obs - offset;
		if(!is.null(x_obs))  data.EStep$obs = data.EStep$obs - drop(x_obs %*% param$w_obs);
		if(!is.null(x_b)){
			gx = matrix(NA, nrow=nItems, ncol=nCategories);
			for(k in 1:nCategories) gx[,k] = drop(x_b %*% param$w_b[k,]);
			data.EStep$obs = data.EStep$obs - gx[ik];
		}
		data.EStep$var = data$var * param$var_obs_adj;
		
		factor = hierarchical_smoothing_1D_2Levels(
			data=data.EStep, q=param$q, param$var_b, param$var_a, nItems=nItems, verbose=verbose, debug=debug
		);

		if(!is.null(x_b)){
			# factor$b$mean needs to be corrected
			for(k in 1:nCategories){
				factor$b$mean[,k] = factor$b$mean[,k] + drop(x_b %*% param$w_b[k,]);
			}
		}
		
		if(verbose > 0){
			Eloglik = Eloglik.hierModel.1Dim2Levels(data=data, a=factor$a, b=factor$b, param=param, lambda=lambda, x_obs=x_obs, x_b=x_b, offset=offset);
			cat("    E-Step     end: E[loglik] = ",Eloglik," + const\n",sep="");
			if(!is.null(data.test)){
				if(smooth.test) factor.test = smooth.hierModel.1Dim2Levels(data=data.test, param=param, feature=feature.test, offset=offset.test)
				else            factor.test = factor;
				rmse.weighted = predict.hierModel.1Dim2Levels(
						data=data.test, factor=factor.test, param=param, feature=feature.test, offset=offset.test
				)$rmse.weighted;
				cat("                      RMSE.wt = ",rmse.weighted,"\n",sep="");
			}
		}
		
		if(fix.param) return(list(factor=factor, param=param));
		
		param = MStep.hierModel.1Dim2Levels(
			data=data, factor=factor, lambda=lambda, q.fixed=q.fixed, x_obs=x_obs, x_b=x_b, offset=offset,
			data.test=data.test, feature.test=feature.test, offset.test=offset.test, smooth.test=smooth.test, 
			fit.var_a=is.null(var_a.fixed), ignore.var.w_b=ignore.var.w_b,
			param.prev=param, verbose=verbose
		);

		if(verbose > 0){
			if(!is.null(data.test)){
				if(smooth.test) factor.test = smooth.hierModel.1Dim2Levels(data=data.test, param=param, feature=feature.test, offset=offset.test)
				else            factor.test = factor;
				rmse.weighted = predict.hierModel.1Dim2Levels(
						data=data.test, factor=factor.test, param=param, feature=feature.test, offset=offset.test
				)$rmse.weighted;
				cat("                      RMSE.wt = ",rmse.weighted,"\n",sep="");
				cat("        (M-step needs another E-step to reflect the benefit)\n");
			}
		}
	}
	
	out = list(factor=factor, param=param);
	return(out);
}

get.EStep.data <- function(data, param, feature, offset){
	ik = cbind(data$item, data$category);
	nItems = if(!is.null(feature$x_b)) nrow(feature$x_b) else max(data$item);
	nCategories = max(data$category);
	data.EStep = data;
	data.EStep$obs = data$obs;
	if(!is.null(offset)) data.EStep$obs = data.EStep$obs - offset;
	if(!is.null(feature$x_obs))  data.EStep$obs = data.EStep$obs - drop(feature$x_obs %*% param$w_obs);
	if(!is.null(feature$x_b)){
		gx = matrix(NA, nrow=nItems, ncol=nCategories);
		for(k in 1:nCategories) gx[,k] = drop(feature$x_b %*% param$w_b[k,]);
		data.EStep$obs = data.EStep$obs - gx[ik];
	}
	data.EStep$var = data$var * param$var_obs_adj;
	return(data.EStep);
}

correct.factor <- function(factor, param, feature, x_b){
	if(!is.null(x_b)){
		# factor$b$mean needs to be corrected
		for(k in 1:ncol(factor$b$mean)){
			factor$b$mean[,k] = factor$b$mean[,k] + drop(x_b %*% param$w_b[k,]);
		}
	}
	return(factor);
}

###
### M-Step for hierarchical smoothing: 1-dimensional, 2 levels
###
MStep.hierModel.1Dim2Levels <- function(
	data, factor, lambda, q.fixed=NULL, x_obs=NULL, x_b=NULL, offset=NULL,
	data.test=NULL, feature.test=NULL, offset.test=NULL, smooth.test=FALSE, fit.var_a=TRUE,
	ignore.var.w_b=FALSE, 
	param.prev=NULL, verbose=0
){
	a = factor$a;  b = factor$b;
	nObs = nrow(data);
	nItems = dim(b$mean)[1];
	nCategories = dim(b$mean)[2];
	if(max(data$item) > nItems) stop("max(data$item) > nItems");
	if(max(data$category) > nCategories) stop("max(data$category) > nCategories");
	check_size(a$mean, size=nItems);
	check_size(a$var,  size=nItems);
	check_size(b$var,  size=c(nItems, nCategories));
	check_size(b$cov,  size=c(nItems, nCategories));
	if(!is.null(x_b) && nrow(x_b) != nItems) stop("nrow(x_b) != nItems");
	ik = cbind(data$item, data$category);
	
	if(is.null(param.prev)){
		out = list(var_obs_adj=NA, q=rep(NA, nCategories), var_b=rep(NA, nCategories), var_a=NA)
		if(!is.null(x_obs)) out$w_obs = rep(NA, ncol(x_obs));
		if(!is.null(x_b))   out$w_b   = matrix(NA, nrow=nCategories, ncol=ncol(x_b));
	}else{
		out = param.prev; check_size(out$q, nCategories); check_size(out$var_b, nCategories);
		if(!is.null(x_obs)) check_size(out$w_obs, ncol(x_obs));
		if(!is.null(x_b))   check_size(out$w_b, c(nCategories, ncol(x_b)));
	}
	
	if(verbose > 0){
		if(is.null(param.prev)) stop("is.null(param.prev)");
		Eloglik = Eloglik.hierModel.1Dim2Levels(data=data, a=a, b=b, param=out, lambda=lambda, x_obs=x_obs, x_b=x_b, offset=offset);
		cat("    M-Step initial: E[loglik] = ",Eloglik," + const\n",sep="");
	}
	
	# Fit w_obs
	y = data$obs - b$mean[ik];
	if(!is.null(offset)) y = y - offset;
	if(!is.null(x_obs)){
		if(nrow(x_obs) != nObs) stop("nrow(x_obs) != nObs");
		model = lm(y ~ x_obs - 1, weights=1/data$var);
		cat("w_obs\n");
		print(summary(model));
		out$w_obs = model$coefficients;
		out$w_obs[is.na(out$w_obs)] = 0;
		err = y - drop(x_obs %*% out$w_obs);
	}else{
		err = y;
	}
	if(verbose > 1){
		Eloglik = Eloglik.hierModel.1Dim2Levels(data=data, a=a, b=b, param=out, lambda=lambda, x_obs=x_obs, x_b=x_b, offset=offset);
		cat("        before var: E[loglik] = ",Eloglik," + const\n",sep="");
	}
	
	# Fit var_obs_adj
	out$var_obs_adj = sum((err^2 + b$var[ik]) / data$var) / nObs;
	
	if(verbose > 0){
		Eloglik = Eloglik.hierModel.1Dim2Levels(data=data, a=a, b=b, param=out, lambda=lambda, x_obs=x_obs, x_b=x_b, offset=offset);
		cat("        after  obs: E[loglik] = ",Eloglik," + const\n",sep="");
		if(verbose > 1 && !is.null(data.test)){
			data.EStep = get.EStep.data(data=data, param=out, feature=list(x_obs=x_obs,x_b=x_b), offset=offset);
			factor.test = hierarchical_smoothing_1D_2Levels(
					data=data.EStep, q=out$q, out$var_b, out$var_a, nItems=nItems, verbose=0, debug=0
			);
			factor.test = correct.factor(factor=factor.test, param=out, x_b=x_b);
			#if(smooth.test) factor.test = smooth.hierModel.1Dim2Levels(data=data.test, param=out, feature=feature.test, offset=offset.test)
			#else            factor.test = factor;
			rmse.weighted = predict.hierModel.1Dim2Levels(
					data=data.test, factor=factor.test, param=out, feature=feature.test, offset=offset.test
			)$rmse.weighted;
			cat("                      RMSE.wt = ",rmse.weighted,"\n",sep="");
		}
	}
	
	# Fit w_b
	item_id_set = split(data$item, data$category);
	if(any(names(item_id_set) != c(1:nCategories))) stop("names(item_id_set) != c(1:nCategories)");
	for(k in 1:nCategories){
		items_k = unique(item_id_set[[k]]);
		str(items_k);
		if(length(items_k) == 0){ warning("No observation for category ",k); next; }
		
		a.mean   = a$mean[items_k];
		b_k.mean = b$mean[items_k,k];
		V_k      = sum(b$var[items_k,k]);
		W_k      = sum(b$cov[items_k,k]);
		V        = sum(a$var[items_k]);
		
		if(!is.null(q.fixed)) out$q[k] = q.fixed[k];
		if(!is.null(q.fixed) && is.null(x_b)){
			gx = 0;
		}else{
			x_k = x_b[items_k,,drop=FALSE];
			if(is.null(q.fixed)){
				X = cbind(x_k, a.mean);
				prior.mean = c(rep(0, length(x_b[1,])),      (lambda+W_k)/(lambda+V));
				prior.prec = c(rep(lambda, length(x_b[1,])),              (lambda+V));
			}else{
				X = x_k;
				prior.mean = rep(0, ncol(x_b));
				prior.prec = rep(lambda, ncol(x_b));
			}
			Sigma.inv = diag(prior.prec, ncol=length(prior.prec));
			
			if(ignore.var.w_b){
				# for testing purposes (this is not correct)
				beta = lm(b_k.mean ~ X - 1, model=FALSE, qr=FALSE)$coefficients;
				beta[is.na(beta)] = 0;
			}else{
				# this is the right formula
				beta = drop(solve(t(X) %*% X + Sigma.inv) %*% (t(X) %*% b_k.mean + Sigma.inv %*% prior.mean));
			}
			
			if(is.null(q.fixed)){
				out$q[k] = beta[length(beta)];
				if(!is.null(x_b)) out$w_b[k,] = beta[1:(length(beta)-1)];
			}else                 out$w_b[k,] = beta;
			
			if(!is.null(x_b)) gx = drop(x_k %*% out$w_b[k,])
			else              gx = 0;
		}
		q_k = out$q[k];
		out$var_b[k] = (sum( (b_k.mean - gx - q_k*a.mean)^2 ) + V_k - 2*W_k*q_k + V*(q_k^2)) / length(items_k);
		
		if(verbose > 2){
			Eloglik = Eloglik.hierModel.1Dim2Levels(data=data, a=a, b=b, param=out, lambda=lambda, x_obs=x_obs, x_b=x_b, offset=offset);
			cat("        after b[",k,"]: E[loglik] = ",Eloglik," + const\n",sep="");
		}
	}
	
	if(verbose > 0){
		Eloglik = Eloglik.hierModel.1Dim2Levels(data=data, a=a, b=b, param=out, lambda=lambda, x_obs=x_obs, x_b=x_b, offset=offset);
		cat("        after    b: E[loglik] = ",Eloglik," + const\n",sep="");
		if(verbose > 1 && !is.null(data.test)){
			data.EStep = get.EStep.data(data=data, param=out, feature=list(x_obs=x_obs,x_b=x_b), offset=offset);
			factor.test = hierarchical_smoothing_1D_2Levels(
					data=data.EStep, q=out$q, out$var_b, out$var_a, nItems=nItems, verbose=0, debug=0
			);
			factor.test = correct.factor(factor=factor.test, param=out, x_b=x_b);
			#if(smooth.test) factor.test = smooth.hierModel.1Dim2Levels(data=data.test, param=out, feature=feature.test, offset=offset.test)
			#else            factor.test = factor;
			rmse.weighted = predict.hierModel.1Dim2Levels(
					data=data.test, factor=factor.test, param=out, feature=feature.test, offset=offset.test
			)$rmse.weighted;
			cat("                      RMSE.wt = ",rmse.weighted,"\n",sep="");
		}
	}
	
	# Fit var_a
	if(fit.var_a){
		items = unique(data$item);
		out$var_a = sum((a$mean[items])^2 + a$var[items]) / length(items);
		if(verbose > 0){
			Eloglik = Eloglik.hierModel.1Dim2Levels(data=data, a=a, b=b, param=out, lambda=lambda, x_obs=x_obs, x_b=x_b, offset=offset);
			cat("        after    a: E[loglik] = ",Eloglik," + const\n",sep="");
			cat("                 var_a = ",out$var_a,"\n",sep="");
			if(verbose > 1 && !is.null(data.test)){
				data.EStep = get.EStep.data(data=data, param=out, feature=list(x_obs=x_obs,x_b=x_b), offset=offset);
				factor.test = hierarchical_smoothing_1D_2Levels(
						data=data.EStep, q=out$q, out$var_b, out$var_a, nItems=nItems, verbose=0, debug=0
				);
				factor.test = correct.factor(factor=factor.test, param=out, x_b=x_b);
				#if(smooth.test) factor.test = smooth.hierModel.1Dim2Levels(data=data.test, param=out, feature=feature.test, offset=offset.test)
				#else            factor.test = factor;
				rmse.weighted = predict.hierModel.1Dim2Levels(
						data=data.test, factor=factor.test, param=out, feature=feature.test, offset=offset.test
				)$rmse.weighted;
				cat("                      RMSE.wt = ",rmse.weighted,"\n",sep="");
			}
		}
	}
	return(out);
}
###
### E[log(likelihood)]
###
Eloglik.hierModel.1Dim2Levels <- function(
	data, 
	a, b, param, lambda=0, x_obs=NULL, x_b=NULL, offset=NULL
){
	w_obs=param$w_obs; var_obs_adj=param$var_obs_adj; q=param$q; w_b=param$w_b; var_b=param$var_b; var_a=param$var_a;
	nObs = nrow(data);
	nItems = dim(b$mean)[1];
	nCategories = dim(b$mean)[2];
	if(max(data$item) > nItems) stop("max(data$item) > nItems");
	if(max(data$category) > nCategories) stop("max(data$category) > nCategories");
	check_size(a$mean, size=nItems);
	check_size(a$var,  size=nItems);
	check_size(b$var,  size=c(nItems, nCategories));
	check_size(b$cov,  size=c(nItems, nCategories));
	
	ik = cbind(data$item, data$category);
	
	err = data$obs - b$mean[ik];
	if(!is.null(offset)) err = err - offset;
	if(!is.null(x_obs))  err = err - drop(x_obs %*% w_obs);
	loglik.obs = - (1/2) * sum( (err^2 + b$var[ik])/(data$var * var_obs_adj) + log(data$var * var_obs_adj));
	
	loglik = rep(0, nCategories);
	item_id_set = split(data$item, data$category);
	if(any(names(item_id_set) != c(1:nCategories))) stop("names(item_id_set) != c(1:nCategories)");
	for(k in 1:nCategories){
		items_k = unique(item_id_set[[k]]);
		a.mean   = a$mean[items_k];
		b_k.mean = b$mean[items_k,k];
		V_k      = sum(b$var[items_k,k]);
		W_k      = sum(b$cov[items_k,k]);
		V        = sum(a$var[items_k]);
		x_k      = x_b[items_k,,drop=FALSE];
		
		if(!is.null(x_b)){ gx = drop(x_k %*% w_b[k,]); g2 = sum(w_b[k,]^2);}
		else             { gx = 0;                     g2 = 0;}
		
		err = b_k.mean - gx - q[k]*a.mean;
		reg = lambda*((q[k]-1)^2 + g2);
		loglik[k] = - (1/2) * ( (sum(err^2) + V_k - 2*W_k*q[k] + V*(q[k]^2) + reg)/var_b[k]  +  length(items_k)*log(var_b[k]));
	}
	loglik.b = sum(loglik);
	
	items = unique(data$item);
	loglik.a = - (1/2) * ( sum(a$mean[items]^2 + a$var[items])/var_a  +  length(items)*log(var_a) );
	
	return(loglik.obs + loglik.b + loglik.a);
}
###
### Prediction
###
predict.hierModel.1Dim2Levels <- function(
	data, factor, param, feature=NULL, offset=NULL
){
	nObs = nrow(data);
	nItems = nrow(factor$b$mean);
	nCategories = ncol(factor$b$mean);
	x_obs = feature$x_obs;
	x_b   = feature$x_b;
	if(is.null(data$var)) data$var = 1;
	###
	### Sanity check
	###
	if(max(data$item) > nItems) stop("max(data$item) > nItems");
	if(max(data$category) > nCategories) stop("max(data$category) > nCategories");
	check_size(factor$b$mean, size=c(nItems, nCategories));
	if(!is.null(x_obs)) check_size(param$w_obs, ncol(x_obs));
	if(!is.null(x_obs) && nrow(x_obs) != nObs) stop("nrow(x_obs) != nObs");
	
	ik = cbind(data$item, data$category);
	
	pred = factor$b$mean[ik];
	if(!is.null(offset)) pred = pred + offset;
	if(!is.null(x_obs))  pred = pred + drop(x_obs %*% param$w_obs);
	
	rmse = sqrt(mean((pred - data$obs)^2));
	rmse.weighted = sqrt(sum((pred - data$obs)^2 / data$var) / sum(1/data$var));
	out = list(prediction=pred, rmse=rmse, rmse.weighted=rmse.weighted);
	return(out);
}
###
### smoothing with fitted parameters
###
###    output: list(a=list(mean,var), b=list(mean,var,cov));
###
smooth.hierModel.1Dim2Levels <- function(
		data, param, feature=NULL, offset=NULL, verbose=0, debug=0
){
	ans = fit.hierModel.1Dim2Levels(
			data=data, lambda=0, nEM.iter=1, feature=feature, factor=NULL, param=param, q.fixed=NULL, offset=offset,
			init.var=1, init.q=1, init.factor.a.sd=1, init.factor.b.sd=1, fix.param=TRUE,
			data.test=NULL, feature.test=NULL, offset.test=NULL, verbose=verbose, debug=debug
	);
	return(ans$factor);
}

###
### Baseline models
###
fit.a.only.1Dim2Levels <- function(
	data, offset=NULL
){
	y = data$obs;
	if(!is.null(offset)) y = y - offset;
	nObs = nrow(data);
	nItems = max(data$item);
	nCategories = max(data$category);
	global.mean = mean(y);
	agg = aggregate(y, list(item=data$item), mean);
	a = list(mean = rep(global.mean, nItems));
	a$mean[agg$item] = agg$x;
	b = list(mean = matrix(global.mean, nrow=nItems, ncol=nCategories));
	for(k in 1:nCategories){
		b$mean[,k] = a$mean;
	}
	return(list(a=a, b=b));
}
fit.b.only.1Dim2Levels <- function(
	data, offset=NULL
){
	y = data$obs;
	if(!is.null(offset)) y = y - offset;
	out = fit.a.only.1Dim2Levels(data=data, offset=offset);
	agg = aggregate(y, list(item=data$item, category=data$category), mean);
	indices = cbind(agg$item, agg$category);
	out$b$mean[indices] = agg$x;
	return(out);
}


###
### Hierarchical smoothing: multi-dimensional, 2 levels
###
###		Very slow ... for debug only :)
###
###		data = list(item, category, obs, var);
###
### 		(item[n], category[n], obs[n,], var[n])
###               i            k
###      This tuple specifies that we have observation obs[n,] for item i in category k
###      with observation variance var[n] * I.
###      obs[n,] can be multi-dimensional
###      An item can be in multiple categories.
### 
###	IMPORTANT NOTE: Item and category indices start from 1 (NOT 0)
###
### Model:
### 	 obs[n,] ~ N(mean = C[k,,] %*% b[i,k,],  var =   var[n] * I)
###      b[i,k,] ~ N(mean = A[k,,] %*% a[i,],    var = var_b[k] * I)
###      a[i,]   ~ N(mean = 0,                   var = var_a    * I)
###
### Output:
###		a$mean[i,]  =   E[a[i,] | obs]
###		a$var[i,,]  = Var[a[i,] | obs]
###	  b$mean[i,k,]  =   E[b[i,k,] | obs]
###	  b$var[i,k,,]  = Var[b[i,k,] | obs]
###   b$cov[i,k,,]  = Cov[b[i,k,], a[i,] | obs]
###
hierarchical_smoothing_2Levels.forDebug <- function(
	data, A, C=NULL, var_b, var_a, nItems, ignore.nodes.without.obs=TRUE, verbose=0, debug=0
){
	y = data$obs;
	if(is.vector(y)) y = array(y, dim=c(length(y), 1));
	if(is.vector(A)) A = array(A, dim=c(length(A), 1, 1));
	nCategories = dim(A)[1];
	dim.b = dim(A)[2];
	dim.a = dim(A)[3];
	dim.y = dim(y)[2];
	if(is.null(C)){
		if(dim.y != dim.b) stop("is.null(C) and dim.y != dim.b");
		C = array(0, dim=c(nCategories, dim.y, dim.b));
		for(k in 1:nCategories) C[k,,] = diag(dim.y);
	}
	a = list(
			mean=array(NA, dim=c(nItems, dim.a)),
			var =array(NA, dim=c(nItems, dim.a, dim.a))
	);
	b = list(
			mean=array(NA, dim=c(nItems, nCategories, dim.b)),
			var =array(NA, dim=c(nItems, nCategories, dim.b, dim.b)),
			cov =array(NA, dim=c(nItems, nCategories, dim.b, dim.a))
	);
	    b_ik_ik = array(NA, dim=c(nItems, nCategories, dim.b));
	Gamma_ik_ik = array(NA, dim=c(nItems, nCategories, dim.b, dim.b));
	     a_i_ik = array(NA, dim=c(nItems, nCategories, dim.a));
	 Gamma_i_ik = array(NA, dim=c(nItems, nCategories, dim.a, dim.a));
	  	 B      = array(NA, dim=c(nItems, nCategories, dim.a, dim.b));

	Sigma_pa = var_a * diag(dim.a);
	Sigma_pa.inv = (1/var_a) * diag(dim.a);
	for(i in 1:nItems){
		sum.for.a.mean = rep(0, dim.a);
		sum.for.a.var  = array(0, dim=c(dim.a, dim.a))
		for(k in 1:nCategories){
			Sigma_ik = A[k,,] %*% Sigma_pa %*% t(A[k,,]) + (var_b[k] * diag(dim.b));
			Sigma_ik.inv = solve(Sigma_ik);
			B[i,k,,] = Sigma_pa %*% t(A[k,,]) %*% Sigma_ik.inv;
			B_ik = B[i,k,,];
			R_ik = Sigma_pa - B_ik %*% A[k,,] %*% Sigma_pa;
			
			this = (data$item == i & data$category == k);
			
			y.this = y[this,];  var.this = data$var[this];
			sum.for.b.mean = rep(0, dim.b);
			sum.for.b.var  = array(0, dim=c(dim.b, dim.b));
			if(length(y.this) > 0){ for(n in 1:length(y.this)){
				sum.for.b.mean = sum.for.b.mean + (t(C[k,,]) %*% y.this[n]) / var.this[n];
				sum.for.b.var  = sum.for.b.var  + (t(C[k,,]) %*% C[k,,]) / var.this[n];
				# cat("   i=",i,",k=",k,":  sum.for.b.mean=",sum.for.b.mean,", sum.for.b.var=",sum.for.b.var,"\n",sep="");
			}}
			
			Gamma_ik_ik[i,k,,] = solve( Sigma_ik.inv + sum.for.b.var );
			b_ik_ik[i,k,] = Gamma_ik_ik[i,k,,] %*% sum.for.b.mean;
			
			a_i_ik[i,k,] = B_ik %*% b_ik_ik[i,k,];
			Gamma_i_ik[i,k,,] = B_ik %*% Gamma_ik_ik[i,k,,] %*% B_ik + R_ik;
			
			if(length(y.this) == 0 && ignore.nodes.without.obs) next;
			
			sum.for.a.mean = sum.for.a.mean + solve(Gamma_i_ik[i,k,,]) %*% a_i_ik[i,k,];
			sum.for.a.var  = sum.for.a.var  + solve(Gamma_i_ik[i,k,,]) - Sigma_pa.inv;
		}
		cat("   i=",i,":  sum.for.a.mean=",sum.for.a.mean,", sum.for.a.var=",sum.for.a.var,"\n",sep="");
		a$var[i,,] = solve( Sigma_pa.inv + sum.for.a.var );
		a$mean[i,] = a$var[i,,] %*% sum.for.a.mean;
		
		for(k in 1:nCategories){
			G_ik_ik = Gamma_ik_ik[i,k,,];
			G_i_ik  = Gamma_i_ik[i,k,,];
			G_i_ik.inv = solve(G_i_ik);
			B_ik = B[i,k,,];  B_ik.t = t(B_ik);
			b$mean[i,k,] = b_ik_ik[i,k,] + G_ik_ik %*% B_ik.t %*% G_i_ik.inv %*% (a$mean[i,] - a_i_ik[i,k,]);
			b$var[i,k,,] = G_ik_ik + G_ik_ik %*% B_ik.t %*% G_i_ik.inv %*% (a$var[i,,] - G_i_ik) %*% G_i_ik.inv %*% B_ik %*% G_ik_ik;
			b$cov[i,k,,] = G_ik_ik %*% B_ik.t %*% G_i_ik.inv %*% a$var[i,,];
		}
	}
	out = list(a=a, b=b);
	return(out);
}
