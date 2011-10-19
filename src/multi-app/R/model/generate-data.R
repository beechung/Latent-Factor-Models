### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### Generate data for testing/debugging/simulation
### NOTE: It's possible that some user IDs have no data in all applications
###
generate.GaussianData <- function(
	nUsers, nApps,
	nFeatures,      # nFeatures[k]: #features in application k
	nItems,         # nItems[k]:    #items    in application k
	nGlobalFactors, # number of global factors per user
	nLocalFactors,  # nLocalFactors[k]: # factors per user in application k
	frac.missing,   # frac.missing[k]:  fraction of missing observations in application k
	identity.A=FALSE, identity.B=FALSE, # whether A or B is an identity matrix
	is.y.logistic=FALSE, # whether to use a logistic model for y
	is.x.logistic=FALSE, # whether to use a logistic model for x
	x.bias=0, y.bias=0,
	A.sd=1, B.sd=1, b.sd=1, alpha.sd=1, beta.sd=1,
	var_x=0.05, var_y=0.05, var_z=1, # can be vectors of size nApps
	var_u=1,
	w_x.mean=1, w_x.var=0, w_y.mean=1, w_y.var=0
){
	if(length(nFeatures) == 1) nFeatures = rep(nFeatures, nApps);
	if(length(nItems)    == 1) nItems    = rep(nItems,    nApps);
	if(length(nLocalFactors) == 1) nLocalFactors = rep(nLocalFactors, nApps);
	if(length(frac.missing)  == 1) frac.missing  = rep(frac.missing,  nApps);
	if(length(var_x) == 1) var_x = rep(var_x, nApps);
	if(length(var_y) == 1) var_y = rep(var_y, nApps);
	if(length(var_z) == 1) var_z = rep(var_z, nApps);
	
	if(identity.A && any(nLocalFactors != nGlobalFactors)) stop("nLocalFactors != nGlobalFactors");
	if(identity.B && any(nLocalFactors != nFeatures)) stop("nLocalFactors != nFeatures");
	
	# generate parameters
	A = list(); B = list(); b = list(); alpha = list(); beta = list();
	for(k in 1:nApps){
		if(identity.A) A[[k]] = 1
		else           A[[k]] = matrix(rnorm(nLocalFactors[k]*nGlobalFactors, mean=0, sd=A.sd), nrow=nLocalFactors[k]);
		if(identity.B) B[[k]] = 1
		else           B[[k]] = matrix(rnorm(nFeatures[k]*nLocalFactors[k], mean=0, sd=B.sd), nrow=nFeatures[k]);
		b[[k]] = x.bias + rnorm(nFeatures[k], mean=0, sd=b.sd);
		alpha[[k]] = y.bias + rnorm(nItems[k], mean=0, sd=alpha.sd);
		beta[[k]]  = matrix(rnorm(nItems[k]*nLocalFactors[k], mean=0, sd=beta.sd), nrow=nItems[k]);
	}
	param = list(A=A, B=B, b=b, alpha=alpha, beta=beta, var_x=var_x, var_y=var_y, var_z=var_z, var_u=var_u)
	
	# generate observations
	feature  = NULL;
	response = NULL;
	for(k in 1:nApps){
		# generate features
		nPickedUsers = ceiling(nUsers * (1-frac.missing[k]));
		f = data.frame(
			user = rep(sample.int(n=nUsers, size=nPickedUsers, replace=FALSE), each=nFeatures[k]),
			app  = as.integer(k),
			index= rep(as.integer(1:nFeatures[k]), times=nPickedUsers),
			x    = NA
		);
		if(w_x.var > 0) f$w = rgamma.mean.var(n=nrow(f), mean=w_x.mean, var=w_x.var);
		feature = rbind(feature, f);
		# generate response
		r = data.frame(
			user = sample.int(n=nUsers, size=nPickedUsers*nItems[k], replace=TRUE),
			app  = as.integer(k),
			item = rep(as.integer(1:nItems[k]), each=nPickedUsers),
			y    = NA
		);
		if(w_y.var > 0) r$w = rgamma.mean.var(n=nrow(r), mean=w_y.mean, var=w_y.var);
		response = rbind(response, r);
	}
	
	output = generate.GaussianObs(
		feature=feature, response=response, param=param,
		is.x.logistic=is.x.logistic, is.y.logistic=is.y.logistic
	);
	return(output);
}

### Input:
###   feature  = data.frame(user, app, index, x, w)
###   response = data.frame(user, app, item,  y, w)
###   param    = list(A, B, b, alpha, beta, var_x, var_y, var_z, var_u)
### This function fill in feature$x and response$y based on param and 
### return list(feature, response, param, factor)
generate.GaussianObs <- function(
	feature, response, param, is.x.logistic, is.y.logistic
){
	nApps = length(param$A);
	if(length(param$B) != nApps) stop("length(param$B) != nApps");
	if(length(param$b) != nApps) stop("length(param$b) != nApps");
	if(length(param$alpha) != nApps) stop("length(param$alpha) != nApps");
	if(length(param$beta ) != nApps) stop("length(param$beta)  != nApps");
	if(length(param$var_x) != nApps) stop("length(param$var_x) != nApps");
	if(length(param$var_y) != nApps) stop("length(param$var_y) != nApps");
	if(length(param$var_z) != nApps) stop("length(param$var_z) != nApps");
	
	# generate global factors u
	nUsers = max(feature$user, response$user);
	nGlobalFactors = ncol(param$A[[1]]);
	u = matrix(rnorm(nUsers*nGlobalFactors, mean=0, sd=sqrt(param$var_u)), nrow=nUsers);
	# generate local factors z
	z = list();
	for(k in 1:nApps){
		z[[k]] = if(length(param$A[[k]])==1) u  *    param$A[[k]] 
	 	         else                        u %*% t(param$A[[k]]);
		z[[k]] = z[[k]] + rnorm(length(z[[k]]), mean=0, sd=sqrt(param$var_z[k]));
	}
	# generate observations
	if(is.x.logistic){
		temp = predict.x.from.z(feature=feature, param=param, z=z, add.noise=FALSE);
		x.prob = 1/(1+exp(-temp));
		feature$x = rbinom(n=nrow(feature),size=1,prob=x.prob);
		feature$w = NULL;
	}else{
		feature$x = predict.x.from.z(feature=feature, param=param, z=z, add.noise=TRUE);
	}
	if(is.y.logistic){
		temp = predict.y.from.z(response=response, param=param, z=z, add.noise=FALSE);
		y.prob = 1/(1+exp(-temp));
		response$y = rbinom(n=nrow(response),size=1,prob=y.prob);
		response$w = NULL;
	}else{
		response$y = predict.y.from.z(response=response, param=param, z=z, add.noise=TRUE);
	}
	
	output = list(feature=feature, response=response, param=param,
				  factor=list(u=u, z=z));

	if(is.x.logistic) output$x.prob = x.prob;
	if(is.y.logistic) output$y.prob = y.prob;
	
	return(output);
}
