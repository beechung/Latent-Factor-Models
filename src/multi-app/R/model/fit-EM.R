### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### Use the EM algorithm to fit the model
### NOTE:
###   * ridge.lambda[c("A", "B", "beta")] are the numbers to be added to the
###     diagonal when fitting the regression to get A, B, and beta
###
fit.EM <- function(
	# Input data
	feature,  # data.frame(user, app, index, x, w)
	response, # data.frame(user, app, item,  y, w)
	# Model setup
	nGlobalFactors,     # num of global factors per user
	nLocalFactors=NULL, # num of local  factors per user (may be a vector of length: #app)
	identity.A=FALSE, identity.B=FALSE, # whether A/B is a identity matrix
	is.y.logistic=FALSE, # whether to use a logistic model for y
	is.x.logistic=FALSE, # whether to use a logistic model for x
	# Test data (optional)
	test.feature=NULL,
	test.response=NULL,
	# Model-fitting parameters
	nIter=30, # num of EM iterations
	ridge.lambda=c(A=0, B=0, beta=0), # lambda values for ridge regression for A, B, beta
	keep.users.without.obs.in.E.step=FALSE,
	keep.users.without.obs.in.M.step=FALSE,
	fix.var_u=1, # Fix var_u to some number (set to NULL if you do not want to fix it)
	# Initialization parameters
	param=NULL, # directly set all the parameters (if not NULL, the following will not be used)
	var_x=1, var_y=1, var_z=1, # initial var (may be vectors of length: #app)
	var_u=1,                   # initial var (length: 1)
	A.sd=1, B.sd=1, beta.sd=1, # A_k ~ N(mean=0, sd=A.sd), and so on.
	# Output options
	out.level=0,  # out.level=1: Save the model in out.dir/model.last and out.dir/model.minTestLoss
	out.dir=NULL, # out.level=2: Save the model after each iteration i in out.dir/model.i
	out.overwrite=FALSE, # whether to allow overwriting existing files
	# Debugging options
	debug=0, verbose=0, use.C=TRUE,
	show.marginal.loglik=FALSE
){
	option = 0;
	if(keep.users.without.obs.in.E.step) option = option + 1;
	if(keep.users.without.obs.in.M.step) option = option + 2;
	option = as.integer(option);
	
	if(!is.null(test.response) && out.level <= 0 && verbose <= 0) stop("!is.null(test.response) && out.level <= 0 && verbose <= 0");
	if(out.level > 0){
		if(is.null(out.dir)) stop("out.dir = NULL");
		if(file.exists(paste(out.dir,"/model.last",sep="")) && !out.overwrite) stop(out.dir," already exists!!");
		if(!file.exists(out.dir)) dir.create(out.dir, recursive=TRUE, mode="0755");
	}
	
	if(is.x.logistic){
		if(!is.null(feature$w)) stop("When is.x.logistic is TRUE, you cannot specify feature$w")
		feature = init.obs.logistic(feature,  target="x");
	}
	if(is.y.logistic){
		if(!is.null(response$w)) stop("When is.y.logistic is TRUE, you cannot specify response$w")
		response = init.obs.logistic(response, target="y");
	}
	
	if(verbose > 0) cat("INITIALIZE THE PARAMETERS\n");
	begin.time.entire = proc.time();
	begin.time = proc.time();
	model = init.simple(
		feature=feature, response=response, nGlobalFactors=nGlobalFactors,
		param=param, nLocalFactors=nLocalFactors, 
		var_x=var_x, var_y=var_y, var_z=var_z, var_u=var_u, 
		identity.A=identity.A, identity.B=identity.B,
		is.x.logistic=is.x.logistic, is.y.logistic=is.y.logistic,
		A.sd=A.sd, B.sd=B.sd, beta.sd=beta.sd
	);
	if(!is.null(fix.var_u)) model$param$var_u = as.double(fix.var_u);
	time.used = proc.time() - begin.time;
	if(verbose > 0) cat("time used: ",time.used[3]," sec\n",sep="");
	size = check.syntax.all(feature=feature, response=response, param=model$param, factor=model$factor, check.indices=TRUE);
	if(!is.null(test.response)){
		check.syntax.all(feature=test.feature, response=test.response, param=model$param, factor=model$factor, check.indices=TRUE, test.data=TRUE);
	}
	if(verbose >= 2){
		cat("--------------------------------------------------------\n",
			"   Problem Dimensionality:\n",
			"--------------------------------------------------------\n",sep="");
		print(size);
	}
	
	if(verbose > 0) cat("START  THE EM-PROCEDURE\n");
	buffer = NULL; buffer.addr = NULL;
	
	prediction = NULL;
	best.model = NULL;  best.testLoss = Inf;
	time.pred  = NA;
	test.loss.y = rep(NA, nIter+1);  test.loss.x = rep(NA, nIter+1);
	marginal.loglik = rep(NA, nIter+1);
	
	if(!is.null(test.response)){
		begin.time = proc.time();
		prediction = predict.x.and.y(feature=test.feature, response=test.response, param=model$param, factor=model$factor);
		time.pred = proc.time() - begin.time;
		test.loss.y[1] = prediction$loss.y;
		test.loss.x[1] = prediction$loss.x;
		if(prediction$loss.y < best.testLoss){
			best.testLoss = prediction$loss.y;
			best.model = deepCopy(model);
		}
	}
	ans = output.results(
		method="EM", feature=feature, response=response, model=model, prediction=prediction,
		out.dir=out.dir, out.level=out.level,
		minTestLoss=best.testLoss, iter=0, show.marginal.loglik=show.marginal.loglik,
		TimeEM=time.used, TimeTest=time.pred, verbose=verbose,
		other=NULL, name="model"
	);
	marginal.loglik[1] = ans$marginal.loglik;
	
	for(iter in seq_len(nIter)){
		if(verbose > 0){
			cat("--------------------------------------------------------\n",
				"   EM-Iteration: ",iter,"\n",
				"--------------------------------------------------------\n",sep="");
		}
		begin.time = proc.time();
		
		# variational approx for logistic model
		if(is.x.logistic){
			feature$w = get.var.logistic(     obs=feature, param=model$param, target="x", verbose=verbose);
			feature$x = get.response.logistic(obs=feature, param=model$param, target="x", verbose=verbose);
			model$param$var_x[] = 1;
		}
		if(is.y.logistic){
			response$w = get.var.logistic(     obs=response, param=model$param, target="y", verbose=verbose);
			response$y = get.response.logistic(obs=response, param=model$param, target="y", verbose=verbose);
			model$param$var_y[] = 1;
		}
		
		if(use.C){
			model.addr = get.address(model);
			buffer = fit.EM.one.iteration.C(
				feature=feature, response=response, param=model$param, factor=model$factor,
				ridge.lambda=ridge.lambda, 
				option=option, debug=debug, verbose=verbose, buffer=buffer
			);
			# Now, model$param and model$factor contain the result after one EM iteration
			# (this is call by reference, not the behavior of regular R functions)
		
			# sanity check
			temp = get.address(model);
			if(is.diff(model.addr,temp,precision=0)) stop("model address changed!!");
			if(is.null(buffer.addr)) buffer.addr = get.address(buffer);
			temp = get.address(buffer);
			if(is.diff(buffer.addr,temp,precision=0)) stop("buffer address changed!!");
			
		}else{
			model = fit.EM.one.iteration.R(feature=feature, response=response, param=model$param)
		}
		if(!is.null(fix.var_u)) model$param$var_u = as.double(fix.var_u);
		
		# variational approx for logistic model
		if(is.x.logistic){
			mean.score = predict.x.from.z(feature=feature, param=model$param, z=model$factor$z, add.noise=FALSE);
			model$param = update.param.logistic(param=model$param, mean.score=mean.score, var.score=model$factor$var.x.score, target="x");
		}
		if(is.y.logistic){
			mean.score = predict.y.from.z(response=response, param=model$param, z=model$factor$z, add.noise=FALSE);
			model$param = update.param.logistic(param=model$param, mean.score=mean.score, var.score=model$factor$var.y.score, target="y");
		}
		
		time.used = proc.time() - begin.time;
		if(verbose > 0) cat("time used: ",time.used[3]," sec\n",sep="");
		
		if(!is.null(test.response)){
			begin.time = proc.time();
			prediction = predict.x.and.y(feature=test.feature, response=test.response, param=model$param, factor=model$factor);
			time.pred = proc.time() - begin.time;
			test.loss.y[iter+1] = prediction$loss.y;
			test.loss.x[iter+1] = prediction$loss.x;
			if(prediction$loss.y < best.testLoss){
				best.testLoss = prediction$loss.y;
				best.model = deepCopy(model);
			}
		}
		ans = output.results(
				method="EM", feature=feature, response=response, model=model, prediction=prediction,
				out.dir=out.dir, out.level=out.level,
				minTestLoss=best.testLoss, iter=iter, show.marginal.loglik=show.marginal.loglik,
				TimeEM=time.used, TimeTest=time.pred, verbose=verbose,
				other=NULL, name="model"
		);
		marginal.loglik[iter+1] = ans$marginal.loglik;
	}
		
	time.used = proc.time() - begin.time.entire;
	if(verbose > 0) cat("END OF THE EM-PROCEDURE\nTotal time used: ",time.used[3]," sec\n",sep="");
	
	out = list(model=model, model.min.test.loss=best.model, test.loss.x=test.loss.x, test.loss.y=test.loss.y, marginal.loglik=marginal.loglik);
	
	return(out);
}

###
### Initialization
###
init.simple <- function(
	feature, response, nGlobalFactors, 
	param=NULL, # directly set all the parameters
	nLocalFactors=NULL,        # which may be an array of length: nApps
	var_x=1, var_y=1, var_z=1, # which may be arrays of length: nApps
	var_u=1, identity.A=FALSE, identity.B=FALSE,
	is.y.logistic=FALSE, # whether to use a logistic model for y
	is.x.logistic=FALSE, # whether to use a logistic model for x
	A.sd=1, B.sd=1, beta.sd=1
){
	check.syntax.obs(feature=feature, response=response);
	nApps  = max(feature$app,  response$app);
	nUsers = max(feature$user, response$user);
	
	if(!is.null(param)){
		size = check.syntax.param(param);
		if(nApps != size$nApps) stop("Input param has ",param$nApps," applications, but the input data has ",nApps);
		if(nGlobalFactors != size$nGlobalFactors) stop("Input param has ",param$nGlobalFactors," global factors per user, but you specify ",nGlobalFactors);
		if(is.null(nLocalFactors)){
			nLocalFactors = size$nLocalFactors;
		}else{
			if(length(nLocalFactors) == 1) nLocalFactors = rep(nLocalFactors, nApps);
			if(length(nLocalFactors) != nApps) stop("length(nLocalFactors) != nApps");
			if(any(nLocalFactors != size$nLocalFactors)) stop("Input param has different nLocalFactors than your specification");
		}
		z = list();
		for(k in seq_len(nApps)){
			z[[k]] = matrix(0.0, nrow=nUsers, ncol=nLocalFactors[k]);
		}
		factor = list(u=matrix(0.0, nrow=nUsers, ncol=nGlobalFactors), z=z);
		out = list(param=param, factor=factor);
		return(out);
	}
	
	if(!identity.B && is.null(nLocalFactors)) stop("Please specify nLocalFactors");
	if(is.null(nLocalFactors)) nLocalFactors = rep(NA, nApps);
	
	if(length(nLocalFactors) == 1) nLocalFactors = rep(nLocalFactors, nApps);
	if(length(var_x) == 1) var_x = rep(var_x, nApps);
	if(length(var_y) == 1) var_y = rep(var_y, nApps);
	if(length(var_z) == 1) var_z = rep(var_z, nApps);
	if(length(nLocalFactors) != nApps) stop("length(nLocalFactors) != nApps");
	if(length(var_x) != nApps) stop("length(var_x) != nApps");
	if(length(var_y) != nApps) stop("length(var_y) != nApps");
	if(length(var_z) != nApps) stop("length(var_z) != nApps");
	
	A = list(); B = list(); b = list(); alpha = list(); beta = list();
	
	temp.f = tapply(seq_len(length(feature$app)),  list(feature$app),  FUN=c, simplify=FALSE);
	temp.r = tapply(seq_len(length(response$app)), list(response$app), FUN=c, simplify=FALSE);
	
	for(k in seq_len(nApps)){
		m = as.character(k); # IMPORTANT: use string to access temp.f and temp.r
		
		select = temp.f[[m]];
		obs = feature[select,];
		if(length(select) == 0){
			if(identity.B) nLocalFactors[k] = 0;
			B[[k]] = matrix(0.0, nrow=0, ncol=nLocalFactors[k]);
			b[[k]] = rep(0.0, 0);
		}else if(identity.B){
			nFeatures = max(obs$index);
			if(is.na(nLocalFactors[k])) nLocalFactors[k] = nFeatures
			else if(nLocalFactors[k] != nFeatures) stop("Misspecification of nLocalFactors[",k,"] = ",nLocalFactors[k]," with identity.B");
			B[[k]] = 1.0
			b[[k]] = rep(0.0, nFeatures);
		}else{
			nFeatures = max(obs$index);
			B[[k]] = matrix(rnorm(nFeatures*nLocalFactors[k],sd=B.sd), nrow=nFeatures, ncol=nLocalFactors[k]);
			b[[k]] = rep(0.0, nFeatures);
			agg = aggregate(obs$x, by=list(by=obs$index), FUN=mean);
			b[[k]][agg$by] = agg$x;
		}
		
		if(identity.A){
			A[[k]] = 1;
			if(nLocalFactors[k] != nGlobalFactors) stop("Misspecification of nLocalFactors[",k,"] = ",nLocalFactors[k]," with identity.A");
		}else{
			A[[k]] = matrix(rnorm(nLocalFactors[k]*nGlobalFactors,sd=A.sd), nrow=nLocalFactors[k], ncol=nGlobalFactors);	
		}
		
		select = temp.r[[m]];
		obs = response[select,];
		if(length(select) == 0){
			beta[[k]]  = matrix(0.0, nrow=0, ncol=nLocalFactors[k]);
			alpha[[k]] = rep(0.0, 0);
		}else{
			nItems = max(obs$item);
			beta[[k]]  = matrix(rnorm(nItems*nLocalFactors[k],sd=beta.sd), nrow=nItems, ncol=nLocalFactors[k]);
			alpha[[k]] = rep(0.0, nItems);
			agg = aggregate(obs$y, by=list(by=obs$item), FUN=mean);
			alpha[[k]][agg$by] = agg$x;
		}
	}
	param = list(A=A, B=B, b=b, alpha=alpha, beta=beta, var_x=var_x, var_y=var_y, var_z=var_z, var_u=var_u,
	             is.x.logistic=is.x.logistic, is.y.logistic=is.y.logistic);

	z = list();
	for(k in seq_len(nApps)){
		z[[k]] = matrix(0.0, nrow=nUsers, ncol=nLocalFactors[k]);
	}
	factor = list(u=matrix(0.0, nrow=nUsers, ncol=nGlobalFactors), z=z);
	
	if(is.x.logistic){
		param = init.param.logistic(param, feature, value=1.0, target="x");
		factor$var.x.score = rep(1.0, nrow(feature));
	}
	if(is.y.logistic){
		param = init.param.logistic(param, response, value=1.0, target="y");
		factor$var.y.score = rep(1.0, nrow(response));
	}	
	
	out = list(param=param, factor=factor);
	return(out);
}

###
### One EM iteration
###		Run E-Step once, M-step once
###	Input:
###		feature = data.frame(user, app, index, x, w)
###    response = data.frame(user, app, item,  y, w)
###       param = list(A, B, b, alpha, beta, var_x, var_y, var_z, var_u)
### Output: list(param, factor)
###		param: The parameter values after the M-step
###    factor=list(u, z): The posterior mean of u and z
###    if(is.x.logistic) factor$var.x.score[ikm] = Var[B_{k,m} z_{ik}]  
###    if(is.y.logistic) factor$var.y.score[ijk] = Var[beta_{jk} z_{ik}]
fit.EM.one.iteration.R <- function(feature, response, param){
	stop("Please implement this function");
}
### IMPORTANT NOTE:
###   * In the C version, to reduce the memory footprint, the content of the input 
###     param and factor will be changed to the updated values.
###     This is call by reference, and is NOT the behavior of regular R functions.
###     buffer is the temp space (which can be NULL) and is also call by reference.
###   * ridge.lambda[c("A", "B", "beta")] are the numbers to be added to the
###     diagonal when fitting the regression to get A, B, and beta
###   * buffer = list(A, B, b, alpha, beta, z)
###     contains the packed version of A, B, b, ... for C/C++ functions
fit.EM.one.iteration.C <- function(
	feature, response, param, factor,
	ridge.lambda=c(A=0, B=0, beta=0),
	buffer=NULL, option=0, debug=0, verbose=0
){
	size = check.syntax.all(feature=feature, response=response, param=param, factor=factor);
	
	if(is.null(buffer)) buffer = list()
	else                check_names(buffer, "buffer", required=c("A", "B", "b", "alpha", "beta", "z"));

	if(length(ridge.lambda) != 3) stop("length(ridge.lambda) != 3");
	if(any(c("A", "B", "beta") != names(ridge.lambda))) stop("ridge.lambda must be a named vector with names: A, B, beta (in this order)");
	
	for(name in c("A", "B", "b", "alpha", "beta")){
		if(is.null(buffer[[name]])) buffer[[name]] = pack.list.of.matrices(param[[name]])
		else                  pack.list.of.matrices(param[[name]], output=buffer[[name]]);
		check_type_size(buffer[[name]]$data, "double", NA);
		check_type_size(buffer[[name]]$dim,  "int",    NA);
	}
	for(name in c("z")){
		if(is.null(buffer[[name]])) buffer[[name]] = pack.list.of.matrices(factor[[name]])
		else                  pack.list.of.matrices(factor[[name]], output=buffer[[name]]);
		check_type_size(buffer[[name]]$data, "double", NA);
		check_type_size(buffer[[name]]$dim,  "int",    NA);
	}
	for(name in c("var_x", "var_y", "var_z")){
		check_type_size(param[[name]], "double", size$nApps);
	}
	check_type_size(param[["var_u"]], "double", 1);
	check_type_size(factor[["u"]], "double", c(size$nUsers, size$nGlobalFactors));

	for(name in c("user", "app", "index")){
		check_type_size(feature[[name]], "int", nrow(feature));
	}
	check_type_size(feature[["x"]], "double", nrow(feature));
	if(is.null(feature[["w"]])){
		feature_has_w = as.integer(0);
	}else{
		feature_has_w = as.integer(1);
		check_type_size(feature[["w"]], "double", nrow(feature));
	}

	for(name in c("user", "app", "item")){
		check_type_size(response[[name]], "int", nrow(response));
	}
	check_type_size(response[["y"]], "double", nrow(response));
	if(is.null(response[["w"]])){
		response_has_w = as.integer(0);
	}else{
		response_has_w = as.integer(1);
		check_type_size(response[["w"]], "double", nrow(response));
	}
	
	for(name in c("nFeatures", "nItems", "nLocalFactors")){
		check_type_size(size[[name]], "int", size$nApps)
	}
	
	check_type_size(factor[["var.y.score"]], "double", length(factor[["var.y.score"]]));
	check_type_size(factor[["var.x.score"]], "double", length(factor[["var.x.score"]]));
	
	.C("EM_one_iteration",
		# Parameters: INPUT and OUTPUT
		buffer[["A"]]$data, buffer[["A"]]$dim,  buffer[["B"]]$data, buffer[["B"]]$dim,  buffer[["b"]]$data, buffer[["b"]]$dim,
		buffer[["alpha"]]$data, buffer[["alpha"]]$dim,  buffer[["beta"]]$data, buffer[["beta"]]$dim,
		param[["var_x"]], param[["var_y"]], param[["var_z"]], param[["var_u"]],
		# Posterior mean of factors: OUTPUT
		factor[["u"]],  buffer[["z"]]$data, buffer[["z"]]$dim,
		# Posterior variance for logistic regression: OUTPUT
		factor[["var.x.score"]], as.integer(length(factor[["var.x.score"]])),
		factor[["var.y.score"]], as.integer(length(factor[["var.y.score"]])),
		# Feature table: INPUT
		feature[["user"]], feature[["app"]], feature[["index"]], feature[["x"]], feature[["w"]],
		as.integer(nrow(feature)), feature_has_w,
		# Response table: INPUT
		response[["user"]], response[["app"]], response[["item"]], response[["y"]], response[["w"]],
		as.integer(nrow(response)), response_has_w,
		# Ridge regression parameters: INPUT
		as.double(ridge.lambda), as.integer(length(ridge.lambda)),
		# Size information: INPUT
		as.integer(size$nApps), as.integer(size$nUsers), as.integer(size$nGlobalFactors),
		size$nFeatures, size$nItems, size$nLocalFactors,
		# Others
		as.integer(option), as.integer(verbose), as.integer(debug),
		DUP=FALSE
	);
	
	for(name in c("A", "B", "b", "alpha", "beta")){
		unpack.list.of.matrices(buffer[[name]], output=param[[name]]);
	}
	unpack.list.of.matrices(buffer[["z"]], output=factor[["z"]]);
	
	return(buffer);
}

