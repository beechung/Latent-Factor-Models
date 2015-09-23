### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
###
### Author: Bee-Chung Chen

### ---------------------------------------------------------------------------
###          MODEL
### ---------------------------------------------------------------------------
###    Observation: (src.id, dst.id, src.context, dst.context, edge.context, y);
###                       i       j          k_s,         k_d,          k_e
###	         y[ijk] ~ N(fScore[ijk] + x_obs[ijk,]'b,  var_y[ijk])
###     fScore[ijk] = alpha[i,k_s] + beta[j,k_d] + gamma[k_e] + sum(v[i,] * v[j,] * w[k_e,])  if u is null
###                                                       ... + sum(u[i,] * v[j,] * w[k_e,])  if u is not null
###      alpha[i,k] ~ N(q[k]*alpha_global[i] + x_src[i,]'g0[,k],  var_alpha[k])
### alpha_global[i] ~ N(0,                                        var_alpha_global=1)
###       beta[j,k] ~ N(r[k]*beta_global[j]  + x_dst[j,]'d0[,k],  var_beta[k])
###  beta_global[j] ~ N(0,                                        var_beta_global=1)
###	       gamma[k] ~ N(x_ctx[k,]'h0, var_gamma)
###           u[i,] ~ N(x_src[i,]'G,  var_u)
###           v[j,] ~ N(x_dst[j,]'D,  var_v)
###           w[k,] ~ N(x_ctx[k,]'H,  var_w)
###
### SPECIAL CASE: param$nLocalFactors is not NULL
###	   * Model: fScore[ijk] = alpha[i,k_s] + beta[j,k_d] + gamma[k_e] + sum(v[i,,k_e] * v[j,,k_e]) if u is NULL
###	                                                              ... + sum(u[i,,k_e] * v[j,,k_e]) if u is not NULL
###                 u[i,,k] ~ N(x_src[i,]' G[,,k],  var_u[k])
###                 v[j,,k] ~ N(x_dst[j,]' D[,,k],  var_v[k])
###    * Data setup:
###        nFactors = nLocalFactors*nEdgeContexts
###        factor$u = matrix(u, nrow=nSrcNodes, ncol=nLocalFactors*nEdgeContexts)
###        factor$v = matrix(v, nrow=nDstNodes, ncol=nLocalFactors*nEdgeContexts)
###        factor$w = array(0.0, dim=c(nEdgeContexts,nFactors));
###        for(k in 1:nEdgeContexts) factor$w[k, (k-1)*nLocalFactors + (1:nLocalFactors) ] = 1;
###        param$var_u = rep(0.0, nFactors);
###        for(k in 1:nEdgeContexts) param$var_u[ (k-1)*nLocalFactors + (1:nLocalFactors) ] = var_u[k];
###        param$var_v = rep(0.0, nFactors);
###        for(k in 1:nEdgeContexts) param$var_v[ (k-1)*nLocalFactors + (1:nLocalFactors) ] = var_v[k];
###        param$var_w = 0;
###
###   NOTE: (1) If you want an INTERCEPT in regression, please add the intercept column in
###             feature${x_obs, x_src, x_dst, x_ctx} by yourself.
###             You can also set feature$x_obs[] = 0 to obtain a zero-mean model.
###         (2) To DISABLE fitting a FACTOR, set param$var_FACTOR == NULL. Both sampling for the
###             FACTOR and regression for the FACTOR will be disabled.
###         (3) Object format for data.train and data.test:
###             data.train = list(
###		           obs = data.frame(src.id, dst.id, src.context, dst.context, edge.context, y),
###                feature = list(x_src, x_dst, x_ctx, x_obs)
###             )
###             data.test has the same format.
###
### IMPORTANT NOTE:
###   (1) If we predict using E[u_i]E[v_j], we should also use E[u_i]E[v_j] in fitting (instead of E[u_i v_j])
###
### OUTPUT:
###   output$model.last:        Parameter estimates at the end of the last iteration
###   output$model.minTestLoss: Parameter estimates with the lowest test-set loss computed using the test dataset
###                             if is.null(test.obs), then this output will NOT be available.
###
fit.multicontext <- function(
	init.model, # Initial model = list(factor, param)
	nSamples,   # Number of samples drawn in each E-step: could be a vector of size nIter.
	nBurnIn,    # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
	nIter=NULL, # Number of EM iterations
	data.train=NULL, # Training data
	data.test=NULL,  # Test data (optional)
	is.logistic=FALSE,
	out.level=0,  # out.level=1: Save the factor & parameter values to out.dir/model.last and out.dir/model.minTestLoss
	out.dir=NULL, # out.level=2: Save the factor & parameter values of each iteration i to out.dir/model.i
	out.overwrite=FALSE,
	debug=0,      # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
	verbose=0,    # Set to 0 to disable console output; Set to 100 to print everything to the console
	verbose.E=verbose,
	verbose.M=verbose,
	use.C=TRUE,   # Whether to use the C implementation (R implementation does not have full functionalities)
	error.handler=stop, # You can change it to warning so that the program won't stop
	output.at.end.of.EStep=FALSE,
	rm.factors.without.obs.in.loglik=TRUE,
	ridge.lambda=c(b=1, g0=1, d0=1, h0=1, G=1, D=1, H=1), # (or a scalar) Add diag(lambda) to X'X in linear regression
	approx.interaction=TRUE, # predict E[uv] as E[u]E[v].
	zero.mean=rep(0,0),  # zero.mean["alpha"] = TRUE  ->  g = 0, etc.
	fix.var=NULL,        # fix.var[["u"]] = n -> var_u = n (NULL -> default, list() -> fix no var)
	max.nObs.for.b=NULL, # maximum number of observations to be used to fit b
	# The following five are for backward competibility when data.train=NULL and/or data.test=NULL
	IDs=NULL,
	obs=NULL,          # Training data: Observation table
	feature=NULL,      #                Features
	test.obs=NULL,     # Test data: Observations for testing
	test.feature=NULL  #            Features for testing
){
	factor = init.model$factor;
	param  = init.model$param;
	param$approx.interaction = approx.interaction;
	if(approx.interaction) test.obs.for.Estep = NULL
	else                   test.obs.for.Estep = test.obs;
	param$is.logistic = is.logistic;

	# setup obs, feature, test.obs, test.feature
	if(!is.null(data.train)){
		if(!all(c("obs", "feature") %in% names(data.train))) stop("Please check input parameter 'data.train' when calling function fit.multicontext or run.multicontext: data.train$obs and data.train$feature cannot be NULL");
		if(!is.null(obs)) stop("When calling function fit.multicontext or run.multicontext, if you already specified 'data.train', then you should set 'obs=NULL'");
		if(!is.null(feature)) stop("When calling function fit.multicontext or run.multicontext, if you already specified 'data.train', then you should set 'feature=NULL'");
		obs=data.train$obs;
		feature=data.train$feature;
		data.train$obs = NULL;
		data.train$feature = NULL;
	}else{
		if(is.null(obs) || is.null(feature)) stop("Please specify input parameter 'data.train' when calling function fit.multicontext or run.multicontext");
	}
	if(!is.null(data.test)){
		if(!all(c("obs", "feature") %in% names(data.test))) stop("Please check input parameter 'data.test' when calling function fit.multicontext or run.multicontext: data.test$obs and data.test$feature cannot be NULL");
		if(!is.null(test.obs)) stop("When calling function fit.multicontext or run.multicontext, if you already specified 'data.test', then you should set 'test.obs=NULL'");
		if(!is.null(test.feature)) stop("When calling function fit.multicontext or run.multicontext, if you already specified 'data.test', then you should set 'test.feature=NULL'");
		test.obs=data.test$obs;
		test.feature=data.test$feature;
		if(is.null(IDs)) IDs = data.test$IDs;
	}else{
		if(( is.null(test.obs) && !is.null(test.feature)) ||
		   (!is.null(test.obs) &&  is.null(test.feature)))
	  		stop("If you want to supply test data to the fitting code, please specify input parameter 'data.test' when calling function fit.multicontext or run.multicontext");
	}

	# Sanity check
	if(out.level > 0 && is.null(out.dir)) stop("Please specify input parameter 'out.dir' when calling function fit.multicontext or run.multicontext with out.level > 0");
	if(out.level > 0 && file.exists(out.dir) && !out.overwrite) stop("Output directory '",out.dir,"' EXISTS!!  Please remove the directory or specify a new directory for the input parameter out.dir.");

	has.u = TRUE;
	if(is.null(factor$u)){
		has.u = FALSE;
		if(!is.null(factor$v) && any(obs$src.id == obs$dst.id)) stop("When init.model$factor$u is disabled (i.e., init.model$factor$u = NULL), self ratings are not allowed. Please remove the rows in the input parameter 'obs' that have obs$src.id == obs$dst.id.");
		if(!is.null(param$var_u)) stop("When init.model$factor$u is NULL, init.model$param$var_u must also be NULL. Please set init.model$param$var_u to NULL.");
	}
	if(is.null(nIter) && length(nSamples) == 1) stop("Please specify input parameter 'nIter' when calling function fit.multicontext or run.multicontext.");
	if(is.null(nIter)) nIter = length(nSamples);
	if(length(nSamples)!=nIter){
		if(length(nSamples) != 1) stop("Please check input parameters 'nSamples' and 'nIter' when calling function fit.multicontext or run.multicontext. When nSamples is a vector of more than one element, the length of the vector must equal nIter (i.e., length(nSamples) == nIter). In this case, nSamples[i] specifies the number of Gibbs samples in the i-th EM iteration. You can just set nIter=NULL to fix the problem (this would force the number of EM iteration = length(nSamples)).");
		nSamples = rep(nSamples,nIter);
	}
	if(length(nBurnIn)!=nIter){
		if(length(nBurnIn) != 1) stop("Please check input parameter 'nBurnIn' when calling function fit.multicontext or run.multicontext. nBurnIn should either be a scalar (i.e., length(nBurnIn) == 1) or a vector with length equal to the length of nSamples (i.e., length(nBurnIn) = length(nSamples))");
		nBurnIn = rep(nBurnIn,nIter);
	}
	if(!is.null(test.obs)){
		if(is.null(test.feature)) stop("Please check input parameters 'test.obs' and 'test.feature' when calling function fit.multicontext or run.multicontext. Notice that test.obs != NULL, but test.feature is NULL. Please either set test.obs = NULL or specify test.feature.");
	}
	if(is.null(names(ridge.lambda)) && length(ridge.lambda) == 1){
		ridge.lambda = rep(ridge.lambda, 7);
		names(ridge.lambda) = c("b", "g0", "d0", "h0", "G", "D", "H");
	}
	if(is.null(fix.var)){
		# default setting
		fix.var = list();
		if(!is.null(factor$u)) fix.var[["u"]] = 1;
		if(!is.null(factor$w)){
			if(is.null(param$nLocalFactors)) fix.var[["w"]] = 1
			else                             fix.var[["w"]] = 0;
		}
	}
	if(!is.null(param$nLocalFactors)){
		nEdgeContexts = ncol(factor$v)/param$nLocalFactors;
		if(is.null(feature$x_ctx)) feature$x_ctx = array(0.0, dim=c(nEdgeContexts,0));
		if(!is.null(test.feature) && is.null(test.feature$x_ctx)) test.feature$x_ctx = array(0.0, dim=c(nEdgeContexts,0));
	}

	if(!all(names(ridge.lambda) %in% c("b","g0","d0","h0","G","D","H")))  stop("Please check input parameter 'ridge.lambda' when calling function fit.multicontext or run.multicontext. names(ridge.lambda) is not specified correctly. The only allowable names are 'b', 'g0', 'd0', 'h0', 'G', 'D', 'H'.");
	if(!all(names(zero.mean) %in% c("alpha","beta","gamma","u","v","w"))) stop("Please check input parameter 'zero.mean' when calling function fit.multicontext or run.multicontext. names(zero.mean) is not specified correctly. The only allowable names are 'alpha', 'beta', 'gamma', 'u', 'v', 'w'.");
	if(!all(names(fix.var) %in% c("alpha","beta","gamma","u","v","w")))   stop("Please check input parameter 'fix.var' when calling function fit.multicontext or run.multicontext. names(fix.var) is not specified correctly. The only allowable names are 'alpha', 'beta', 'gamma', 'u', 'v', 'w'.");

	warning.any.not.in(c("alpha", "beta"), names(factor), "Please check input parameter 'init.model' when calling function fit.multicontext or run.multicontext. The following components in init.model$factor are required: ", stop=TRUE);

	# Initialize obs for the logistic model
	obs      = init.obs(obs=obs,      is.logistic=is.logistic);
	test.obs = init.obs(obs=test.obs, is.logistic=is.logistic);
	# now obs$response stores the input response values

	if(!is.integer(obs$src.id)) obs$src.id = as.integer(obs$src.id);
	if(!is.integer(obs$dst.id)) obs$dst.id = as.integer(obs$dst.id);
	if(!is.null(obs$src.context)  && !is.integer(obs$src.context))  obs$src.context  = as.integer(obs$src.context);
	if(!is.null(obs$dst.context)  && !is.integer(obs$dst.context))  obs$dst.context  = as.integer(obs$dst.context);
	if(!is.null(obs$edge.context) && !is.integer(obs$edge.context)) obs$edge.context = as.integer(obs$edge.context);

	subset.info = NULL;
	if(rm.factors.without.obs.in.loglik) subset.info = get.subset.info(obs, output.any=!has.u);

	size = syncheck.multicontext.spec(factor=factor, obs=obs, feature=feature, param=param,
			warning=10, print=TRUE);
	check.obs.feature(obs, feature, nSrcContexts=size$nSrcContexts, nDstContexts=size$nDstContexts, nEdgeContexts=size$nEdgeContexts);
	if(!is.null(test.obs)){
		nEdgeContexts = NA;
		if(!is.null(param$nLocalFactors)) nEdgeContexts = size$nEdgeContexts;
		check.obs.feature(test.obs, test.feature, nSrcContexts=size$nSrcContexts, nDstContexts=size$nDstContexts, nEdgeContexts=nEdgeContexts);
	}
	if(!is.null(test.obs.for.Estep)){
		temp = feature; temp$x_obs = test.feature$x_obs;
		check.obs.feature(test.obs.for.Estep, temp, nSrcContexts=size$nSrcContexts, nDstContexts=size$nDstContexts, nEdgeContexts=size$nEdgeContexts);
	}

	factor = deepCopy(factor); # Make a copy (if not, the C code will modify the input factor values)

	if(verbose >= 1) cat("================= START fit.MCEM =====================================\n",sep="");
	if(verbose >= 2) print(size);

	begin.time = proc.time();

	CD.loglik = rep(NA, nIter+1);
	E.loglik  = rep(NA, nIter+1);
	loglik = get.logLikelihood(obs, factor, feature, param, subset.info=subset.info, verbose=verbose, is.logistic=is.logistic, prefix=" Initial");
	CD.loglik[1] = loglik$CD;
	E.loglik[1]  = loglik$E;

	prediction = NULL;  TestLoss = NULL;  minTestLoss = NULL;  model.minTestLoss = NULL;
	if(!is.null(test.obs)){
		TestLoss = rep(NA, nIter+1); # TestLoss records the loss in the test set of each iteration
		prediction = predict.multicontext(model=list(factor=factor, param=param), obs=test.obs, feature=test.feature, is.logistic=is.logistic);
		minTestLoss = prediction$test.loss;
		TestLoss[1] = minTestLoss;
		# Factor & paramter estimates with min TestLoss
		model.minTestLoss = list(factor=factor, param=param);
	}
	time.used = proc.time() - begin.time;

	if(verbose >= 1 && !is.null(test.obs)) cat("      test loss:    ", minTestLoss, " (",time.used[3]," sec)\n",sep="");

	output.to.dir(
			out.dir=out.dir, factor=factor, param=param, IDs=IDs,
			prediction=prediction, loglik=loglik$CD,
			minTestLoss=minTestLoss, nSamples=nSamples, iter=0, out.level=out.level, out.overwrite=out.overwrite,
			TimeEStep=0, TimeMStep=0, TimeTest=time.used[3], verbose=verbose, name="model"
	);

	for(iter in 1:nIter){

		# Create response (for the logistic model)
		response = generate.response(obs=obs, param=param, is.logistic=is.logistic, verbose=verbose);
		obs$y = response$y;
		param$var_y = response$var_y;

		if(verbose >= 1){
			cat("---------------------------------------------------------\n",
				"        Iteration ",iter,"\n",
				"---------------------------------------------------------\n",
				"start E-STEP\n",sep="");
		}
		b.time = proc.time();

		###
		### E-STEP
		###
		if(use.C){
			mc_e = MCEM_EStep.multicontext.C(
					factor=factor, obs=obs, feature=feature, param=param,
					nSamples=nSamples[iter], nBurnIn=nBurnIn[iter],
					test.obs=test.obs.for.Estep,
					debug=debug, verbose=verbose.E
			);
			# IMPORTANT NOTE: factor is call-by-reference in MCEM_EStep.multicontext.C
			# Now, factor = mc_e$mean with the same memory address!!!
		}else{
			mc_e = MCEM_EStep.R(
					factor, obs, feature, param, nSamples[iter], nBurnIn=nBurnIn[iter],
					debug=debug, verbose=verbose.E
			);
		}
		factor = mc_e$mean;

		time.used.1 = proc.time() - b.time;

		# verbose
		if(verbose >= 1 || output.at.end.of.EStep){
			if(verbose >= 1) cat("end   E-STEP (used ",time.used.1[3]," sec)\n", sep="");
			loglik.E = get.logLikelihood(obs, factor, feature, param, factor.var=mc_e$var, factor.cov=mc_e$cov, subset.info=subset.info, verbose=verbose, is.logistic=is.logistic, prefix="        ");
			if(!is.null(test.obs) && (verbose >= 2 || output.at.end.of.EStep)){
				b.time.TestLoss = proc.time();
				prediction = predict.multicontext(model=list(factor=factor, param=param), obs=test.obs, feature=test.feature, is.logistic=is.logistic, fScore=mc_e$test.fScore$mean);
				time.used.TestLoss = proc.time() - b.time.TestLoss;
				if(verbose >= 2) cat("      test loss:    ", prediction$test.loss, " (",time.used.TestLoss[3]," sec)\n",sep="");
			}
		}

		if(output.at.end.of.EStep){
			output.to.dir(
				out.dir=out.dir, factor=factor, param=param, IDs=IDs,
				prediction=prediction, loglik=loglik.E$CD,
				minTestLoss=minTestLoss, nSamples=nSamples, iter=iter, out.level=out.level, out.overwrite=out.overwrite,
				TimeEStep=time.used.1[3], TimeMStep=0, TimeTest=time.used.TestLoss[3], verbose=verbose,
				other=list(mc_e=mc_e), name="model-end-of-E"
			);
		}

		###
		### M-STEP
		###
		b.time = proc.time();
		if(verbose >= 1) cat("start M-STEP\n",sep="");
		param.new = MCEM_MStep.multicontext(
				factor.mean=mc_e$mean, factor.var=mc_e$var, factor.cov=mc_e$cov,
				obs=obs, feature=feature, param=param, subset.info=subset.info,
				ridge.lambda=ridge.lambda, zero.mean=zero.mean, fix.var=fix.var,
				is.logistic=is.logistic, max.nObs.for.b=max.nObs.for.b,
				debug=debug, verbose=verbose.M
		);
		param = update.param(factor.mean=mc_e$mean, factor.var=mc_e$var, param=param.new, x_obs=feature$x_obs, obs=obs, is.logistic=is.logistic);

		time.used.2 = proc.time() - b.time;

		if(verbose >= 1) cat("end   M-STEP (used ",time.used.2[3]," sec)\n", sep="");

		b.time.test = proc.time();

		loglik = get.logLikelihood(obs, factor, feature, param, factor.var=mc_e$var, factor.cov=mc_e$cov, subset.info=subset.info, verbose=verbose, is.logistic=is.logistic, prefix="        ");
		CD.loglik[iter+1] = loglik$CD;
		E.loglik[iter+1]  = loglik$E;
		if(verbose >= 1){
			if(loglik$E < loglik.E$E){
				cat("\nNOTICE: E[loglik] decreases after the M step (may be due to regularized regression)\n\n");
				if(all(ridge.lambda == 0) && !approx.interaction) stop("E[loglik] decreases after the M step");
				if(debug > 0) warning("E[loglik] decreases after the M step ",loglik.E$E," -> ",loglik$E," (may be due to regularized regression)");
			}
		}

		if(!is.null(test.obs)){
			b.time.TestLoss = proc.time();
			prediction = predict.multicontext(model=list(factor=factor, param=param), obs=test.obs, feature=test.feature, is.logistic=is.logistic, fScore=mc_e$test.fScore$mean);
			TestLoss[iter+1] = prediction$test.loss;
			time.used.TestLoss = proc.time() - b.time.TestLoss;
		}
		time.used.3 = proc.time() - b.time.test;

		if(verbose >= 1){
			cat("  training loss:       ", attr(loglik$CD,"loss"), "\n",sep="");
			if(!is.null(test.obs))
			cat("      test loss:       ", prediction$test.loss, " (",time.used.TestLoss[3]," sec)\n",sep="");
		}

		###
		### Update the model.minTestLoss model if the TestLoss decreases
		###
		if(!is.null(test.obs) && TestLoss[iter+1] < minTestLoss){
			minTestLoss = TestLoss[iter+1];
			model.minTestLoss$factor=deepCopy(factor);
			model.minTestLoss$param=param;

			if(verbose >= 2) cat("TestLoss decreases!\n");
		}

		output.to.dir(
				out.dir=out.dir, factor=factor, param=param, IDs=IDs,
				prediction=prediction, loglik=loglik$CD,
				minTestLoss=minTestLoss, nSamples=nSamples, iter=iter, out.level=out.level, out.overwrite=out.overwrite,
				TimeEStep=time.used.1[3], TimeMStep=time.used.2[3], TimeTest=time.used.3[3], verbose=verbose, name="model",
				data.train=data.train
		);
	}

	output = list(
			CD.loglik   = CD.loglik,
			E.loglik    = E.loglik,
			test.loss   = TestLoss,
			minTestLoss = minTestLoss,
			model.minTestLoss = model.minTestLoss,
			model.last = list(factor=factor, param=param)
	);

	time.used = proc.time() - begin.time;
	if(verbose >= 1){
		if(!is.null(test.obs)) cat("Minimum test-set loss: ",minTestLoss,"\n",sep="");
		cat("Total time: ",time.used[3]," sec.\n",
			"=================== END fit.MCEM =====================================\n",sep="");
	}

	return(output);
}
# run the model including initialization
#   setting = data.frame(name, nFactors, has.u, has.gamma, nLocalFactors, is.logistic)
#                                        T/F    T/F        0 or n         T/F
#   Output dir is out.dir_name
#   data.train = list(
#		obs = data.frame(src.id, dst.id, src.context, dst.context, edge.context, y),
#       feature = list(x_src, x_dst, x_ctx, x_obs)
#   )
#   data.test has the same format
#
run.multicontext <- function(
	data.train=NULL, # Training data
	setting,    # Model setting
	nSamples,   # Number of samples drawn in each E-step: could be a vector of size nIter.
	nBurnIn,    # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
	nIter=NULL, # Number of EM iterations
	data.test=NULL, # Test data (optional)
	return.models=FALSE,
	approx.interaction=TRUE, # predict E[uv] as E[u]E[v].
	reg.algo=NULL,     # The regression algorithm to be used in the M-step (NULL => linear regression)
	reg.control=NULL,  # The control paramter for reg.algo
	# initialization parameters
	var_alpha=1, var_beta=1, var_gamma=1,
	var_v=1, var_u=1, var_w=1, var_y=NULL,
	relative.to.var_y=FALSE, var_alpha_global=1, var_beta_global=1,
	# others
	out.level=0,  # out.level=1: Save the factor & parameter values to out.dir/model.last and out.dir/model.minTestLoss
	out.dir=NULL, # out.level=2: Save the factor & parameter values of each iteration i to out.dir/model.i
	out.overwrite=FALSE,
	debug=0,      # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
	verbose=0,    # Set to 0 to disable console output; Set to 100 to print everything to the console
	verbose.E=verbose,
	verbose.M=verbose,
	use.C=TRUE,   # Whether to use the C implementation (R implementation does not have full functionalities)
	error.handler=stop, # You can change it to warning so that the program won't stop
	output.at.end.of.EStep=FALSE,
	rm.factors.without.obs.in.loglik=TRUE,
	ridge.lambda=1, # Add diag(lambda) to X'X in linear regression
	zero.mean=rep(0,0), # zero.mean["alpha"] = TRUE  ->  g = 0, etc.
	fix.var=NULL,       # fix.var[["u"]] = n -> var_u = n (NULL -> default, list() -> fix no var)
	max.nObs.for.b=NULL,# maximum number of observations to be used to fit b
	rnd.seed.init=NULL, rnd.seed.fit=NULL,
	# The following five are for backward competibility when data.train=NULL and/or data.test=NULL
	IDs=NULL,
	obs=NULL,          # Training data: Observation table
	feature=NULL,      #                Features
	test.obs=NULL,     # Test data: Observations for testing
	test.feature=NULL  #            Features for testing
){
	if(length(unique(setting$name)) != nrow(setting)) stop("Please check input parameter 'setting' when calling function run.multicontext: setting$name must be a column of unique identifiers");
	if(!out.overwrite){
		for(k in 1:nrow(setting)){
			temp = paste(out.dir,"_",setting[k,"name"],sep="");
			file = paste(temp,"/model.last",sep="");
			if(file.exists(file)) stop("Please check input parameters 'out.overwrite' and 'out.dir' when calling function run.multicontext. Directory '",temp,"' already exists. Please either remove the directory or set out.overwrite = TRUE to overwrite the existing directory.");
		}
	}
	if(!all(setting$has.u %in% c(TRUE, FALSE))) stop("Please check input parameter 'setting' when calling function run.multicontext: setting$has.u should be either TRUE or FALSE.");
	if(!all(setting$has.gamma %in% c(TRUE, FALSE))) stop("Please check input parameter 'setting' when calling function run.multicontext: setting$has.gamma should be either TRUE or FALSE.");
	if(!all(setting$is.logistic %in% c(TRUE, FALSE))) stop("Please check input parameter 'setting' when calling function run.multicontext: setting$is.logistic should be either TRUE or FALSE.");
	if(!all(setting$nFactors >= 0)) stop("Please check input parameter 'setting' when calling function run.multicontext: setting$nFactors = ",setting$nFactors,", which should be >= 0.");
	if(!all(setting$nLocalFactors >= 0)) stop("Please check input parameter 'setting' when calling function run.multicontext: setting$nLocalFactors = ",setting$nLocalFactors,", which should be >= 0.");
	nEdgeContexts = if(is.null(obs$edge.context)) 0 else max(obs$edge.context);
	if(any(setting$nLocalFactors != 0 & setting$nFactors != setting$nLocalFactors * nEdgeContexts)) stop("Please check input parameter 'setting' when calling function run.multicontext: setting$nFactors must = setting$nLocalFactors * max(obs$edge.context).");
	if(any(!setting$has.u) && !is.null(IDs$SrcIDs) && !is.null(IDs$DstIDs)){
		if(length(IDs$SrcIDs) != length(IDs$DstIDs) || any(IDs$SrcIDs != IDs$DstIDs))
			stop("Please check input parameters 'setting' and 'obs' when calling function run.multicontext: Some row of setting$has.u is FALSE, but obs is indexed in a way that does not support this option. Please call function indexData with src.dst.same=TRUE before calling function run.multicontext.");
	}

	out = list(summary=setting);
	out$summary$best.CD.loglik = NA;
	out$summary$last.CD.loglik = NA;
	out$summary$best.test.loss = NA;
	out$summary$last.test.loss = NA;
	out$summary$time.used = NA;
	out$test.loss = list();
	out$CD.loglik = list();
	if(return.models) out$model = list();
	for(k in 1:nrow(setting)){
		name = as.character(setting[k,"name"]);
		if(!is.null(rnd.seed.init)) set.seed(rnd.seed.init);
		begin.time = proc.time();
		init = init.simple.random(
			data.train=data.train, obs=obs, feature=feature,
			nFactors=setting[k,"nFactors"], has.u=setting[k,"has.u"], has.gamma=setting[k,"has.gamma"],
			nLocalFactors=setting[k,"nLocalFactors"], is.logistic=setting[k,"is.logistic"],
			var_alpha=var_alpha, var_beta=var_beta, var_gamma=var_gamma,
			var_v=var_v, var_u=var_u, var_w=var_w, var_y=var_y,
			relative.to.var_y=relative.to.var_y, var_alpha_global=var_alpha_global, var_beta_global=var_beta_global
		);
		init$param$reg.algo = reg.algo;
		init$param$reg.control = reg.control;

		if(!is.null(rnd.seed.fit)) set.seed(rnd.seed.fit);
		ans = fit.multicontext(
			data.train=data.train, data.test=data.test,
			obs=obs, feature=feature, init.model=init, nSamples=nSamples, nBurnIn=nBurnIn, nIter=nIter,
			test.obs=test.obs, test.feature=test.feature,
			IDs=IDs, is.logistic=setting[k,"is.logistic"],
			out.level=out.level, out.dir=paste(out.dir,"_",name,sep=""), out.overwrite=out.overwrite,
			debug=debug,
			verbose=verbose, verbose.E=verbose.E, verbose.M=verbose.M,
			use.C=use.C, error.handler=error.handler,
			output.at.end.of.EStep=output.at.end.of.EStep,
			rm.factors.without.obs.in.loglik=rm.factors.without.obs.in.loglik,
			ridge.lambda=ridge.lambda, approx.interaction=approx.interaction,
			zero.mean=zero.mean, fix.var=fix.var, max.nObs.for.b=max.nObs.for.b
		);
		time.used = proc.time() - begin.time;
		out$summary$time.used[k] = time.used[3];
		out$summary$best.CD.loglik[k] = max(ans$CD.loglik);
		out$summary$last.CD.loglik[k] = ans$CD.loglik[length(ans$CD.loglik)];
		out$CD.loglik[[name]] = ans$CD.loglik;
		if(!is.null(ans$test.loss)){
			out$summary$best.test.loss[k] = min(ans$test.loss);
			out$summary$last.test.loss[k] = ans$test.loss[length(ans$test.loss)];
			out$test.loss[[name]] = ans$test.loss;
		}
		if(return.models) out$model[[name]] = ans$model.last;
	}
	return(out);
}


###
### MCEM_EStep.C. See MCEM_EStep_multicontext(...) in src/C/factor_model_multicontext.c
###
### INPUT:  factor  = list(alpha, beta, gamma, u, v, w); # Initial factor values
###         obs     = data.frame(src.id, dst.id, src.context, dst.context, edge.context, y);
###         feature = list(x_obs, x_src, x_dst, x_ctx);
###         param   = list(b, g0, d0, h0, G, D, H, q, r,
###                        var_y, var_alpha, var_alpha_global, var_beta_global, var_beta, var_gamma, var_u, var_v, var_w);
###
### OUTPUT: mean = list(alpha, alpha_global, beta, gamma, beta_global, u, v, w, fScore);
###         var  = list(alpha, alpha_global, beta, gamma, beta_global, u, v, w, fScore);
###         cov  = list(alpha, beta);
###         test.fScore = list(mean, var) # fScore for test cases
###
### NOTE: factor$alpha and factor$beta cannot be NULL!!
###       Set factor$gamma = NULL to disable gamma
###	      All src nodes, dst nodes, etc., in test.obs must be in factor and feature.
###
MCEM_EStep.multicontext.C <- function(
	factor, obs, feature, param, nSamples, nBurnIn=1, test.obs=NULL,
	debug=0, verbose=0
){
	size = syncheck.multicontext.spec(factor=factor, obs=obs, feature=feature, param=param);

	if(is.null(test.obs))             nTestObs = as.integer(0)
	else if(is.data.frame(test.obs))  nTestObs = as.integer(nrow(test.obs))
	else stop("test.obs should be either NULL or a data frame.");

	if(is.null(factor$alpha)) stop("factor$alpha cannot be null");
	if(is.null(factor$beta))  stop("factor$beta cannot be null");
	if(nrow(obs) == 0) stop("No observation data (i.e., nrow(obs) == 0)");

	# Prepare the feature-based prior
	xb = reg.predict(model=param$b, x=feature$x_obs, algo=param$reg.algo);
	x_src.g0 = NULL; x_dst.d0 = NULL; x_ctx.h0 = NULL; x_src.G = NULL; x_dst.D = NULL; x_ctx.H = NULL;

	if(!is.null(param$g0)) x_src.g0 = reg.predict(model=param$g0, x=feature$x_src, algo=param$reg.algo, ncol=size$nSrcContexts);
	if(!is.null(param$d0)) x_dst.d0 = reg.predict(model=param$d0, x=feature$x_dst, algo=param$reg.algo, ncol=size$nDstContexts);
	if(!is.null(param$h0)) x_ctx.h0 = reg.predict(model=param$h0, x=feature$x_ctx, algo=param$reg.algo);
	if(!is.null(param$G))   x_src.G = reg.predict(model=param$G,  x=feature$x_src, algo=param$reg.algo, ncol=size$nFactors);
	if(!is.null(param$D))   x_dst.D = reg.predict(model=param$D,  x=feature$x_dst, algo=param$reg.algo, ncol=size$nFactors);
	if(!is.null(param$H))   x_ctx.H = reg.predict(model=param$H,  x=feature$x_ctx, algo=param$reg.algo, ncol=size$nFactors);

	# Allocate space for the output
	if(is.null(factor$fScore)) factor$fScore = rep(double(1),size$nObs);
	if(is.null(factor$alpha_global) && size$nSrcContexts > 1) factor$alpha_global = rep(double(1),size$nSrcNodes);
	if(is.null(factor$beta_global)  && size$nDstContexts > 1) factor$beta_global  = rep(double(1),size$nDstNodes);
	sampvar = list();
	sampvar$fScore = rep(double(1), size$nObs);
	sampvar$alpha = array(double(1),dim=c(size$nSrcNodes, size$nSrcContexts));
	sampvar$beta  = array(double(1),dim=c(size$nDstNodes, size$nDstContexts));
	if(!is.null(factor$gamma)) sampvar$gamma = rep(double(1),size$nEdgeContexts);
	if(!is.null(factor$u)) sampvar$u = array(double(1),dim=c(size$nSrcNodes, size$nFactors));
	if(!is.null(factor$v)) sampvar$v = array(double(1),dim=c(size$nDstNodes, size$nFactors));
	if(!is.null(factor$w)) sampvar$w = array(double(1),dim=c(size$nEdgeContexts, size$nFactors));
	if(size$nSrcContexts > 1) sampvar$alpha_global = rep(double(1),size$nSrcNodes);
	if(size$nDstContexts > 1) sampvar$beta_global  = rep(double(1),size$nDstNodes);
	sampcov = list();
	if(size$nSrcContexts > 1) sampcov$alpha = array(double(1),dim=c(size$nSrcNodes,size$nSrcContexts));
	if(size$nDstContexts > 1) sampcov$beta  = array(double(1),dim=c(size$nDstNodes,size$nDstContexts));

	test.fScore = list();
	if(nTestObs > 0){
		test.fScore$mean = rep(double(1), nTestObs);
		test.fScore$var  = rep(double(1), nTestObs);
	}

	obs.y = drop(obs$y - xb);

    #  dim = {0:nObs, 1:nAlpha, 2:nBeta, 3:nrowU, 4:nrowV, 5:nAlphaContexts, 6:nBetaContexts,
    #         7:nEdgeContexts, 8:nFactors,
    #         9:nVar_y, 10:nVar_alpha, 11:nVar_beta, 12:nVar_u, 13:nVar_v, 14:nVar_w
    #        15:nVar_alpha_global, 16:nVar_beta_global, 17:nAlpha_prior, 18:nBeta_prior, 19:nGamma, 20:nVar_gamma}
    #		0               1                      2
	dim = c(nObs=size$nObs, nAlpha=size$nSrcNodes, nBeta=size$nDstNodes,
	#       3                 4
	        nrowU=size$nrowU, nrowV=size$nDstNodes,
	#       5                                 6
	        nAlphaContexts=size$nSrcContexts, nBetaContexts=size$nDstContexts,
    #       7                                 8
	        nEdgeContexts=size$nEdgeContexts, nFactors=size$nFactors,
	#       9                   10                          11
			nVar_y=size$nVar_y, nVar_alpha=size$nVar_alpha, nVar_beta=size$nVar_beta,
	#       12                  13                  14
			nVar_u=size$nVar_u, nVar_v=size$nVar_v, nVar_w=size$nVar_w,
	#       15                                        16
			nVar_alpha_global=size$nVar_alpha_global, nVar_beta_global=size$nVar_beta_global,
	#       17                                         18
			nAlpha_prior=as.integer(length(x_src.g0)), nBeta_prior=as.integer(length(x_dst.d0)),
	#       19                                       20
			nGamma=as.integer(length(factor$gamma)), nVar_gamma=as.integer(length(param$var_gamma)),
	#       21
			nTestObs=nTestObs);

	check_type_size(factor$alpha, "double", dim[c("nAlpha","nAlphaContexts")]);
	check_type_size(factor$beta,  "double", dim[c("nBeta", "nBetaContexts")]);
	check_type_size(factor$gamma, "double", dim[c("nGamma")]);
	check_type_size(factor$u, "double", dim[c("nrowU", "nFactors")]);
	check_type_size(factor$v, "double", dim[c("nrowV", "nFactors")]);
	check_type_size(factor$w, "double", dim[c("nEdgeContexts", "nFactors")]);
	check_type_size(factor$alpha_global, "double", dim[c("nAlpha")], isNullOK=(size$nSrcContexts==1));
	check_type_size(factor$beta_global,  "double", dim[c("nBeta")],  isNullOK=(size$nDstContexts==1));
	check_type_size(factor$fScore,  "double", dim[c("nObs")]);
	check_type_size(sampvar$fScore, "double", dim[c("nObs")]);
	check_type_size(sampvar$alpha,        "double", dim[c("nAlpha","nAlphaContexts")]);
	check_type_size(sampvar$alpha_global, "double", dim[c("nAlpha")],                  isNullOK=(size$nSrcContexts==1));
	check_type_size(sampcov$alpha,        "double", dim[c("nAlpha","nAlphaContexts")], isNullOK=(size$nSrcContexts==1));
	check_type_size(sampvar$beta,         "double", dim[c("nBeta", "nBetaContexts")]);
	check_type_size(sampvar$beta_global,  "double", dim[c("nBeta")],                   isNullOK=(size$nDstContexts==1));
	check_type_size(sampcov$beta,         "double", dim[c("nBeta", "nBetaContexts")],  isNullOK=(size$nDstContexts==1));
	check_type_size(sampvar$gamma,        "double", dim[c("nGamma")]);
	check_type_size(sampvar$u,            "double", dim[c("nrowU", "nFactors")]);
	check_type_size(sampvar$v,            "double", dim[c("nrowV", "nFactors")]);
	check_type_size(sampvar$w,            "double", dim[c("nEdgeContexts", "nFactors")]);
	check_type_size(obs$src.id, "int", dim["nObs"]);
	check_type_size(obs$dst.id, "int", dim["nObs"]);
	check_type_size(obs$src.context,  "int", dim["nObs"], isNullOK=(size$nSrcContexts==1));
	check_type_size(obs$dst.context,  "int", dim["nObs"], isNullOK=(size$nDstContexts==1));
	check_type_size(obs$edge.context, "int", dim["nObs"], isNullOK=(size$nEdgeContexts==0));
	check_type_size(param[["q"]], "double", dim[c("nAlphaContexts")], isNullOK=(dim["nAlphaContexts"]==1));
	check_type_size(param[["r"]], "double", dim[c("nBetaContexts")],  isNullOK=(dim["nBetaContexts"]==1));
	check_type_size(obs.y,   "double", dim["nObs"]);
	check_type_size(x_src.g0, "double", dim[c("nAlpha", "nAlphaContexts")]);
	check_type_size(x_dst.d0, "double", dim[c("nBeta",  "nBetaContexts")]);
	check_type_size(x_ctx.h0, "double", dim[c("nGamma")]);
	check_type_size(x_src.G,  "double", dim[c("nrowU",  "nFactors")]);
	check_type_size(x_dst.D,  "double", dim[c("nrowV",  "nFactors")]);
	check_type_size(x_ctx.H,  "double", dim[c("nEdgeContexts", "nFactors")]);
	check_type_size(param$var_y,     "double", dim["nVar_y"]);
	check_type_size(param$var_alpha, "double", dim["nVar_alpha"]);
	check_type_size(param$var_beta,  "double", dim["nVar_beta"]);
	check_type_size(param$var_gamma, "double", dim["nVar_gamma"]);
	check_type_size(param$var_u,     "double", dim["nVar_u"]);
	check_type_size(param$var_v,     "double", dim["nVar_v"]);
	check_type_size(param$var_w,     "double", dim["nVar_w"]);
	check_type_size(param$var_alpha_global, "double", dim["nVar_alpha_global"]);
	check_type_size(param$var_beta_global,  "double", dim["nVar_beta_global"]);
	check_type_size(test.obs$src.id, "int", dim["nTestObs"], isNullOK=(nTestObs == 0));
	check_type_size(test.obs$dst.id, "int", dim["nTestObs"], isNullOK=(nTestObs == 0));
	check_type_size(test.obs$src.context,  "int", dim["nTestObs"], isNullOK=(nTestObs == 0 || size$nSrcContexts==1));
	check_type_size(test.obs$dst.context,  "int", dim["nTestObs"], isNullOK=(nTestObs == 0 || size$nDstContexts==1));
	check_type_size(test.obs$edge.context, "int", dim["nTestObs"], isNullOK=(nTestObs == 0 || size$nEdgeContexts==0));

	var_w = param$var_w;
	if(all(var_w == 0)) dim["nVar_w"] = as.integer(0);
	if(!is.integer(dim)) stop("!is.integer(dim)");

	ans = .Call("MCEM_EStep_multicontext_Call",
    	# INPUT (initial factor values) & OUTPUT (Monte Carlo mean of factor values)
    	factor$alpha, factor$beta, factor$gamma, factor$u, factor$v, factor$w,
        # OUTPUT
    	factor$alpha_global, factor$beta_global, factor$fScore, sampvar$fScore,
    	sampvar$alpha, sampvar$alpha_global, sampcov$alpha,
    	sampvar$beta,  sampvar$beta_global,  sampcov$beta,
		sampvar$gamma, sampvar$u, sampvar$v, sampvar$w,
		test.fScore$mean, test.fScore$var,
    	# INPUT
    	as.integer(nSamples), as.integer(nBurnIn),
    	obs$src.id, obs$dst.id, obs$src.context, obs$dst.context, obs$edge.context,
    	param[["q"]], param[["r"]],
    	obs.y, x_src.g0, x_dst.d0, x_ctx.h0, x_src.G, x_dst.D, x_ctx.H,
    	param$var_y, param$var_alpha, param$var_alpha_global, param$var_beta, param$var_beta_global,
		param$var_gamma, param$var_u, param$var_v, var_w,
		test.obs$src.id, test.obs$dst.id, test.obs$src.context, test.obs$dst.context, test.obs$edge.context,
		dim, as.integer(length(dim)),
    	# OTHER
    	as.integer(debug), as.integer(verbose)
	);

	output = list(
		mean=factor, var=sampvar, cov=sampcov
	);
	if(nTestObs > 0) output$test.fScore = test.fScore;
	return(output);
}


###
### Simple random init
###	  alpha, beta, gamma, u, v, w are generated by N(0,sd)
###   b: set to global bias
###   g0, d0, h0, G, D, H are all 0
###
init.simple.random <- function(
	data.train=NULL, obs=NULL, feature=NULL, nFactors, has.u, has.gamma,
	nLocalFactors=NULL, is.logistic=FALSE,
	var_alpha=1, var_beta=1, var_gamma=1,
	var_v=1, var_u=1, var_w=1, var_y=1,
	relative.to.var_y=FALSE, var_alpha_global=1, var_beta_global=1
){
	# setup obs, feature, test.obs, test.feature
	if(!is.null(data.train)){
		if(!all(c("obs", "feature") %in% names(data.train))) stop("Please check input parameter 'data.train' when calling function fit.multicontext or run.multicontext: data.train$obs and data.train$feature cannot be NULL");
		if(!is.null(obs)) stop("When calling function fit.multicontext or run.multicontext, if you already specified 'data.train', then you should set 'obs=NULL'");
		if(!is.null(feature)) stop("When calling function fit.multicontext or run.multicontext, if you already specified 'data.train', then you should set 'feature=NULL'");
		obs=data.train$obs;
		feature=data.train$feature;
		data.train$obs = NULL;
		data.train$feature = NULL;
	}else{
		if(is.null(obs) || is.null(feature)) stop("Please specify input parameter 'data.train' when calling function fit.multicontext or run.multicontext");
	}

	nObs = nrow(obs); nSrcNodes = nrow(feature$x_src); nDstNodes = nrow(feature$x_dst);
	nSrcFeatures  = ncol(feature$x_src);
	nDstFeatures  = ncol(feature$x_dst);
	nCtxFeatures  = if(is.null(feature$x_ctx))    0 else ncol(feature$x_ctx);
	nSrcContexts  = if(is.null(obs$src.context))  1 else max(obs$src.context);
	nDstContexts  = if(is.null(obs$dst.context))  1 else max(obs$dst.context);
	nEdgeContexts = if(is.null(obs$edge.context)) 0 else max(obs$edge.context);

	if(!is.null(nLocalFactors) && nLocalFactors == 0) nLocalFactors = NULL;

	if(!has.u && nSrcNodes != nDstNodes) stop("When has.u = FALSE, the number of source nodes should be the same as the number of destination nodes (i.e., nrow(feature$x_src) == nrow(feature$x_dst)).");
	if(has.gamma && nEdgeContexts == 0) stop("When has.gamma = TRUE, obs$edge.context cannot be NULL");
	if(!is.null(nLocalFactors) && nFactors != nLocalFactors*nEdgeContexts) stop("nFactors != nLocalFactors*nEdgeContexts");

	b = rep(0, ncol(feature$x_obs));
	if(is.logistic){
		if(is.null(var_y)) var_y = 1;
		y = obs$y;
		y.values = unique(y);
		if(length(y.values) != 2) stop("When is.logistic = TRUE, the response should be binary: obs$y belongs to either {0, 1} or {-1, 1}.");
		y.values = sort(y.values);
		if(all(y.values == c(-1, 1))) y[y == -1] = 0;
		if(any(y.values != c( 0, 1))) stop("When is.logistic = TRUE, the response should be binary: obs$y belongs to either {0, 1} or {-1, 1}.");
		y.mean = mean(y);
		bias = log(y.mean / (1-y.mean));
		if(all(feature$x_obs[,1] == 1)){
			b[1] = bias;
		}else if(all(feature$x_obs[,ncol(feature$x_obs)] == 1)){
			b[ncol(feature$x_obs)] = bias;
		}
	}else{
		bias = mean(obs$y);
		if(all(feature$x_obs[,1] == 1)){
			b[1] = bias;
		}else if(all(feature$x_obs[,ncol(feature$x_obs)] == 1)){
			b[ncol(feature$x_obs)] = bias;
		}
		if(is.null(var_y)) var_y = var(obs$y);
	}
	if(relative.to.var_y){
		var_alpha=var_alpha*var_y; var_beta=var_beta*var_y;
		var_alpha_global=var_alpha_global*var_y; var_beta_global=var_beta_global*var_y;
		var_gamma=var_gamma*var_y; var_v=var_v*var_y; var_u=var_u*var_y; var_w=var_w*var_y;
	}

	factor = list(
		alpha = matrix(rnorm(nSrcNodes*nSrcContexts, sd=sqrt(var_alpha)), nrow=nSrcNodes),
		beta  = matrix(rnorm(nDstNodes*nDstContexts, sd=sqrt(var_beta )), nrow=nDstNodes)
	);
	if(nSrcContexts > 1) factor$alpha_global = rnorm(nSrcNodes, sd=sqrt(var_alpha_global));
	if(nDstContexts > 1) factor$beta_global  = rnorm(nDstNodes, sd=sqrt(var_beta_global));
	if(has.gamma) factor$gamma = rnorm(nEdgeContexts, sd=sqrt(var_gamma));
	if(nFactors > 0){
		if(has.u) factor$u = matrix(rnorm(nSrcNodes*nFactors, sd=sqrt(var_u)), nrow=nSrcNodes);
		factor$v = matrix(rnorm(nDstNodes*nFactors, sd=sqrt(var_v)), nrow=nDstNodes);
		if(nEdgeContexts > 0){
			if(is.null(nLocalFactors)){
				factor$w = matrix(rnorm(nEdgeContexts*nFactors, sd=sqrt(var_w)), nrow=nEdgeContexts);
			}else{
				factor$w = array(0.0, dim=c(nEdgeContexts,nFactors));
				for(k in 1:nEdgeContexts) factor$w[k, select.factor.indices(k,nEdgeContexts,nFactors) ] = 1;
			}
		}
	}

	if(length(var_alpha) == 1 && nSrcContexts > 1) var_alpha = rep(var_alpha, nSrcContexts);
	if(length(var_beta)  == 1 && nDstContexts > 1) var_beta  = rep(var_beta,  nDstContexts);

	param = list(
		b=b,
		g0=matrix(0.0,nrow=nSrcFeatures,ncol=nSrcContexts),
		d0=matrix(0.0,nrow=nDstFeatures,ncol=nDstContexts)
	);
	if(has.gamma) param$h0 = rep(0.0, nCtxFeatures);
	if(nFactors > 0){
		if(has.u) param$G = matrix(0.0,nrow=nSrcFeatures,ncol=nFactors);
		param$D = matrix(0.0,nrow=nDstFeatures,ncol=nFactors);
		if(nEdgeContexts > 0) param$H = matrix(0.0,nrow=nCtxFeatures,ncol=nFactors);
	}
	if(nSrcContexts > 1) param[["q"]] = rep(1.0, nSrcContexts);
	if(nDstContexts > 1) param[["r"]] = rep(1.0, nDstContexts);

	param$var_y = var_y;  param$var_alpha = var_alpha;  param$var_beta = var_beta;
	if(nSrcContexts > 1) param$var_alpha_global = var_alpha_global;
	if(nDstContexts > 1) param$var_beta_global = var_beta_global;
	if(has.gamma) param$var_gamma = var_gamma;
	if(nFactors > 0){
		n = if(is.null(nLocalFactors)) 1 else nFactors;
		if(has.u) param$var_u = rep(var_u, n);
		param$var_v = rep(var_v, n);
		if(nEdgeContexts > 0) param$var_w = if(is.null(nLocalFactors)) var_w else 0;
	}

	if(is.logistic) param$xi = rep(1.0, nrow(obs));
	if(!is.null(nLocalFactors)) param$nLocalFactors = nLocalFactors;

	ans = list(factor=factor, param=param);
}
