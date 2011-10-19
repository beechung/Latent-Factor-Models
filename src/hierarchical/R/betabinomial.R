### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

### Fit a Beta-Binomial DAG model
###
### Input: obs.data = data.frame(pNodeID, tID, nPositives, nTrials)
###					  pNodeID: ID of a population node
###					  tID: ID of a treatment
###        pop.dag  = data.frame(nodeID, parentID)
###		   tre.dag  = data.frame(nodeID, parentID)
###		   tre.map  = data.frame(tID, tNodeID)
###
###		   var_logitMean and var_logSize are arrays of values to be tried.
###		   (var_logitMean[i], var_logSize[i]) is the ith possible pair.
###
### NODE ID ORDERING:
###   * Root should have ID = 0.
###   * Population nodes and treatment nodes are labeled separately; i.e., both start from 0.
###   * For any given node, every one of its ancestors must have ID < the ID of the given node.
###
### INPUT OBSERVATION DATA ORDERING:
###   * Observation data table: obs.data${pNodeID, tID, nPositives, nTrials}
###   * This table should include all pNodeIDs.  For a non-leaf pNodeID, nPositives and nTrials
###     are sum(nPositives) and sum(nTrials) over the leaves under it.
###   * This table must be ordered by pNodeID
###   * cubeAggregate.pDAG(...) can be used to create this table
###
### Output:
###	  param[pnode+1,tnode+1,] = (mean, size) of the Beta distribution for (pnode, tnode)
###   status[pnode+1,tnode+1] = 1: Moment estimate
###								2: Posterior mode
###								0: Low support (no estimate)
###							   -1: Invalid (the node was not reached)
###							   -2: Should compute posterior mode, but only the moment estimate is ready
###							   -3: Posterior mode computation failed
### Options:
###		option=0: moment estimates only
###		option=1: postierior mode (use the likelihood criterion for selecting immediate parents)
###		option=2: postierior mode (use fixed priors for all nodes, specified by root_priorMean and root_priorSize)
###
fit.BetaBinomialDAG <- function(
		obs.data, # observation data: (pNodeID, tID, nPositives, nTrials)
		pop.dag,  # population DAG: (nodeID, parentID)
		tre.dag,  # treatment  DAG: (nodeID, parentID)
		tre.map,  # treatment ID to node mapping: (tID, tNodeID)
		obs.tune=NULL,  # tuning data: (pNodeID, tID, nPositives, nTrials)
		obs.test=NULL,  #   test data: (pNodeID, tID, nPositives, nTrials)
		option=1,       # -1: naive estimates;  0: moment estimates only;  1: also estimate postierior modes
		criterion=0,    # 0: invalid;  1: entire loglikelihood;  2: frontier loglikelihood;  3: leaf loglikelihood
		ME_nRefines=0,  # number of refinements for the moment estimator
		root_priorMean=NULL, # prior for the mean of Beta at root
		root_priorSize=NULL, # prior for the size of Beta at root
		var_logitMean=NULL,  # variance of logit(mu) from its parent node
		var_logSize=NULL,    # variance of log(size) from its parent node
		ME_threshold_numObs=NULL,    # Use moment estimates for notes with
		ME_threshold_numTrials=NULL, # #obs >= ME_threshold_numObs and avg(#trial) >= ME_threshold_numTrials
		nTrials_ForSmall=100, nObs_ForLarge=1e6, epsilon=1e-10, stepSize=1, 
		minMean=1e-10, maxMean=1-1e-10, minSize=1e-10, maxSize=1e10,
		minNumObs=3, minTotalTrials=10, minTotalPostives=0,
		maxIter1=20, nLnsrchStep1=0, maxIter2=100, nLnsrchStep2=10,
		out.level=0,   # out.level=1: save the current best model in <out.file>.model
		out.file=NULL, # out.level=2: save all the tuning models in <out.file>.model.<i>
		out.overwrite=FALSE,
		verbose=0, debug=0
){
	if(option %in% c(0, -1)){
		root_priorMean=-1; root_priorSize=-1; var_logitMean=-1; var_logSize=-1;
		ME_threshold_numObs=0; ME_threshold_numTrials=0;
	}
	if(length(var_logitMean) != length(var_logSize)) stop("var_logitMean and var_logSize must have the same length");
	if(is.null(root_priorMean)) stop("Please specify root_priorMean");
	if(is.null(root_priorSize)) stop("Please specify root_priorSize");
	if(is.null(var_logitMean)) stop("Please specify var_logitMean");
	if(is.null(var_logSize)) stop("Please specify var_logSize");
	if(is.null(ME_threshold_numObs)) stop("Please specify ME_threshold_numObs");
	if(is.null(ME_threshold_numTrials)) stop("Please specify ME_threshold_numTrials");
	if(length(var_logitMean) > 1 && is.null(obs.tune)) stop("When length(var_logitMean) > 1, you must specify obs.tune");
	if(length(var_logitMean) > 1 && !(criterion %in% 1:3)) stop("When length(var_logitMean) > 1, you must specify criterion (1:entire loglikelihood, 2:frontier loglikelihood, 3:leaf loglikelihood)");
	if(!is.null(out.file) && out.level <=0) stop("When you specify out.file, you must specify out.level > 0");
	
	if(out.level > 0){
		if(is.null(out.file)) stop("out.level > 0 && is.null(out.file)");
		file.model   = paste(out.file, ".model", sep="");
		file.summary = paste(out.file, ".summary", sep="");
		if(!out.overwrite && (file.exists(file.model) || file.exists(file.summary))) stop("Output file exists");
	}
	best.model = NULL;
	best.loglik = -Inf;
	smry = data.frame(var_logitMean=var_logitMean, var_logSize=var_logSize, time.used=NA, nLowSupp=NA, nFailed=NA);
	
	if(!is.null(obs.tune)){
		smry$tune.loglik.entire = NA;
		smry$tune.loglik.frontier = NA;
		smry$tune.loglik.leaf = NA;
	}
	if(!is.null(obs.test)){
		smry$test.loglik.entire = NA;
		smry$test.loglik.frontier = NA;
		smry$test.loglik.leaf = NA;
	}
	
	for(i in 1:length(var_logitMean)){
		
		time.begin = proc.time();
		model = fit.BetaBinomialDAG.C(
			obs.data=obs.data, pop.dag=pop.dag, tre.dag=tre.dag, tre.map=tre.map,
			option=option, ME_nRefines=ME_nRefines, root_priorMean=root_priorMean, root_priorSize=root_priorSize, 
			var_logitMean=var_logitMean[i], var_logSize=var_logSize[i],
			ME_threshold_numObs=ME_threshold_numObs, ME_threshold_numTrials=ME_threshold_numTrials,
			nTrials_ForSmall=nTrials_ForSmall, nObs_ForLarge=nObs_ForLarge, epsilon=epsilon, stepSize=stepSize, 
			minMean=minMean, maxMean=maxMean, minSize=minSize, maxSize=maxSize,
			minNumObs=minNumObs, minTotalTrials=minTotalTrials, minTotalPostives=minTotalPostives,
			maxIter1=maxIter1, nLnsrchStep1=nLnsrchStep1, maxIter2=maxIter2, nLnsrchStep2=nLnsrchStep2,
			verbose=verbose, debug=debug
		);
		time.used = proc.time() - time.begin;
		smry[i,"time.used"] = time.used[3];
		smry[i,"nLowSupp"]  = model$numLowSupportNodes;
		smry[i,"nFailed"]   = model$numFailedNodes;
		
		if(!is.null(obs.tune)){
			loglik.tune = loglik.BetaBinomialDAG(
				model=model, obs.data=obs.tune, pop.dag=pop.dag, tre.dag=tre.dag, tre.map=tre.map,
				option=1, output.nodeLoglik=TRUE,
				nTrials_ForSmall=nTrials_ForSmall, nObs_ForLarge=nObs_ForLarge, 
				minMean=minMean, maxMean=maxMean, minSize=minSize, maxSize=maxSize, avg.loglik=TRUE,
				verbose=verbose, debug=debug
			);
			smry[i,"tune.loglik.entire"]   = loglik.tune$avg.entire.loglik;
			smry[i,"tune.loglik.frontier"] = loglik.tune$avg.frontier.loglik;
			smry[i,"tune.loglik.leaf"]     = loglik.tune$avg.leaf.loglik;
		}
		if(!is.null(obs.test)){
			loglik.test = loglik.BetaBinomialDAG(
					model=model, obs.data=obs.test, pop.dag=pop.dag, tre.dag=tre.dag, tre.map=tre.map,
					option=1, output.nodeLoglik=TRUE,
					nTrials_ForSmall=nTrials_ForSmall, nObs_ForLarge=nObs_ForLarge, 
					minMean=minMean, maxMean=maxMean, minSize=minSize, maxSize=maxSize, avg.loglik=TRUE,
					verbose=verbose, debug=debug
			);
			smry[i,"test.loglik.entire"]   = loglik.test$avg.entire.loglik;
			smry[i,"test.loglik.frontier"] = loglik.test$avg.frontier.loglik;
			smry[i,"test.loglik.leaf"]     = loglik.test$avg.leaf.loglik;
		}
		
		model.changed = FALSE;
		if(criterion %in% c(1,2,3)){
			if(criterion == 1){ name = "tune.loglik.entire";}
			else if(criterion == 2){ name = "tune.loglik.frontier";}
			else if(criterion == 3){ name = "tune.loglik.leaf";}
			if(smry[i,name] > best.loglik){
				best.model = model;
				best.loglik = smry[i,name];
				model.changed = TRUE;
			}
		}else if(length(var_logitMean)==1){
			best.model = model;
			best.loglik = NA;
			model.changed = TRUE;
		}else stop("criterion=",criterion,", which must be 1, 2, 3");
		
		if(out.level > 0){
			append = (i>1);
			write.table(smry[i,], file=file.summary, row.names=FALSE, col.names=!append, append=append, quote=FALSE, sep="\t");
			if(model.changed) save(model, file=file.model);
			if(out.level >= 2) save(model, file=paste(file.model,".",i,sep=""));
		}
	}
	best.model$summary = smry;
	return(best.model);
}

sample.moment.est <- function(
	obs.data, tre.map
){
	d = merge(obs.data, tre.map);
	d$mean = d$nPositives / d$nTrials;
	d$var  = d$nPositives / d$nTrials;
	agg = table.aggregate(d, group.by=c("pNodeID","tNodeID"),attr=c("mean", "var"),fn=c(mean, var), use.joined.column=FALSE);
	agg$size = size.beta(mu=agg$mean, var=agg$var);
	param = array(NA, dim=c(max(d$pNodeID)+1, max(d$tNodeID)+1, 2));
	param[cbind(agg$pNodeID+1,agg$tNodeID+1,1)] = agg$mean;
	param[cbind(agg$pNodeID+1,agg$tNodeID+1,2)] = agg$size;
	return(param);
}

fit.BetaBinomialDAG.R <- function(
	obs.data.full, # observation data: (pNodeID, tNodeID, nPositives, nTrials)
	pop.dag,  # population DAG: (nodeID, parentID)
	tre.dag,  # treatment  DAG: (nodeID, parentID)
	option=1, # 0: moment estimates only;  1: also estimate postierior modes
	ME_nRefines=0,  # number of refinements for the moment estimator
	root_priorMean, # prior for the mean of Beta at root
	root_priorSize, # prior for the size of Beta at root
	var_logitMean,  # variance of logit(mu) from its parent node
	var_logSize,    # variance of log(size) from its parent node
	ME_threshold_numObs,    # Use moment estimates for notes with
	ME_threshold_numTrials, # #obs >= ME_threshold_numObs and avg(#trial) >= ME_threshold_numTrials
	nTrials_ForSmall=0, nObs_ForLarge=1e6, epsilon=1e-10, stepSize=1, 
	maxIter1=20, nLnsrchStep1=0, maxIter2=100, nLnsrchStep2=10, 
	verbose=0, debug=0,
	use.C=FALSE
){
	maxPNode = max(pop.dag$nodeID);
	maxTNode = max(tre.dag$nodeID);
	param  = array(0.0, dim=c(maxPNode+1, maxTNode+1, 2));
	for(pnode in 0:maxPNode){
		for(tnode in 0:maxTNode){
			obs = obs.data.full[obs.data.full$pNodeID==pnode & obs.data.full$tNodeID==tnode,];
			x = obs$nPositives;
			n = obs$nTrials;
			# moment estimate
			me = betabinom.MoM(x=x, n=n, w=1, refine=ME_nRefines);
			me.mu = me[1];  
			rho=me[2]; # if(rho < 1e-20) rho = 1e-20; if(rho > 1 - 1e-20) rho = 1 - 1e-20;
			me.gamma = 1/rho - 1;
			if(verbose >= 5) cat("(",pnode,", ",tnode,"): ME=c(",me.mu,", ",me.gamma,")\n",sep="");
			if((nrow(obs) >= ME_threshold_numObs && mean(obs$nTrials) >= ME_threshold_numTrials) ||
			   	option == 0){
				param[pnode+1,tnode+1,] = c(me.mu, me.gamma);
				next;
			}
			# find the best prior
			if(pnode==0 && tnode==0){
				mu0 = root_priorMean;  gamma0 = root_priorSize;
			}else{
				PaP = pop.dag$parentID[pop.dag$nodeID == pnode];
				PaT = tre.dag$parentID[tre.dag$nodeID == tnode];
				if(length(PaP)+length(PaT) == 0) stop("length(PaP)+length(PaT) == 0");
				loglik = rep(NA, length(PaP)+length(PaT));
				mu_tmp = rep(NA, length(PaP)+length(PaT));
				gamma_tmp = rep(NA, length(PaP)+length(PaT));
				i = 0;
				for(p in PaP){
					theta = param[p+1,tnode+1,];
					i = i+1;
					mu_tmp[i] = theta[1];
					gamma_tmp[i] = theta[2];
					loglik[i] = loglik.BetaBin.mu.gamma(param=theta, pos=x, size=n, mu0=theta[1], gamma0=theta[2], var_mu=var_logitMean, var_gamma=var_logSize);
				}
				for(t in PaT){
					theta = param[pnode+1,t+1,];
					i = i+1;
					mu_tmp[i] = theta[1];
					gamma_tmp[i] = theta[2];
					loglik[i] = loglik.BetaBin.mu.gamma(param=theta, pos=x, size=n, mu0=theta[1], gamma0=theta[2], var_mu=var_logitMean, var_gamma=var_logSize);
				}
				index = which.max(loglik);
				mu0 = mu_tmp[index];
				gamma0 = gamma_tmp[index];
			}
			# find posterior mode
			if(use.C){
				pm = fit.BetaBin(
						pos=x, size=n, mu0=mu0, gamma0=gamma0, var_mu=var_logitMean, var_gamma=var_logSize,
						threshold=nTrials_ForSmall, nSamplesForLarge=nObs_ForLarge, epsilon=epsilon, stepSize=stepSize, 
						maxIter1=maxIter1, nLnsrchStep1=nLnsrchStep1, maxIter2=maxIter2, nLnsrchStep2=nLnsrchStep2, 
						verbose=verbose+10, debug=debug
				);
			}else{
				init.gamma = if(me.gamma == Inf || me.gamma == 0) gamma0 else me.gamma;
				temp = nlm(deriv.negLoglik.BetaBin.logitMu.logGamma, c(logit(me.mu),log(init.gamma)), pos=x, size=n, mu0=mu0, gamma0=gamma0, var_mu=var_logitMean, var_gamma=var_logSize, check.analyticals=TRUE);
				pm = c(1/(1+exp(-temp$estimate[1])), exp(temp$estimate[2]));
			}
			param[pnode+1,tnode+1,] = pm;
		}
	}
	return(param);
}

loglik.BetaBinomialDAG.R <- function(
		model,    # list(param, status)
		obs.data.full, # observation data: (pNodeID, tNodeID, nPositives, nTrials)
		pop.dag,  # population DAG: (nodeID, parentID)
		tre.dag,  # treatment  DAG: (nodeID, parentID)
		verbose=0, debug=0
){
	maxPNode = max(pop.dag$nodeID);
	maxTNode = max(tre.dag$nodeID);
	ans = array(0.0, dim=c(maxPNode+1, maxTNode+1));
	for(pnode in 0:maxPNode){
		for(tnode in 0:maxTNode){
			if(model$status[pnode+1,tnode+1] <= 0) next;
			obs = obs.data.full[obs.data.full$pNodeID==pnode & obs.data.full$tNodeID==tnode,];
			x = obs$nPositives;
			n = obs$nTrials;
			mu = model$param[pnode+1,tnode+1,1];
			gamma = model$param[pnode+1,tnode+1,2];
			
			temp = 0;
			for(i in 1:length(n)){
				temp = temp + testLoglik.BetaBin.mu.gamma(c(mu,gamma), x[i], n[i]);
			}
			
			ans[pnode+1,tnode+1] = temp / length(n);
		}
	}
	out = list(node.loglik=ans);
	return(out);
}

pbetabinom.mu.gamma <- function(x, n, mu, gamma){
	alpha = mu * gamma;  beta = gamma - alpha;
	logP=lgamma(beta+n)+lgamma(gamma)-lgamma(beta)-lgamma(gamma+n);
	
	sum=1.0;
	Tr=1.0;
	bn=-beta-n;
	
	if(x > 0){
		for(i in 0:(x-1)) {
			r=i;
			rp=(i+1);
			Tr = Tr * ((r+alpha)*(r-n))/(rp*(bn+rp));
			sum = sum+Tr;
		}
	}
	logP = logP + log(sum);
	return(exp(logP));
}

loglik.BetaBin.mu.gamma <- function(param, pos, size, mu0, gamma0, var_mu, var_gamma){
	mu = param[1];  gamma = param[2];  alpha = mu * gamma;  beta = gamma - alpha;
	x = pos; n = size;
	ans = (- (logit(mu) - logit(mu0))^2  / (2 * var_mu)
		   - (log(gamma) - log(gamma0))^2 / (2 * var_gamma)
		   + sum(
				   lgamma(x+alpha) + lgamma(n-x+beta) + lgamma(gamma)
		   		 - lgamma(alpha)   - lgamma(beta)     - lgamma(n+gamma)
		     )
	);
	return(ans);
}
testLoglik.BetaBin.mu.gamma <- function(param, pos, size){
	mu = param[1];  gamma = param[2];  alpha = mu * gamma;  beta = gamma - alpha;
	x = pos; n = size;
	ans = ( lchoose(size,pos)
			+ sum(
					lgamma(x+alpha) + lgamma(n-x+beta) + lgamma(gamma)
				  - lgamma(alpha)   - lgamma(beta)     - lgamma(n+gamma)
			)
	);
	return(ans);
}

gr.BetaBin.mu.gamma <- function(param, pos, size, mu0, gamma0, var_mu, var_gamma){
	mu = param[1];  gamma = param[2];  alpha = mu * gamma;  beta = gamma - alpha;
	x = pos; n = size;
	gr_mu = (- (logit(mu) - logit(mu0)) * (1/var_mu) * (1/mu + 1/(1-mu))
			 + gamma * sum(digamma(x+alpha) - digamma(n-x+beta) - digamma(alpha) + digamma(beta))
	);
	gr_gamma = (- (log(gamma) - log(gamma0)) / (gamma * var_gamma)
				+ sum(
						digamma(x+alpha)*mu + digamma(n-x+beta)*(1-mu) + digamma(gamma)
					  - digamma(alpha)*mu   - digamma(beta)*(1-mu)     - digamma(n+gamma)
				  )
	);
	return(c(gr_mu, gr_gamma));
}

H.BetaBin.mu.gamma <- function(param, pos, size, mu0, gamma0, var_mu, var_gamma){
	mu = param[1];  gamma = param[2];  alpha = mu * gamma;  beta = gamma - alpha;
	x = pos; n = size;
	H_mu = (- (logit(mu) - logit(mu0)) * (1/var_mu) * (-1/(mu^2) + 1/(1-mu)^2)
			- (1/var_mu) * (1/mu + 1/(1-mu))^2
			+ gamma^2 * sum(trigamma(x+alpha) + trigamma(n-x+beta) - trigamma(alpha) - trigamma(beta))
	);
	H_gamma = (- (1 - log(gamma) + log(gamma0)) / (gamma^2 * var_gamma)
			   + sum(
					   trigamma(x+alpha)*mu^2 + trigamma(n-x+beta)*(1-mu)^2 + trigamma(gamma)
			   		 - trigamma(alpha)*mu^2   - trigamma(beta)*(1-mu)^2     - trigamma(n+gamma)
			     )
	);
	H_muga = sum(
				digamma(x+alpha) - digamma(n-x+beta) - digamma(alpha) + digamma(beta)
			  + trigamma(x+alpha)*alpha - trigamma(n-x+beta)*beta - trigamma(alpha)*alpha + trigamma(beta)*beta
	);
	return(cbind(c(H_mu, H_muga), c(H_muga, H_gamma)));
}

deriv.prior.BetaBin.mu.gamma <- deriv(
	y ~ -(log(mu)-log(1-mu)-log(mu0)+log(1-mu0))^2/(2*var_mu) - (log(gamma)-log(gamma0))^2/(2*var_gamma),
	c("mu", "gamma"), c("mu", "gamma", "mu0", "gamma0", "var_mu", "var_gamma"), hessian=TRUE);

deriv.singleObs.BetaBin.mu.gamma <- deriv(
	y ~ lgamma(x + mu*gamma) + lgamma(n - x + (1-mu)*gamma) + lgamma(gamma) - lgamma(mu*gamma) - lgamma((1-mu)*gamma) - lgamma(n + gamma),
	c("mu", "gamma"), c("mu", "gamma", "x", "n"), hessian=TRUE);

deriv.BetaBin.mu.gamma <- function(param, pos, size, mu0, gamma0, var_mu, var_gamma){
	prior = deriv.prior.BetaBin.mu.gamma(mu=param[1], gamma=param[2], mu0=mu0, gamma0=gamma0, var_mu=var_mu, var_gamma=var_gamma);
	obs   = deriv.singleObs.BetaBin.mu.gamma(mu=param[1], gamma=param[2], x=pos, n=size);
	out = prior[1] + sum(obs);
	attr(out, "gradient") = attr(prior,"gradient")[1,] + apply(attr(obs,"gradient"),2,sum);
	attr(out, "hessian")  = attr(prior,"hessian")[1,,] + apply(attr(obs,"hessian"),c(2,3),sum);
	return(out);
}

deriv.negLoglik.BetaBin.mu.gamma <- function(param, pos, size, mu0, gamma0, var_mu, var_gamma){
	out = deriv.BetaBin.mu.gamma(param, pos, size, mu0, gamma0, var_mu, var_gamma);
	out = - out;
	attr(out, "gradient") = -attr(out, "gradient");
	attr(out, "hessian")  = -attr(out, "hessian");
	return(out);
}

deriv.prior.BetaBin.logitMu.logGamma <- deriv(
		y ~ -(logitMu-log(mu0)+log(1-mu0))^2/(2*var_mu) - (logGamma-log(gamma0))^2/(2*var_gamma),
		c("logitMu", "logGamma"), c("logitMu", "logGamma", "mu0", "gamma0", "var_mu", "var_gamma"), hessian=TRUE);

deriv.singleObs.BetaBin.logitMu.logGamma <- deriv(
		y ~ lgamma(x + (1/(1+exp(-logitMu)))*exp(logGamma)) + lgamma(n - x + (1-(1/(1+exp(-logitMu))))*exp(logGamma)) + lgamma(exp(logGamma)) - lgamma((1/(1+exp(-logitMu)))*exp(logGamma)) - lgamma((1-(1/(1+exp(-logitMu))))*exp(logGamma)) - lgamma(n + exp(logGamma)),
		c("logitMu", "logGamma"), c("logitMu", "logGamma", "x", "n"), hessian=TRUE);

deriv.BetaBin.logitMu.logGamma <- function(param, pos, size, mu0, gamma0, var_mu, var_gamma){
	prior = deriv.prior.BetaBin.logitMu.logGamma(logitMu=param[1], logGamma=param[2], mu0=mu0, gamma0=gamma0, var_mu=var_mu, var_gamma=var_gamma);
	obs   = deriv.singleObs.BetaBin.logitMu.logGamma(logitMu=param[1], logGamma=param[2], x=pos, n=size);
	out = prior[1] + sum(obs);
	attr(out, "gradient") = attr(prior,"gradient")[1,] + apply(attr(obs,"gradient"),2,sum);
	attr(out, "hessian")  = attr(prior,"hessian")[1,,] + apply(attr(obs,"hessian"),c(2,3),sum);
	return(out);
}

deriv.negLoglik.BetaBin.logitMu.logGamma <- function(param, pos, size, mu0, gamma0, var_mu, var_gamma){
	out = deriv.BetaBin.logitMu.logGamma(param, pos, size, mu0, gamma0, var_mu, var_gamma);
	out = - out;
	attr(out, "gradient") = -attr(out, "gradient");
	attr(out, "hessian")  = -attr(out, "hessian");
	return(out);
}


genData.BetaBin <- function(nObs, nTrials, mu0, var_mu, gamma0, var_gamma, type="Poisson"){
	logit_mu  = rnorm(n=1, mean=log(mu0)-log(1-mu0), sd=sqrt(var_mu));
	log_gamma = rnorm(n=1, mean=log(gamma0), sd=sqrt(var_gamma));
	mu    = 1/(1+exp(-logit_mu));
	gamma = exp(log_gamma);
	rho = rho.beta(size=gamma);
	n = nTrials;
	if(type == "Poisson"){
		n = 2+rpois(nObs, nTrials-2);
	}else if(type != "fixed") error("Unknown type: ",type);
	x = rbetabin(nObs, size=n, prob=mu, rho=rho);
	out = list(param=c(mu=mu, gamma=gamma), obs=data.frame(positive=as.integer(x), size=as.integer(n)));
	return(out);
}

random.integer <- function(min, max){
	return(as.integer(floor(runif(1,min=min, max=max+1-1e-20))));
}

###
### When generating the groundtruth (mu, gamma) for a node,
### randomly pick a parent u, and draw 
###		 logit(mu) ~ N(logit(mu_u),  var_logitMu)
###		log(gamma) ~ N(log(gamma_u), var_logGamma)
###
genData.BetaBinDAG.simple <- function(
	pop.dag, tre.dag,
	nObs.perLeaf, nTrials, mu0, gamma0, var_logitMu, var_logGamma, type="Poisson",
	leaves.only=FALSE, verbose=0
){
	pop.dag=data.frame(nodeID=as.integer(pop.dag$nodeID), parentID=as.integer(pop.dag$parentID));
	tre.dag=data.frame(nodeID=as.integer(tre.dag$nodeID), parentID=as.integer(tre.dag$parentID));
	
	pop.internals = unique(pop.dag$parentID);
	tre.internals = unique(tre.dag$parentID);
	pop.leaves = unique(pop.dag$nodeID);
	pop.leaves = pop.leaves[!(pop.leaves %in% pop.internals)];
	tre.leaves = unique(tre.dag$nodeID);
	tre.leaves = tre.leaves[!(tre.leaves %in% tre.internals)];
	
	nPNodes = max(pop.leaves)+1;
	nTNodes = max(tre.leaves)+1;
	
	tre.map = data.frame(
		tID=as.integer(0:(length(tre.leaves)*nObs.perLeaf-1)),
		tNodeID=as.integer(rep(tre.leaves,each=nObs.perLeaf))
	);
	
	param = array(NA, dim=c(nPNodes, nTNodes, 2));
	param[1,1,] = c(mu0, gamma0);
	
	for(pNodeID in 0:(nPNodes-1)){
		pParent = pop.dag$parentID[pop.dag$nodeID == pNodeID];
		for(tNodeID in 0:(nTNodes-1)){
			if(verbose > 1) cat("  process (",pNodeID,", ",tNodeID,"):\n",sep="");
			if(pNodeID == 0 && tNodeID == 0) next;
			tParent = tre.dag$parentID[tre.dag$nodeID == tNodeID];
			if(verbose > 2) cat("    pParents: ",paste(pParent,collapse=", "),"\n",
								"    tParents: ",paste(tParent,collapse=", "),"\n",sep="");
			if(length(pParent) == 0 || (length(tParent) != 0 && runif(1)>0.5)){
				p = tParent[random.integer(min=1, max=length(tParent))];
				prior = param[pNodeID+1,p+1,];
				if(verbose > 1) cat("    pick: (",pNodeID,", ",p,") ",sep="")
			}else{
				p = pParent[random.integer(min=1, max=length(pParent))];	
				prior = param[p+1,tNodeID+1,];
				if(verbose > 1) cat("    pick: (",p,", ",tNodeID,") ",sep="")
			}
			if(verbose > 1) cat("with param: ",paste(prior,collapse=", "),"\n",sep="")
			logit_mu  = rnorm(n=1, mean=log(prior[1])-log(1-prior[1]), sd=sqrt(var_logitMu));
			log_gamma = rnorm(n=1, mean=log(prior[2]), sd=sqrt(var_logGamma));
			mu    = 1/(1+exp(-logit_mu));
			gamma = exp(log_gamma);
			param[pNodeID+1,tNodeID+1,] = c(mu, gamma);
			if(verbose > 1) cat("    mu=",mu,", gamma=",gamma,"\n",sep="");
		}
	}

	num = length(pop.leaves) * nrow(tre.map);
	obs.pNodeID    = integer(num);
	obs.tID        = integer(num); 
	obs.nPositives = integer(num);
	obs.nTrials    = integer(num);
	k = 0;
	
	if(verbose > 1) cat("Generate observations:\n");

	for(pNodeID in pop.leaves){
		for(i in 1:length(tre.leaves)){
			tNodeID = tre.leaves[i];
			tIDs  = (i-1)*nObs.perLeaf + (0:(nObs.perLeaf-1));
			mu    = param[pNodeID+1,tNodeID+1,1];
			gamma = param[pNodeID+1,tNodeID+1,2];
			rho = rho.beta(size=gamma);
			n = nTrials;
			if(type == "Poisson"){
				n = 2+rpois(nObs.perLeaf, nTrials-2);
			}else if(type != "fixed") error("Unknown type: ",type);
			x = rbetabin(nObs.perLeaf, size=n, prob=mu, rho=rho);
			
			obs.pNodeID[k+(1:nObs.perLeaf)] = as.integer(pNodeID);
			obs.tID[k+(1:nObs.perLeaf)] = as.integer(tIDs);
			obs.nPositives[k+(1:nObs.perLeaf)] = as.integer(x);
			obs.nTrials[k+(1:nObs.perLeaf)] = as.integer(n);
			
			if(verbose > 1){
				cat("  Data (",pNodeID,", ",tNodeID,"):  mu=",mu,"  gamma=",gamma,"\n",sep="");
				temp = rbind(n=n, x=x);
				print(temp);
			}
			
			k = k+nObs.perLeaf;
		}
	}
	if(k != num) stop("error");
	
	obs = data.frame(pNodeID=obs.pNodeID, tID=obs.tID, nPositives=obs.nPositives, nTrials=obs.nTrials);
	
	if(leaves.only){
		obs.data = obs[order(obs$pNodeID, obs$tID),];
	}else{
		obs.data = cubeAggregate.pDAG(obs, pop.dag);
	}

	out = list(
		pop.dag=pop.dag,
		tre.dag=tre.dag,
		tre.map=tre.map,
		obs.data=obs.data,
		param=param
	);
	return(out);
}

betabinom.Naive <- function(x, n){
	p = x/n;
	pq = p * (1-p);
	mu = mean(p);
	rho = var(p) / mu*(1-mu);
	return(c(mu=mu, rho=rho));
}

betabinom.PM <- function(x, n, mu0, gamma0, var_logitMu, var_logGamma){
	ans = fit.BetaBin(pos=x, size=n, mu0=mu0, gamma0=gamma0, var_mu=var_logitMu, var_gamma=var_logGamma);
	return(c(mu=ans[1], rho=1/(1+ans[2])));
}

betabinom.MLE <- function(x, n, epsilon=1e-10, trace=FALSE){
	fit = vglm(cbind(x,n-x) ~ 1, betabinomial, trace=trace, epsilon=epsilon);
	return(Coef(fit));
}

betabinom.MoM <- function(x, n, w=1, refine=0){
	k = length(x);
	if(length(w) == 1) w = rep(w, k);
	w_total = sum(w);
	p_i   = x/n;
	p_hat = sum(w * p_i) / w_total;
	S     = sum( w * (p_i - p_hat)^2 ) * (k-1) / k;
	
	pq = p_hat * (1 - p_hat);
	one_minus_w_i_over_w = 1 - (w / w_total);
	part1 = sum( (w/n) * one_minus_w_i_over_w );
	part2 = sum( w * one_minus_w_i_over_w );
	rho = (S - pq * part1) / (pq * (part2 - part1));
	
	if(rho < 0) rho = 0;
	
	if(refine <= 0){
		return(c(mu=p_hat, rho=rho));
	}else{
		return(betabinom.MoM(x,n,w=n/(1+rho*(n-1)),refine=refine-1));
	}
}


var.beta <- function(alpha=NULL, beta=NULL, mu=NULL, rho=NULL, size=NULL){
	if(!is.null(alpha)){
		if(is.null(beta)) stop("Please specify beta");
		if(!is.null(mu) || !is.null(rho) || !is.null(size)) stop("Please do not specify parameters more than necessary");
		mu = alpha / (alpha + beta);
		rho = 1/(alpha + beta + 1);
	}else if(!is.null(rho)){
		if(is.null(mu)) stop("Please specify mu");
		if(!is.null(alpha) || !is.null(beta) || !is.null(size)) stop("Please do not specify parameters more than necessary");
	}else if(!is.null(size)){
		if(is.null(mu)) stop("Please specify mu");
		if(!is.null(alpha) || !is.null(beta) || !is.null(rho)) stop("Please do not specify parameters more than necessary");
		rho = 1/(size + 1);
	}else stop("Please specify enough parameters!");
	return(mu*(1-mu)*rho);
}

rho.beta <- function(alpha=NULL, beta=NULL, size=NULL, mu=NULL, var=NULL){
	if(!is.null(size)){
		if(!is.null(alpha) || !is.null(beta) || !is.null(mu) || !is.null(var)) stop("Please do not specify parameters more than necessary");
		return(1 / (size+1));
	}else if(!is.null(alpha)){
		if(is.null(beta)) stop("Please specify beta");
		if(!is.null(size) || !is.null(mu) || !is.null(var)) stop("Please do not specify parameters more than necessary");
		return(1 / (alpha+beta+1));
	}else if(!is.null(var)){
		if(is.null(mu)) stop("Please specify mu");
		if(!is.null(alpha) || !is.null(beta) || !is.null(size)) stop("Please do not specify parameters more than necessary");
		ans = var / (mu * (1-mu));
		ans[var==0] = 0;
		return(ans);
	}else stop("Please specify enough parameters!");
}

cv.beta <- function(alpha=NULL, beta=NULL, size=NULL, mu=NULL, var=NULL, rho=NULL){
	if(!is.null(alpha)){
		if(is.null(beta)) stop("Please specify beta");
		if(!is.null(mu) || !is.null(rho) || !is.null(size) || !is.null(var)) stop("Please do not specify parameters more than necessary");
		mu = alpha / (alpha + beta);
		var = var.beta(alpha=alpha, beta=beta);
	}else if(!is.null(size)){
		if(is.null(mu)) stop("Please specify mu");
		if(!is.null(alpha) || !is.null(beta) || !is.null(rho) || !is.null(var)) stop("Please do not specify parameters more than necessary");
		var = var.beta(mu=mu, size=size);
	}else if(!is.null(rho)){
		if(is.null(mu)) stop("Please specify mu");
		if(!is.null(alpha) || !is.null(beta) || !is.null(size) || !is.null(var)) stop("Please do not specify parameters more than necessary");
		var = var.beta(mu=mu, rho=rho);
	}else if(!is.null(var)){
		if(is.null(mu)) stop("Please specify mu");
		if(!is.null(alpha) || !is.null(beta) || !is.null(size) || !is.null(rho)) stop("Please do not specify parameters more than necessary");
	}else stop("Please specify enough parameters!");
	ans = sqrt(var)/mu;
	ans[var == 0] = 0;
	return(ans);
}

size.beta <- function(alpha=NULL, beta=NULL, rho=NULL, mu=NULL, var=NULL){
	if(!is.null(rho)){
		if(!is.null(alpha) || !is.null(beta) || !is.null(mu) || !is.null(var)) stop("Please do not specify parameters more than necessary");
		return(1/rho - 1);
	}else if(!is.null(alpha)){
		if(is.null(beta)) stop("Please specify beta");
		if(!is.null(rho) || !is.null(mu) || !is.null(var)) stop("Please do not specify parameters more than necessary");
		return(alpha+beta);
	}else if(!is.null(var)){
		if(is.null(mu)) stop("Please specify mu");
		if(!is.null(alpha) || !is.null(beta) || !is.null(rho)) stop("Please do not specify parameters more than necessary");
		ans = (mu * (1-mu))/var - 1;
		ans[var==0] = Inf;
		return(ans);
	}else stop("Please specify enough parameters!");
}

###############################################################################
#      BACKUP
###############################################################################

### THIS FUNCTION IS INCORRECT
### When generating the groundtruth (mu, gamma) for a node,
### randomly pick a parent u, and draw 
###		 logit(mu) ~ N(logit(mu_u),  var_logitMu)
###		log(gamma) ~ N(log(gamma_u), var_logGamma)
###
genData.BetaBinDAG.simple.old <- function(
		pop.dag, tre.dag,
		nObs.perLeaf, nTrials, mu0, gamma0, var_logitMu, var_logGamma, type="Poisson",
		leaves.only=FALSE, verbose=0
){
	pop.internals = unique(pop.dag$parentID);
	tre.internals = unique(tre.dag$parentID);
	pop.leaves = unique(pop.dag$nodeID);
	pop.leaves = pop.leaves[!(pop.leaves %in% pop.internals)];
	tre.leaves = unique(tre.dag$nodeID);
	tre.leaves = tre.leaves[!(tre.leaves %in% tre.internals)];
	
	nPNodes = max(pop.leaves)+1;
	nTNodes = max(tre.leaves)+1;
	
	tre.map = data.frame(
			tID=as.integer(0:(length(tre.leaves)*nObs.perLeaf-1)),
			tNodeID=as.integer(rep(tre.leaves,each=nObs.perLeaf))
	);
	
	param = array(NA, dim=c(nPNodes, nTNodes, 2));
	param[1,1,] = c(mu0, gamma0);
	
	for(pNodeID in 0:(nPNodes-1)){
		pParent = pop.dag$parentID[pop.dag$nodeID == pNodeID];
		for(tNodeID in 0:(nTNodes-1)){
			if(verbose > 1) cat("  process (",pNodeID,", ",tNodeID,"):\n",sep="");
			if(pNodeID == 0 && tNodeID == 0) next;
			tParent = tre.dag$parentID[tre.dag$nodeID == tNodeID];
			if(verbose > 2) cat("    pParents: ",paste(pParent,collapse=", "),"\n",
						"    tParents: ",paste(tParent,collapse=", "),"\n",sep="");
			if(length(pParent) == 0 || (length(tParent) != 0 && runif(1)>0.5)){
				p = tParent[random.integer(min=1, max=length(tParent))];
				prior = param[pNodeID+1,p+1,];
				if(verbose > 1) cat("    pick: (",pNodeID,", ",p,") ",sep="")
			}else{
				p = pParent[random.integer(min=1, max=length(pParent))];	
				prior = param[p+1,tNodeID+1,];
				if(verbose > 1) cat("    pick: (",p,", ",tNodeID,") ",sep="")
			}
			if(verbose > 1) cat("with param: ",paste(prior,collapse=", "),"\n",sep="")
			logit_mu  = rnorm(n=1, mean=log(prior[1])-log(1-prior[1]), sd=sqrt(var_logitMu));
			log_gamma = rnorm(n=1, mean=log(prior[2]), sd=sqrt(var_logGamma));
			mu    = 1/(1+exp(-logit_mu));
			gamma = exp(log_gamma);
			param[pNodeID+1,tNodeID+1,] = c(mu, gamma);
			if(verbose > 1) cat("    mu=",mu,", gamma=",gamma,"\n",sep="");
		}
	}
	
	data.nPositives = matrix(integer(1), nrow=nrow(tre.map), ncol=nPNodes);
	data.nTrials    = matrix(integer(1), nrow=nrow(tre.map), ncol=nPNodes);
	
	if(verbose > 1) cat("Generate observations:\n");
	
	for(pNodeID in pop.leaves){
		for(tNodeID in tre.leaves){
			mu    = param[pNodeID+1,tNodeID+1,1];
			gamma = param[pNodeID+1,tNodeID+1,2];
			rho = rho.beta(size=gamma);
			n = nTrials;
			if(type == "Poisson"){
				n = 2+rpois(nObs.perLeaf, nTrials-2);
			}else if(type != "fixed") error("Unknown type: ",type);
			x = rbetabin(nObs.perLeaf, size=n, prob=mu, rho=rho);
			data.nPositives[tre.map$tNodeID == tNodeID, pNodeID+1] = as.integer(x);
			data.nTrials[tre.map$tNodeID == tNodeID, pNodeID+1] = as.integer(n);
			if(verbose > 1){
				cat("  Data (",pNodeID,", ",tNodeID,"):  mu=",mu,"  gamma=",gamma,"\n",sep="");
				temp = rbind(n=n, x=x);
				print(temp);
			} 
		}
	}
	for(pNodeID in sort(pop.internals, decreasing=TRUE)){
		children = pop.dag$nodeID[pop.dag$parentID == pNodeID];
		if(verbose > 1) cat("  process (",pNodeID,", *):\n",
					"    children: ",paste(children,collapse=", "),"\n",sep="");
		data.nPositives[,pNodeID+1] = apply(data.nPositives[,children+1,drop=FALSE], 1, sum);
		data.nTrials[,pNodeID+1] = apply(data.nTrials[,children+1,drop=FALSE], 1, sum);
		
		if(verbose > 1){
			temp = rbind(n=data.nTrials[,pNodeID+1], x=data.nPositives[,pNodeID+1]);
			print(temp);
		} 
	}
	
	obs.data = data.frame(
			pNodeID=rep(0:(nPNodes-1), each=nrow(tre.map)),
			tID=rep(tre.map$tID, times=nPNodes),
			nPositives=as.vector(data.nPositives),
			nTrials=as.vector(data.nTrials)
	);
	
	out = list(
			pop.dag=data.frame(nodeID=as.integer(pop.dag$nodeID), parentID=as.integer(pop.dag$parentID)),
			tre.dag=data.frame(nodeID=as.integer(tre.dag$nodeID), parentID=as.integer(tre.dag$parentID)),
			tre.map=tre.map,
			obs.data=obs.data,
			param=param
	);
	return(out);
}
