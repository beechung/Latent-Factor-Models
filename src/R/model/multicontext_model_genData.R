### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
###
### Author: Bee-Chung Chen

# Generate some normal data for testing
#
#   All feature values (x_obs, x_src, x_dst, x_ctx) follow N(0,1)
#
#	Visually checked 12/20, 2010
#

# Generate a data with just intercept + alpha_i + beta_j without factors or cold-start models, for debugging purpose
genMainEffectData <- function(
  nSrcNodes, nDstNodes, nObs, intercept=1,
  var_y, var_alpha, var_beta, binary.response=FALSE
) {
  alpha = rnorm(nSrcNodes, mean=0, sd=sqrt(var_alpha));
  beta = rnorm(nDstNodes, mean=0, sd=sqrt(var_beta));
  src.id = c(1:nSrcNodes, sample.int(nSrcNodes, nObs-nSrcNodes, replace=TRUE));
  dst.id = c(sample.int(nDstNodes, nObs-nDstNodes, replace=TRUE), 1:nDstNodes);
  output=list();
  output$obs = data.frame(src.id=as.integer(src.id), dst.id=as.integer(dst.id));
  pred.y = intercept + alpha[src.id] + beta[dst.id];
  if(binary.response){
    output$y.prob = 1/(1+exp(-pred.y));
    output$obs$y  = rbinom(n=nObs,size=1,prob=output$y.prob);
  }else{
    output$obs$y =  pred.y + rnorm(nObs, mean=0, sd=sqrt(var_y));
  }

  # Placeholder for features, intercept only
  x_obs = matrix(1, nrow=nObs, ncol=1);
  x_src = matrix(1, nrow=nSrcNodes, ncol=2);
  x_dst = matrix(1, nrow=nDstNodes, ncol=2);

  output$alpha = alpha
  output$beta = beta
  output$intercept = intercept
  output$feature = list(x_src=x_src, x_dst=x_dst, x_obs=x_obs);
  output
}

genNormalData <- function(
    nSrcNodes, nDstNodes, nObs,
	nSrcContexts, nDstContexts, nEdgeContexts, nFactors, has.u, has.gamma, nLocalFactors,
    b, g0, d0, h0=NULL, G=NULL, D=NULL, H=NULL, q=NULL, r=NULL,
	var_y=NULL, var_alpha, var_beta, var_gamma=NULL, var_u=NULL, var_v=NULL, var_w=NULL,
	var_alpha_global=NULL, var_beta_global=NULL,
	has.intercept,
	binary.response=FALSE,
	binary.features=FALSE,
	sparse.matrices=FALSE, index.value.format=FALSE,
	frac.zeroFeatures=0, y.bias=0
){
    if(!is.vector(b))  stop("b should be a vector");

	nObsFeatures = length(b);
	nSrcFeatures = nrow(g0);
	nDstFeatures = nrow(d0);
	nCtxFeatures = nrow(H);

	if(nFactors > 0){
		if(has.u && (is.null(G) || nrow(G) != nSrcFeatures || ncol(G) != nFactors || is.null(var_u))) stop("some problem with G");
		if(is.null(D) || nrow(D) != nDstFeatures || ncol(D) != nFactors || is.null(var_v)) stop("some problem with D");
	}
	if(nEdgeContexts > 0){
		if(is.null(H) || nrow(H) != nCtxFeatures || ncol(H) != nFactors || is.null(var_w)) stop("some problem with H");
	}
	if(nSrcContexts > 1){
		if(is.null(q) || length(q) != nSrcContexts || is.null(var_alpha_global)) stop("some problem with q");
		if(length(var_alpha) == 1) var_alpha = rep(var_alpha, nSrcContexts);
	}
	if(nDstContexts > 1){
		if(is.null(r) || length(r) != nDstContexts || is.null(var_beta_global)) stop("some problem with r");
		if(length(var_beta) == 1) var_beta = rep(var_beta, nDstContexts);
	}
	if(has.gamma){
		if(is.null(h0) || length(h0) != nCtxFeatures || nEdgeContexts <= 0) stop("some problem with h0");
		if(is.null(var_gamma)) stop("some problem with var_gamma");
	}
	if(nLocalFactors > 0){
		if(nLocalFactors*nEdgeContexts != nFactors) stop("nLocalFactors*nEdgeContexts != nFactors");
	}

    if(nObs < nSrcNodes || nObs < nDstNodes) stop("nObs < nSrcNodes || nObs < nDstNodes");

    x_obs = matrix(rnorm(nObs*nObsFeatures), nrow=nObs, ncol=nObsFeatures, dimnames=list(NULL, sprintf("x_obs_%03d", 1:nObsFeatures)));
    x_src = matrix(rnorm(nSrcNodes*nSrcFeatures), nrow=nSrcNodes, ncol=nSrcFeatures,  dimnames=list(NULL, sprintf("x_src_%03d", 1:nSrcFeatures)));
    x_dst = matrix(rnorm(nDstNodes*nDstFeatures), nrow=nDstNodes, ncol=nDstFeatures,  dimnames=list(NULL, sprintf("x_dst_%03d", 1:nDstFeatures)));
	x_ctx = NULL;
	if(nEdgeContexts > 0) x_ctx = matrix(rnorm(nEdgeContexts*nCtxFeatures), nrow=nEdgeContexts, ncol=nCtxFeatures,  dimnames=list(NULL, sprintf("x_ctx_%03d", 1:nCtxFeatures)));

	if(frac.zeroFeatures > 0){
		x_obs[sample.int(length(x_obs), frac.zeroFeatures*length(x_obs))] = 0;
		x_src[sample.int(length(x_src), frac.zeroFeatures*length(x_src))] = 0;
		x_dst[sample.int(length(x_dst), frac.zeroFeatures*length(x_dst))] = 0;
		if(!is.null(x_ctx)) x_ctx[sample.int(length(x_ctx), frac.zeroFeatures*length(x_ctx))] = 0;
	}
	if(binary.features){
		x_obs[,] = x_obs > 0;
		x_src[,] = x_src > 0;
		x_dst[,] = x_dst > 0;
		if(!is.null(x_ctx)) x_ctx[,] = x_ctx > 0;
	}
	if(has.intercept){
		x_obs[,1] = 1;
		x_src[,1] = 1;
		x_dst[,1] = 1;
		if(!is.null(x_ctx)) x_ctx[,1] = 1;
	}

	alpha = as.matrix(x_src %*% g0);
	alpha_global = NULL;
	if(nSrcContexts > 1){
		alpha_global = rnorm(nSrcNodes, mean=0, sd=sqrt(var_alpha_global));
		for(k in 1:nSrcContexts){
			alpha[,k] = alpha[,k] + q[k]*alpha_global + rnorm(nSrcNodes,mean=0,sd=sqrt(var_alpha[k]));
		}
	}else{
		alpha = alpha + rnorm(nSrcNodes, mean=0, sd=sqrt(var_alpha));
	}

	beta = as.matrix(x_dst %*% d0);
	beta_global = NULL;
	if(nDstContexts > 1){
		beta_global = rnorm(nDstNodes, mean=0, sd=sqrt(var_beta_global));
		for(k in 1:nDstContexts){
			beta[,k] = beta[,k] + r[k]*beta_global + rnorm(nDstNodes,mean=0,sd=sqrt(var_beta[k]));
		}
	}else{
		beta = beta + rnorm(nDstNodes, mean=0, sd=sqrt(var_beta));
	}

	if(has.gamma){
		gamma = drop(as.matrix(x_ctx %*% h0) + rnorm(nEdgeContexts, mean=0, sd=sqrt(var_gamma)));
	}

	if(has.u)             u = as.matrix(x_src %*% G) + rnorm(nSrcNodes*nFactors, mean=0, sd=sqrt(var_u));
	if(nFactors > 0)      v = as.matrix(x_dst %*% D) + rnorm(nDstNodes*nFactors, mean=0, sd=sqrt(var_v));
	if(nLocalFactors > 0){
		w = array(0.0, dim=c(nEdgeContexts,nFactors));
		for(k in 1:nEdgeContexts) w[k, select.factor.indices(k,nEdgeContexts,nFactors) ] = 1;
	}else{
		if(nEdgeContexts > 0) w = as.matrix(x_ctx %*% H) + rnorm(nEdgeContexts*nFactors, mean=0, sd=sqrt(var_w));
	}

    src.id = c(1:nSrcNodes, sample.int(nSrcNodes, nObs-nSrcNodes, replace=TRUE));
    dst.id = c(sample.int(nDstNodes, nObs-nDstNodes, replace=TRUE), 1:nDstNodes);
	src.context = NULL;
	if(nSrcContexts > 1) src.context = sample.int(nSrcContexts, nObs, replace=TRUE);
	dst.context = NULL;
	if(nDstContexts > 1) dst.context = sample.int(nDstContexts, nObs, replace=TRUE);
	edge.context = NULL;
	if(nEdgeContexts > 0) edge.context = sample.int(nEdgeContexts, nObs, replace=TRUE);

	output=list();
	output$obs = data.frame(src.id=as.integer(src.id), dst.id=as.integer(dst.id));
	output$obs$src.context  = src.context;
	output$obs$dst.context  = dst.context;
	output$obs$edge.context = edge.context;

	output$feature = list(x_src=x_src, x_dst=x_dst, x_obs=x_obs, x_ctx=x_ctx);

	output$factor  = list(alpha=alpha, beta=beta);
	if(nSrcContexts > 1)  output$factor$alpha_global = alpha_global;
	if(nDstContexts > 1)  output$factor$beta_global  = beta_global;
	if(has.gamma)         output$factor$gamma = gamma;
	if(has.u)             output$factor$u = u;
	if(nFactors > 0)      output$factor$v = v;
	if(nEdgeContexts > 0) output$factor$w = w;

	param = list(b=b, g0=g0, d0=d0);
	if(has.gamma) param$h0 = h0;
	if(nFactors > 0){
		if(has.u) param$G = G;
		param$D = D;
		if(nEdgeContexts > 0) param$H = H;
	}
	if(nSrcContexts > 1) param[["q"]] = q;
	if(nDstContexts > 1) param[["r"]] = r;

	param$var_y = var_y;  param$var_alpha = var_alpha;  param$var_beta = var_beta;
	if(nSrcContexts > 1) param$var_alpha_global = var_alpha_global;
	if(nDstContexts > 1) param$var_beta_global = var_beta_global;
	if(has.gamma) param$var_gamma = var_gamma;
	if(nFactors > 0){
		if(has.u) param$var_u = var_u;
		param$var_v = var_v;
		if(nEdgeContexts > 0) param$var_w = var_w;
	}

	output$param = param;

	pred.y = y.bias + predict.y.from.factors(output$obs, output$factor, output$feature, output$param);
	if(binary.response){
		output$y.prob = 1/(1+exp(-pred.y));
		output$obs$y  = rbinom(n=nObs,size=1,prob=output$y.prob);
	}else{
		output$obs$y =  pred.y + rnorm(nObs, mean=0, sd=sqrt(var_y));
	}

	if(sparse.matrices){
		if(index.value.format) func = matrix.to.index.value
		else                   func = function(x){ return(Matrix(x, sparse=TRUE)); };
		x_obs = func(x_obs);
		x_src = func(x_src);
		x_dst = func(x_dst);
		if(!is.null(x_ctx)) x_ctx = func(x_ctx);
		output$feature = list(x_src=x_src, x_dst=x_dst, x_obs=x_obs, x_ctx=x_ctx);
	}

    return(output);
}

# Visually checked 12/20, 2010
generate.GaussianData <- function(
    nSrcNodes, nDstNodes, nObs,
	nSrcContexts, nDstContexts, nEdgeContexts, nFactors, has.gamma, has.u,
	nObsFeatures, nSrcFeatures, nDstFeatures, nCtxFeatures=0, nLocalFactors=0,
	b.sd=1, g0.sd=1, d0.sd=1, h0.sd=1, G.sd=1, D.sd=1, H.sd=1, q.sd=1, r.sd=1,
	q.mean=5, r.mean=5,
	var_y, var_alpha, var_beta, var_gamma=NULL, var_u=NULL, var_v=NULL, var_w=NULL,
	var_alpha_global=NULL, var_beta_global=NULL,
	has.intercept=TRUE, y.bias=0,
	binary.features=FALSE, binary.response=FALSE,
	sparse.matrices=FALSE, index.value.format=FALSE, frac.zeroFeatures=0
){
	b  = rnorm(nObsFeatures, mean=0, sd=b.sd);
	g0 = matrix(rnorm(nSrcFeatures*nSrcContexts, mean=0, sd=g0.sd),nrow=nSrcFeatures);
	d0 = matrix(rnorm(nDstFeatures*nDstContexts, mean=0, sd=d0.sd),nrow=nDstFeatures);
	G=NULL; D=NULL; H=NULL; h0=NULL;
	if(has.gamma) h0 = rnorm(nCtxFeatures, mean=0, sd=h0.sd);
	if(nFactors > 0){
		if(has.u){ G = matrix(rnorm(nSrcFeatures*nFactors, mean=0, sd=G.sd), nrow=nSrcFeatures); }
		else{ if(nSrcNodes > nDstNodes) stop("has.u && nSrcNodes > nDstNodes"); }
		D = matrix(rnorm(nDstFeatures*nFactors, mean=0, sd=D.sd), nrow=nDstFeatures);
		if(nEdgeContexts > 0) H = matrix(rnorm(nCtxFeatures*nFactors, mean=0, sd=H.sd), nrow=nCtxFeatures);
	}
	q = NULL; r = NULL;
	if(nSrcContexts > 1) q = rnorm(nSrcContexts, mean=q.mean, sd=q.sd);
	if(nDstContexts > 1) r = rnorm(nDstContexts, mean=r.mean, sd=r.sd);

	ans = genNormalData(
	    nSrcNodes=nSrcNodes, nDstNodes=nDstNodes, nObs=nObs,
		nSrcContexts=nSrcContexts, nDstContexts=nDstContexts, nEdgeContexts=nEdgeContexts,
		nFactors=nFactors, nLocalFactors=nLocalFactors, has.u=has.u, has.gamma=has.gamma,
	    b=b, g0=g0, d0=d0, h0=h0, G=G, D=D, H=H, q=q, r=r,
		var_y=var_y, var_alpha=var_alpha, var_beta, var_gamma=var_gamma, var_u=var_u, var_v=var_v, var_w=var_w,
		var_alpha_global=var_alpha_global, var_beta_global=var_beta_global,
		has.intercept=has.intercept, binary.response=binary.response, binary.features=binary.features,
		sparse.matrices=sparse.matrices, frac.zeroFeatures=frac.zeroFeatures, y.bias=y.bias,
		index.value.format=index.value.format
	);
	return(ans);
}

