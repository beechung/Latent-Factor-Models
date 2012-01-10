### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

#####################################################################
### Index The Input Data
#####################################################################
###
### Index src nodes, dst nodes, contexts ...
### INPUT: obs   = data.frame(src_id, dst_id, src_context, dst_context, ctx_id, y);
###                where ctx_id is the ID of the edge context
###                DENSE FORMAT                                      SPARSE FORMAT
###        x_obs = data.frame(feature1, feature2, ...)           or  data.frame(obs_id, index, value)
###        x_src = data.frame(src_id, feature1, feature2, ...)   or  data.frame(src_id, index, value)
###        x_dst = data.frame(dst_id, feature1, feature2, ...)   or  data.frame(dst_id, index, value)
###        x_ctx = data.frame(ctx_id, feature1, feature2, ...)   or  data.frame(ctx_id, index, value)
### 
### Note: obs_id=n means the nth row of the obs table
###
indexData <- function(
	obs,          # observation
	src.dst.same, # whether src nodes and dst nodes are from the same set of nodes
	rm.self.link=FALSE, # whether to remove self link
	x_obs=NULL, x_src=NULL, x_dst=NULL, x_ctx=NULL,
	add.intercept=TRUE,   # whether to add an intercept
	SrcIDs=NULL, DstIDs=NULL, SrcContexts=NULL, DstContexts=NULL, CtxIDs=NULL,
	# E.g., out$x_src[i,] should correspond to SrcIDs[i]
	rm.obs.mismatch=FALSE,# whether to remove observations that cannot find the specified IDs
	other.columns=NULL
){
	if(!is.data.frame(obs)) stop("Please check input parameter 'obs' when calling function indexData or indexTestData: obs should be a data frame.");
	if(nrow(obs) == 0) stop("Please check input parameter 'obs' when calling function indexData or indexTestData: obs contains no data (i.e., nrow(obs) == 0).");
	
	check_names(obs,"obs",required = c("src_id", "dst_id", "y"),prefix="Please check input parameter 'obs' when calling function indexData or indexTestData. ");
	if(!is.null(x_src)){
		if(is.null(x_src$src_id)) stop("Please check input parameter 'x_src' when calling function indexData or indexTestData: src_id must be a column in x_src");
		if(!is.sparse.feature(x_src) && length(unique(x_src$src_id)) != nrow(x_src)) stop("Please check input parameter 'x_src' when calling function indexData or indexTestData: x_src$src_id contains duplicate IDs.");
		if(is.factor(x_src$src_id)) x_src$src_id = as.character(x_src$src_id);
	} 
	if(!is.null(x_dst)){
		if(is.null(x_dst$dst_id)) stop("Please check input parameter 'x_dst' when calling function indexData or indexTestData: dst_id must be a column in x_dst");
		if(!is.sparse.feature(x_dst) && length(unique(x_dst$dst_id)) != nrow(x_dst)) stop("Please check input parameter 'x_dst' when calling function indexData or indexTestData: x_dst$dst_id contains duplicate IDs.");
		if(is.factor(x_dst$dst_id)) x_dst$dst_id = as.character(x_dst$dst_id);
	}
	if(!is.null(x_ctx)){
		if(is.null(x_ctx$ctx_id)) stop("Please check input parameter 'x_ctx' when calling function indexData or indexTestData: ctx_id must be a column in x_ctx");
		if(!is.sparse.feature(x_ctx) && length(unique(x_ctx$ctx_id)) != nrow(x_ctx)) stop("Please check input parameter 'x_ctx' when calling function indexData or indexTestData: x_ctx$ctx_id contains duplicate IDs.");
		if(is.factor(x_ctx$ctx_id)) x_ctx$ctx_id = as.character(x_ctx$ctx_id);
	}
	
	if(!is.null(x_obs) && ncol(x_obs) == 3 && all(c("obs_id", "index", "value") %in% names(x_obs))){
		# x_obs is in the sparse format
		if(!is.numeric(x_obs$obs_id)) stop("Please check input parameter 'x_obs' when calling function indexData or indexTestData: x_obs$obs_id must be numeric, specifying row numbers of table obs.");
		if(any(x_obs$obs_id <= 0)) stop("Please check input parameter 'x_obs' when calling function indexData or indexTestData: x_obs$obs_id must be > 0, specifying row numbers of table obs (starting from 1, not 0).");
		x_obs = get.feature.matrix(x=x_obs, id.colname="obs_id", selected.id=1:nrow(obs), add.intercept=add.intercept, err.prefix="Please check input parameter 'x_obs' when calling function indexData or indexTestData: ", err.x.name="x_obs", err.select.name="1:nrow(obs)");
	}else{
		# x_obs is in the dense format
		if(add.intercept){ regFormula = formula(~.);   default = 1.0; }
		else{              regFormula = formula(~.-1); default = 0.0; }
		if(is.null(x_obs)) x_obs = matrix(default, nrow=nrow(obs), ncol=1)
		else               x_obs = model.matrix(regFormula, x_obs);
	}
	
	if(nrow(obs) != nrow(x_obs)) stop("Please check input parameters 'obs' and 'x_obs' when calling function indexData or indexTestData: nrow(obs) != the number of observations in x_obs. When x_obs is in the dense format, please make sure x_obs[i,] represents the feature vector for obs[i,]. When x_obs is in the sparse format, please make sure x_obs[i,] represents a non-zero feature for obs[x_obs[i,'obs_id'],]");

	for(name in c("src_id", "dst_id", "src_context", "dst_context", "ctx_id")){
		if(!is.null(obs[[name]]) && is.factor(obs[[name]])) obs[[name]] = as.character(obs[[name]]);
	}
	if(rm.self.link){
		if(!src.dst.same) stop("Please check input parameters 'rm.self.link' and 'src.dst.same' when calling function indexData or indexTestData: When rm.self.link = TRUE, src.dst.same must also = TRUE");
		select = obs$src_id != obs$dst_id;
		if(sum(select) != nrow(obs)){
			obs   = obs[select,];
			x_obs = x_obs[select,,drop=FALSE];
		}
	}
	
	if(src.dst.same){
		if(!is.null(SrcIDs)){
			if(is.null(DstIDs)) stop("!is.null(SrcIDs) && is.null(DstIDs)");
			# both non-null
			if(length(SrcIDs) != length(DstIDs)) stop("length(SrcIDs) != length(DstIDs)");
			if(any(SrcIDs != DstIDs)) stop("any(SrcIDs != DstIDs)");
		}else{
			if(!is.null(DstIDs)) stop("is.null(SrcIDs) && !is.null(DstIDs)");
			# both null
			SrcIDs = sort(unique(c(obs$src_id, obs$dst_id)));
			DstIDs = SrcIDs;
		}
	}else{
		if(is.null(SrcIDs)) SrcIDs = sort(unique(obs$src_id));
		if(is.null(DstIDs)) DstIDs = sort(unique(obs$dst_id));
	}
	if(!is.null(obs$src_context) && is.null(SrcContexts)) SrcContexts = sort(unique(obs$src_context));
	if(!is.null(obs$dst_context) && is.null(DstContexts)) DstContexts = sort(unique(obs$dst_context));
	if(!is.null(obs$ctx_id) && is.null(CtxIDs)) CtxIDs = sort(unique(obs$ctx_id));

	out = list();
	
	# Observations
	out$obs = data.frame(
			src.id=match(obs$src_id, SrcIDs),
			dst.id=match(obs$dst_id, DstIDs)
	);
	if(!is.null(obs$src_context)) out$obs$src.context = match(obs$src_context, SrcContexts);
	if(!is.null(obs$dst_context)) out$obs$dst.context = match(obs$dst_context, DstContexts);
	if(!is.null(obs$ctx_id)) out$obs$edge.context = match(obs$ctx_id, CtxIDs);
	out$obs$y = as.double(obs$y);
	
	if(!is.null(other.columns)){
		for(k in 1:length(other.columns)) out$obs[,other.columns[k]] = obs[,other.columns[k]];
	}
	
	# check NA
	select = rep(TRUE, nrow(out$obs));
	in.name  = c("src_id", "dst_id", "src_context", "dst_context", "ctx_id");
	out.name = c("src.id", "dst.id", "src.context", "dst.context", "edge.context");
	id.name  = c("SrcIDs", "DstIDs", "SrcContexts", "DstContexts", "CtxIDs");
	for(i in 1:length(in.name)){
		name = out.name[i];
		if(!is.null(out$obs[[name]])){
			if(rm.obs.mismatch){
				select = select & (!is.na(out$obs[[name]]));
			}else{
				if(any(is.na(out$obs[[name]]))) stop("Please check input parameters 'obs' and '",id.name[i],"' when calling function indexData: Some IDs in obs$",in.name[i]," cannot be found in ",id.name[i]);
			}
		}
	}
	if(rm.obs.mismatch && !all(select)){
		out$obs = out$obs[select,];
		x_obs   = x_obs[select,,drop=FALSE];
	} 
	
	out$IDs = list(SrcIDs=SrcIDs, DstIDs=DstIDs, SrcContexts=SrcContexts, DstContexts=DstContexts, CtxIDs=CtxIDs);
	
	# Features
	out$feature = list(x_obs=x_obs);
	out$feature$x_src = get.feature.matrix(
		x=x_src, id.colname="src_id", selected.id=SrcIDs, add.intercept=add.intercept, 
		err.prefix="Please check input parameters 'x_src' and 'obs' when calling function indexData or indexTestData: ", err.x.name="x_src", err.select.name="obs$src_id"
	);
	out$feature$x_dst = get.feature.matrix(
		x=x_dst, id.colname="dst_id", selected.id=DstIDs, add.intercept=add.intercept, 
		err.prefix="Please check input parameters 'x_dst' and 'obs' when calling function indexData or indexTestData: ", err.x.name="x_dst", err.select.name="obs$dst_id"
	);
	
	if(!is.null(out$obs$edge.context)){
		out$feature$x_ctx = get.feature.matrix(
			x=x_ctx, id.colname="ctx_id", selected.id=CtxIDs, add.intercept=add.intercept,
			err.prefix="Please check input parameters 'x_ctx' and 'obs' when calling function indexData or indexTestData: ", err.x.name="x_ctx", err.select.name="obs$ctx_id"
		);
	}else{
		if(!is.null(x_ctx)) stop("Please check input parameters 'x_ctx' and 'obs' when calling function indexData or indexTestData: When obs$ctx_id = NULL, x_ctx must also = NULL");
	}
	out$add.intercept   = add.intercept;
	out$src.dst.same    = src.dst.same;
	out$rm.self.link    = rm.self.link;
	out$rm.obs.mismatch = rm.obs.mismatch;
	return(out);
}

indexTestData <- function(
	data.train,
	obs, x_obs=NULL, x_src=NULL, x_dst=NULL, x_ctx=NULL,
	other.columns=NULL
){
	if(data.train$src.dst.same){
		SrcIDs = c(data.train$IDs$SrcIDs, setdiff(unique(c(obs$src_id, obs$dst_id)), data.train$IDs$SrcIDs));
		DstIDs = SrcIDs;
	}else{
		SrcIDs = c(data.train$IDs$SrcIDs, setdiff(unique(obs$src_id), data.train$IDs$SrcIDs));
		DstIDs = c(data.train$IDs$DstIDs, setdiff(unique(obs$dst_id), data.train$IDs$DstIDs));
	}
	if(!is.null(obs$src_context)) SrcContexts = c(data.train$IDs$SrcContexts, setdiff(unique(obs$src_context), data.train$IDs$SrcContexts))
	else                          SrcContexts = NULL;
	if(!is.null(obs$dst_context)) DstContexts = c(data.train$IDs$DstContexts, setdiff(unique(obs$dst_context), data.train$IDs$DstContexts))
	else                          DstContexts = NULL;
	if(!is.null(obs$ctx_id)) CtxIDs = c(data.train$IDs$CtxIDs, setdiff(unique(obs$ctx_id), data.train$IDs$CtxIDs))
	else                     CtxIDs = NULL;
	
	data.test = indexData(
		obs=obs, src.dst.same=data.train$src.dst.same, rm.self.link=data.train$rm.self.link,
		x_obs=x_obs, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx,
		add.intercept=data.train$add.intercept,
		SrcIDs=SrcIDs, DstIDs=DstIDs, SrcContexts=SrcContexts, DstContexts=DstContexts, CtxIDs=CtxIDs,
		rm.obs.mismatch=data.train$rm.obs.mismatch, other.columns=data.train$other.columns
	);
	return(data.test);
}

###
### Get subset information
###	
###		subset.info$src.context[[k]]:  a vector of src node IDs that have data in src context k
###		subset.info$dst.context[[k]]:  a vector of dst node IDs that have data in dst context k
###		subset.info$edge.context[[k]]: a vector of node IDs that have data in edge context k
###		subset.info$src.id: a vector of src node IDs that have data in obs
###		subset.info$dst.id: a vector of dst node IDs that have data in obs
###		subset.info$any.id: a vector of node IDs that have data in obs
###
get.subset.info <- function(obs, output.any=FALSE){
	subset.info = list(
			src.id  = unique(obs$src.id),
			dst.id  = unique(obs$dst.id)
	);
	if(output.any) subset.info$any.id = unique(subset.info$src.id, subset.info$dst.id);
	if(!is.null(obs$edge.context)) subset.info$edge.id = unique(obs$edge.context);
	
	if(!is.null(obs$src.context)){
		src = tapply(obs$src.id, list(obs$src.context), FUN=unique, simplify=FALSE);
		if(any(names(src) != 1:length(src))) stop("Some src.context has no data");
		subset.info$src.context = src;
	}else{
		subset.info$src.context = list(subset.info$src.id);
	}
	if(!is.null(obs$dst.context)){
		dst = tapply(obs$dst.id, list(obs$dst.context), FUN=unique, simplify=FALSE);
		if(any(names(dst) != 1:length(dst))) stop("Some dst.context has no data");
		subset.info$dst.context = dst;
	}else{
		subset.info$dst.context = list(subset.info$dst.id);
	}
	if(output.any){
		if(!is.null(obs$edge.context)){
			src = tapply(obs$src.id, list(obs$edge.context), FUN=unique, simplify=FALSE);
			dst = tapply(obs$dst.id, list(obs$edge.context), FUN=unique, simplify=FALSE);
			if(any(names(src) != 1:length(src))) stop("Some edge.context has src nodes");
			if(any(names(dst) != 1:length(dst))) stop("Some edge.context has dst nodes");
			if(length(src) != length(dst)) stop("length(src) != length(dst)");
			subset.info$edge.context = list();
			for(k in 1:length(src)){
				subset.info$edge.context[[k]] = unique(src[[k]], dst[[k]]);
			}
		}else{
			subset.info$edge.context = list(subset.info$any.id);
		}
	}
	return(subset.info);
}

###
### Functions for the Local Interaction Factors
###
select.factor.indices <- function(k, nContexts, nFactors){
	nLocalFactors = as.integer(nFactors / nContexts);
	if(nLocalFactors * nContexts != nFactors) stop("nLocalFactors * nContexts != nFactors");
	if(k < 1) stop("k < 1");
	if(k > nContexts) stop("k > nContexts");
	return((k-1)*nLocalFactors + (1:nLocalFactors));
}
get.context.index <- function(f, nContexts, nFactors){
	nLocalFactors = as.integer(nFactors / nContexts);
	if(nLocalFactors * nContexts != nFactors) stop("nLocalFactors * nContexts != nFactors");
	return(as.integer((f-1)/nLocalFactors)+1);
}

#####################################################################
### Functions to do Sanity Check
#####################################################################

###
### Syntactic check of model specification
###
###   factor  = list(alpha, beta, gamma, u, v, w); # Initial factor values
###   obs     = data.frame(src.id, dst.id, src.context, dst.context, edge.context, y);
###   feature = list(x_obs, x_src, x_dst, x_ctx);
###   param   = list(b, g0, d0, h0, G, D, H, q, r,
###                  var_y, var_alpha, var_alpha_global, var_beta, var_beta_global, var_gamma, var_u, var_v, var_w);
###
syncheck.multicontext.spec <- function(
		factor, obs, feature, param, warning=0, print=FALSE,
		factor.name.required = c("alpha", "beta"),
		factor.name.optional = c("gamma", "u", "v", "w"),
		factor.name.allowed  = c("fScore", "alpha_global", "beta_global"),
		obs.name.required = c("src.id", "dst.id", "y"),
		obs.name.optional = c("src.context", "dst.context", "edge.context"),
		obs.name.allowed  = c("response"),
		feature.name.required = c("x_obs", "x_src", "x_dst"),
		feature.name.optional = c("x_ctx"),
		feature.name.allowed  = c(),
		param.name.required = c("b", "g0", "d0", "var_y", "var_alpha", "var_beta"),
	    param.name.optional = c("h0", "G", "D", "H", "q", "r", "var_alpha_global","var_beta_global", "var_gamma", "var_u", "var_v", "var_w"),
		param.name.allowed  = c("nLocalFactors", "approx.interaction", "xi", "reg.algo", "reg.control", "is.logistic")
){
	factor.name.all = c(factor.name.required, factor.name.optional, factor.name.allowed);
	obs.name.all = c(obs.name.required, obs.name.optional, obs.name.allowed);
	feature.name.all = c(feature.name.required, feature.name.optional, feature.name.allowed);
	param.name.all = c(param.name.required, param.name.optional, param.name.allowed);
	
	warning.any.not.in(names(factor),  factor.name.all, "You specified the following unnecessary components in factor: ",  stop=TRUE);
	warning.any.not.in(names(obs),     obs.name.all,    "You specified the following unnecessary components in obs: ",     stop=TRUE);
	warning.any.not.in(names(feature), feature.name.all,"You specified the following unnecessary components in feature: ", stop=TRUE);
	warning.any.not.in(names(param),   param.name.all,  "You specified the following unnecessary components in param: ",   stop=TRUE);

	warning.any.not.in(factor.name.required,  names(factor),  "You did not specify the following components in factor: ",  stop=TRUE);
	warning.any.not.in(obs.name.required,     names(obs),     "You did not specify the following components in obs: ",     stop=TRUE);
	warning.any.not.in(feature.name.required, names(feature), "You did not specify the following components in feature: ", stop=TRUE);
	warning.any.not.in(param.name.required,   names(param),   "You did not specify the following components in param: ",   stop=TRUE);
	
	if(warning > 2){
		warning.any.not.in(c(factor.name.required,  factor.name.optional),  names(factor),  "You did not specify the following components in factor: ", print=print);
		warning.any.not.in(c(obs.name.required,     obs.name.optional),     names(obs),     "You did not specify the following components in obs: ", print=print);
		warning.any.not.in(c(feature.name.required, feature.name.optional), names(feature), "You did not specify the following components in feature: ", print=print);
		warning.any.not.in(c(param.name.required,   param.name.optional),   names(param),   "You did not specify the following components in param: ", print=print);
	}
	out = list(
			nObs=get.size(nrow(obs), nrow(feature[["x_obs"]])),
			nSrcNodes=get.size(nrow(factor[["alpha"]]), nrow(factor[["u"]]), nrow(feature[["x_src"]])),
			nDstNodes=get.size(nrow(factor[["beta"]]),  nrow(factor[["v"]]), nrow(feature[["x_dst"]])),
			nSrcContexts=get.size(ncol(factor[["alpha"]]), length(param[["q"]])),
			nDstContexts=get.size(ncol(factor[["beta"]]),  length(param[["r"]])),
			nEdgeContexts=get.size(nrow(factor[["w"]]), nrow(feature[["x_ctx"]])),
			nFactors=get.size(ncol(factor[["u"]]), ncol(factor[["v"]]), ncol(factor[["w"]]), ncol(param[["G"]]), ncol(param[["D"]]), ncol(param[["H"]])),
			nObsFeatures=get.size(ncol(feature[["x_obs"]])),
			nSrcFeatures=get.size(ncol(feature[["x_src"]])),
			nDstFeatures=get.size(ncol(feature[["x_dst"]])),
			nCtxFeatures=get.size(ncol(feature[["x_ctx"]])),
			nVar_y=as.integer(length(param[["var_y"]])), 
			nVar_alpha=as.integer(length(param[["var_alpha"]])), 
			nVar_alpha_global=as.integer(length(param[["var_alpha_global"]])), 
			nVar_beta=as.integer(length(param[["var_beta"]])), 
			nVar_beta_global=as.integer(length(param[["var_beta_global"]])), 
			nVar_gamma=as.integer(length(param[["var_gamma"]])),
			nVar_u=as.integer(length(param[["var_u"]])),
			nVar_v=as.integer(length(param[["var_v"]])),
			nVar_w=as.integer(length(param[["var_w"]])),
			nrowU=get.size(nrow(factor[["u"]])),
			nGamma=as.integer(length(factor[["gamma"]]))
	);
	
	if(out$nSrcContexts > 1){
		if(is.null(obs[["src.context"]])) stop("obs$src.context cannot be null");
		if(out$nSrcContexts < max(obs[["src.context"]])) stop("out$nSrcContexts < max(obs$src.context)");
		check.individual("param$q", param[["q"]], "double", out$nSrcContexts, isNullOK=FALSE, stopIfAnyNull=list("factor$alpha"=factor$alpha));
		check.individual("param$var_alpha_global", param$var_alpha_global, "double", list(1,out$nSrcNodes), isNullOK=FALSE);
	}
	if(out$nDstContexts > 1){
		if(is.null(obs[["dst.context"]])) stop("obs$dst.context cannot be null");
		if(out$nDstContexts < max(obs[["dst.context"]])) stop("out$nDstContexts < max(obs$dst.context)");
		check.individual("param$r", param[["r"]], "double", out$nDstContexts, isNullOK=FALSE, stopIfAnyNull=list("factor$beta"=factor$beta));
		check.individual("param$var_beta_global", param$var_beta_global, "double", list(1,out$nDstNodes), isNullOK=FALSE);
	}
	if(out$nEdgeContexts > 1){
		if(is.null(obs[["edge.context"]])) stop("obs$edge.context cannot be null");
		if(out$nEdgeContexts < max(obs[["edge.context"]])) stop("out$nEdgeContexts < max(obs$edge.context)");
		check.individual("param$H", param$H, "model", size=NULL, isNullOK=FALSE, stopIfAnyNull=list("factor$w"=factor$w));
		check.individual("param$var_w", param$var_w, "double", 1, isNullOK=FALSE);
	}
	
	check.individual(
		"factor$alpha", factor$alpha, "double", c(out$nSrcNodes, out$nSrcContexts), isNullOK=FALSE,
		stopIfAnyNull=list("obs$y"=obs$y,"obs$src.id"=obs$src.id,"feature$x_src"=feature$x_src,"param$g0"=param$g0,"param$var_alpha"=param$var_alpha)
	);
	check.individual(
		"factor$beta", factor$beta, "double", c(out$nDstNodes, out$nDstContexts), isNullOK=FALSE,
		stopIfAnyNull=list("obs$y"=obs$y,"obs$dst.id"=obs$dst.id,"feature$x_dst"=feature$x_dst,"param$d0"=param$d0,"param$var_beta"=param$var_beta)
	);
	check.individual(
		"factor$gamma", factor$gamma, "double", out$nEdgeContexts, isNullOK=TRUE,
		stopIfAnyNull=list("obs$y"=obs$y,"obs$edge.context"=obs$edge.context,"feature$x_ctx"=feature$x_ctx,"param$h0"=param$h0,"param$var_gamma"=param$var_gamma)
	);
	check.individual(
		"factor$u", factor$u, "double", c(out$nSrcNodes, out$nFactors), isNullOK=TRUE,
		stopIfAnyNull=list("obs$y"=obs$y,"obs$src.id"=obs$src.id,"feature$x_src"=feature$x_src,"param$G"=param$G,"param$var_u"=param$var_u)
	);
	check.individual(
		"factor$v", factor$v, "double", c(out$nDstNodes, out$nFactors), isNullOK=(out$nFactors==0),
		stopIfAnyNull=list("obs$y"=obs$y,"obs$dst.id"=obs$dst.id,"feature$x_dst"=feature$x_dst,"param$D"=param$D,"param$var_v"=param$var_v)
	);
	check.individual(
		"factor$w", factor$w, "double", c(out$nEdgeContexts, out$nFactors), isNullOK=(out$nEdgeContexts==0),
		stopIfAnyNull=list("obs$y"=obs$y,"obs$edge.context"=obs$edge.context,"feature$x_ctx"=feature$x_ctx,"param$H"=param$H,"param$var_w"=param$var_w)
	);

	check.individual("feature$x_obs", feature$x_obs, c("double", "dgCMatrix"), c(out$nObs, out$nObsFeatures), isNullOK=FALSE, stopIfAnyNull=list("param$b"=param$b));
	check.individual("feature$x_src", feature$x_src, c("double", "dgCMatrix"), c(out$nSrcNodes, out$nSrcFeatures), isNullOK=FALSE, stopIfAnyNull=list("param$g0"=param$g0));
	check.individual("feature$x_dst", feature$x_dst, c("double", "dgCMatrix"), c(out$nDstNodes, out$nDstFeatures), isNullOK=FALSE, stopIfAnyNull=list("param$d0"=param$d0));
	check.individual("feature$x_ctx", feature$x_ctx, c("double", "dgCMatrix"), c(out$nEdgeContexts, out$nCtxFeatures), isNullOK=TRUE, stopIfAnyNull=list("param$H"=param$H));
	
	check.individual("param$b",  param$b,  "model", size=NULL, isNullOK=FALSE, stopIfAnyNull=list("feature$x_obs"=feature$x_obs));
	check.individual("param$g0", param$g0, "model", size=NULL, isNullOK=FALSE, stopIfAnyNull=list("factor$alpha"=factor$alpha));
	check.individual("param$d0", param$d0, "model", size=NULL, isNullOK=FALSE, stopIfAnyNull=list("factor$beta"=factor$beta));
	check.individual("param$h0", param$h0, "model", size=NULL, isNullOK=TRUE,  stopIfAnyNull=list("factor$gamma"=factor$gamma));
	check.individual("param$G",  param$G,  "model", size=NULL, isNullOK=TRUE,  stopIfAnyNull=list("factor$u"=factor$u));
	check.individual("param$D",  param$D,  "model", size=NULL, isNullOK=TRUE,  stopIfAnyNull=list("factor$v"=factor$v));
	
	check.individual("param$var_y", param$var_y, "double", list(1,out$nObs), isNullOK=FALSE, stopIfAnyNull=list("param$b"=param$b));
	check.individual("param$var_alpha", param$var_alpha, "double", out$nSrcContexts, isNullOK=FALSE, stopIfAnyNull=list("factor$alpha"=factor$alpha));
	check.individual("param$var_beta",  param$var_beta,  "double", out$nDstContexts, isNullOK=FALSE, stopIfAnyNull=list("factor$beta"=factor$beta));
	check.individual("param$var_gamma", param$var_gamma, "double", 1, isNullOK=TRUE,  stopIfAnyNull=list("factor$gamma"=factor$gamma));
	check.individual("param$var_u",  param$var_u,  "double", list(1, out$nFactors), isNullOK=TRUE, stopIfAnyNull=list("factor$u"=factor$u));
	check.individual("param$var_v",  param$var_v,  "double", list(1, out$nFactors), isNullOK=TRUE, stopIfAnyNull=list("factor$v"=factor$v));

	if(out$nGamma > 0 && out$nGamma != out$nEdgeContexts)  stop("out$nGamma > 0 && out$nGamma != out$nEdgeContexts");
	if(is.null(factor$u) && out$nFactors > 0 && out$nDstNodes < out$nSrcNodes) stop("factor$u is null && nDstNodes < nSrcNodes");
	
	return(out);
}
warning.any.not.in <- function(a, b, msg, stop=FALSE, print=FALSE){
	if(any(!(a %in% b))){
		temp = a[!(a %in% b)];
		if(stop) stop(msg,paste(temp, collapse=", "))
		else     warning(msg,paste(temp, collapse=", "),call.=FALSE);
		if(print) cat("\nWARNING: ",msg,paste(temp, collapse=", "),"\n\n",sep="");
	}
}
get.size <- function(...){
	temp = c(...);
	if(is.null(temp)) return(as.integer(0));
	s = max(temp);
	if(any(temp[temp != 0] != s)) warning("Not all input elements have the same size");
	return(as.integer(s));
}
check.individual <- function(name, x, type.list, size, isNullOK, stopIfAnyNull=NULL){
	if(is.null(x)){
		if(isNullOK) return(TRUE)
		else stop(name," is null");
	}
	
	type.correct = FALSE;
	for(type in type.list){
		if(type == "double"){
			if(is.double(x)){ type.correct = TRUE; break; }
		}else if(type == "integer" || type == "int"){
			if(is.integer(x)){ type.correct = TRUE; break; }
		}else if(type == "dgCMatrix"){
			if(class(x) == "dgCMatrix"){ type.correct = TRUE; break; }
		}else if(type == "model"){
			type.correct = TRUE; break;
		}else stop("Unknown type: ",type,sep="");
	}
	if(!type.correct) stop("Type of ",name," is not one of ",paste(type.list,collapse=", "));
	
	if(!is.null(stopIfAnyNull)){
		for(i in 1:length(stopIfAnyNull)){
			if(is.null(stopIfAnyNull[[i]])) stop("When ",name," exists, ",names(stopIfAnyNull)[i]," cannot be null");
		}
	}
	
	if(!is.null(size)){
		d = dim(x);
		if(is.null(d)) d = length(x);
		if(is.list(size)){
			dim.correct = FALSE;
			for(i in 1:length(size)){
				if(length(d) == length(size[[i]]) && all(d == size[[i]])){dim.correct = TRUE; break;}
			}
			if(!dim.correct){
				if(length(size) == 1) stop(name," has dimensionality mismatch: (",paste(d,collapse=" x "),") vs (",paste(size[[1]],collapse=" x "),")");
				if(length(size) == 2) stop(name," has dimensionality mismatch: (",paste(d,collapse=" x "),") vs (",paste(size[[1]],collapse=" x "),") or (",paste(size[[2]],collapse=" x "),")");
				if(length(size) == 3) stop(name," has dimensionality mismatch: (",paste(d,collapse=" x "),") vs (",paste(size[[1]],collapse=" x "),") or (",paste(size[[2]],collapse=" x "),") or (",paste(size[[3]],collapse=" x "),")");
				stop(name," has dimensionality mismatch: (",paste(d,collapse=" x "),") vs (",paste(size[[1]],collapse=" x "),") or (",paste(size[[2]],collapse=" x "),") or (",paste(size[[3]],collapse=" x "),") ...");
			}
		}else{
			if(length(d) != length(size) || any(d != size)) stop(name," has dimensionality mismatch: (",paste(d,collapse=" x "),") vs (",paste(size,collapse=" x "),")");
		}
	}
}

check_names <- function(x, display.name, required, should.have=c(), optional=c(), prefix=""){
	names.all = c(required, optional, should.have);
	temp = if(is.data.frame(x)) "column(s)" else "component(s)";
	if(!is.null(optional)) warning.any.not.in(names(x), names.all, paste(prefix,"The following ",temp," should not be in ",display.name,": ", sep=""), stop=TRUE);
	warning.any.not.in(required,  names(x),  paste(prefix,"You must have the following ",temp," in ",display.name,": ",sep=""), stop=TRUE);
	warning.any.not.in(should.have, names(x),paste(prefix,"The following ",temp," are missing in ",display.name,": ",sep=""), stop=FALSE, print=TRUE);
}

check.obs.feature <- function(obs, feature, nSrcContexts=1, nDstContexts=1, nEdgeContexts=0){
	if(nrow(obs) != nrow(feature$x_obs)) stop("nrow(obs) != nrow(feature$x_obs)");
	if(!is.na(nEdgeContexts) && nEdgeContexts > 1){
		if(is.null(obs$edge.context)) stop("is.null(obs$edge.context)");
		if(is.null(feature$x_ctx)) stop("is.null(feature$x_ctx)");
	}
	nSrcNodes = nrow(feature$x_src);
	nDstNodes = nrow(feature$x_dst);
	names   = c( "src.id",  "dst.id", "src.context", "dst.context", "edge.context");
	u.bound = c(nSrcNodes, nDstNodes,  nSrcContexts,  nDstContexts,  nEdgeContexts);
	require = c(     TRUE,      TRUE,nSrcContexts>1,nDstContexts>1,  !is.na(nEdgeContexts) && nEdgeContexts>1);
	for(k in 1:length(names)){
		name = names[k];
		if(is.null(obs[[name]]) && require[k]) stop("obs$",name," is missing");
		if(is.null(obs[[name]])) next;
		if(any(is.na(obs[[name]]))) stop("obs$",name," has NA in it!");
		if(any(obs[[name]] < 1)) stop("obs$",name," has some number < 1");
		if(is.na(u.bound[k])) next;
		if(any(obs[[name]] > u.bound[k])) stop("obs$",name," is out of bound");
	}
}

#####################################################################
### Functions to Make Predictions
#####################################################################

###
### Predict the response using the factors
###
###    src.id > nrow(factor$alpha)   are new src nodes
###    dst.id > nrow(factor$beta)    are new dst nodes
###    edge.context > nrow(factor$w) are new edge contexts
###
predict.y.from.factors <- function(obs, factor, feature, param, ignore.xb=FALSE){
	src.id = obs$src.id;   dst.id = obs$dst.id;
	if(is.null(src.id) || is.null(dst.id)) stop("is.null(src.id) || is.null(dst.id)");
	src.context = obs$src.context; dst.context = obs$dst.context; edge.context = obs$edge.context;
	if(xor(is.null(factor$w),is.null(edge.context))) stop("xor(is.null(factor$w),is.null(edge.context))");
	alpha = get.factors(factor$alpha, param$g0, feature$x_src, reg.algo=param$reg.algo);
	beta  = get.factors(factor$beta,  param$d0, feature$x_dst, reg.algo=param$reg.algo);
	gamma = get.factors(factor$gamma, param$h0, feature$x_ctx, reg.algo=param$reg.algo);
	u     = get.factors(factor$u,     param$G,  feature$x_src, reg.algo=param$reg.algo);
	v     = get.factors(factor$v,     param$D,  feature$x_dst, reg.algo=param$reg.algo);
	w     = get.factors(factor$w,     param$H,  feature$x_ctx, reg.algo=param$reg.algo);
	# xb
	if(ignore.xb) ans = 0
	else          ans = reg.predict(model=param$b, x=feature$x_obs, algo=param$reg.algo);
	# alpha
	if(is.null(src.context)) ans = ans + alpha[src.id,1]
	else                     ans = ans + alpha[cbind(src.id,src.context)];
	# beta
	if(is.null(dst.context)) ans = ans + beta[dst.id,1]
	else                     ans = ans + beta[cbind(dst.id,dst.context)];
	# gamma
	if(!is.null(gamma)){
		if(is.null(edge.context)) stop("!is.null(factor$gamma) && is.null(edge.context)");
		ans = ans + gamma[edge.context];
	}
	# uvw
	temp = c(ncol(u), ncol(v), ncol(w));
	nFactors = max(c(0,temp));
	if(nFactors > 0){
		if(any(temp != nFactors)) stop("factor dimensionality mismatch");
		if(is.null(v))  stop("is.null(factor$v)");
		if(is.null(u))  uvw = v[src.id,,drop=FALSE] * v[dst.id,,drop=FALSE]
		else            uvw = u[src.id,,drop=FALSE] * v[dst.id,,drop=FALSE];
		if(!is.null(w)) uvw = uvw * w[edge.context,,drop=FALSE];
		ans = ans + sum_margin(uvw, 1);
	}
	if(any(is.na(ans))) stop("any(is.na(ans))");
	return(drop(ans));
}

predict.y.from.featuresOnly <- function(obs, factor, feature, param, ignore.xb=FALSE){
	src.id = obs$src.id;   dst.id = obs$dst.id;
	if(is.null(src.id) || is.null(dst.id)) stop("is.null(src.id) || is.null(dst.id)");
	src.context = obs$src.context; dst.context = obs$dst.context; edge.context = obs$edge.context;
	if(xor(is.null(factor$w),is.null(edge.context))) stop("xor(is.null(factor$w),is.null(edge.context))");
	alpha = get.factors.from.features(factor$alpha, param$g0, feature$x_src, reg.algo=param$reg.algo);
	beta  = get.factors.from.features(factor$beta,  param$d0, feature$x_dst, reg.algo=param$reg.algo);
	gamma = get.factors.from.features(factor$gamma, param$h0, feature$x_ctx, reg.algo=param$reg.algo);
	u     = get.factors.from.features(factor$u,     param$G,  feature$x_src, reg.algo=param$reg.algo);
	v     = get.factors.from.features(factor$v,     param$D,  feature$x_dst, reg.algo=param$reg.algo);
	w     = get.factors.from.features(factor$w,     param$H,  feature$x_ctx, reg.algo=param$reg.algo);
	# xb
	if(ignore.xb) ans = 0
	else          ans = reg.predict(model=param$b, x=feature$x_obs, algo=param$reg.algo);
	# alpha
	if(is.null(src.context)) ans = ans + alpha[src.id,1]
	else                     ans = ans + alpha[cbind(src.id,src.context)];
	# beta
	if(is.null(dst.context)) ans = ans + beta[dst.id,1]
	else                     ans = ans + beta[cbind(dst.id,dst.context)];
	# gamma
	if(!is.null(gamma)){
		if(is.null(edge.context)) stop("!is.null(factor$gamma) && is.null(edge.context)");
		ans = ans + gamma[edge.context];
	}
	# uvw
	temp = c(ncol(u), ncol(v), ncol(w));
	nFactors = max(c(0,temp));
	if(nFactors > 0){
		if(any(temp != nFactors)) stop("factor dimensionality mismatch");
		if(is.null(v))  stop("is.null(factor$v)");
		if(is.null(u))  uvw = v[src.id,,drop=FALSE] * v[dst.id,,drop=FALSE]
		else            uvw = u[src.id,,drop=FALSE] * v[dst.id,,drop=FALSE];
		if(!is.null(w)) uvw = uvw * w[edge.context,,drop=FALSE];
		ans = ans + sum_margin(uvw, 1);
	}
	if(any(is.na(ans))) stop("any(is.na(ans))");
	return(drop(ans));
}

###
### Prediction Function
###
###    src.id > nrow(factor$alpha)   are new src nodes
###    dst.id > nrow(factor$beta)    are new dst nodes
###    edge.context > nrow(factor$w) are new edge contexts
###
predict.multicontext <- function(model, obs, feature, is.logistic=FALSE, fScore=NULL){
	nSrcContexts = ncol(model$factor$alpha);
	nDstContexts = ncol(model$factor$beta);
	nEdgeContexts = NA;
	if(!is.null(model$param$nLocalFactors)) nEdgeContexts = ncol(model$factor$v) / model$param$nLocalFactors;
	check.obs.feature(obs, feature, nSrcContexts=nSrcContexts, nDstContexts=nDstContexts, nEdgeContexts=nEdgeContexts);
	if(is.null(fScore)){
		pred.y = predict.y.from.factors(obs=obs, factor=model$factor, feature=feature, param=model$param);
	}else{
		if(length(fScore) != nrow(obs)) stop("length(fScore) != nrow(obs)");
		pred.y = fScore + reg.predict(model=model$param$b, x=feature$x_obs, algo=model$param$reg.algo);
	}
	out = predict.response.from.gaussian(pred.y=pred.y, y=obs$y, param=model$param, is.logistic=is.logistic);
	return(out);
}


# Let F = nrow(factor) and N = nrow(x)
# The factors of cases F+1, ..., N need to be predicted
get.factors <- function(factor, coeff, x, reg.algo=NULL){
	if(is.null(factor)) return(NULL);
	N = nrow(x);
	if(is.matrix(factor)){
		F = nrow(factor); 
	}else if(is.vector(factor)){
		F = length(factor);
	}else stop("unknown type");
	if(N < F) stop("nrow(x) < nrow(factor)");
	if(N == F) return(factor);
	pred = reg.predict(model=coeff, x=x[(F+1):N,,drop=FALSE], algo=reg.algo, ncol=ncol(factor));
	if(is.matrix(factor)) return(rbind(factor, pred));
	if(is.vector(factor)) return(c(factor, pred));
	stop("you should never get here");
}
get.factors.from.features <- function(factor, coeff, x, reg.algo=NULL){
	if(is.null(factor)) return(NULL);
	pred = reg.predict(model=coeff, x=x, algo=reg.algo, ncol=ncol(factor));
	return(pred);
}

#####################################################################
### Likelihood Computation
#####################################################################

get.logLikelihood <- function(
	obs, factor, feature, param, factor.var=NULL, factor.cov=NULL, subset.info=NULL,
	verbose=0, is.logistic=FALSE, prefix="", suffix="", prev.loglik=NULL
){
	b.time = proc.time();
	CD.loglik = complete.logLikelihood(obs, factor, feature, param, is.logistic=is.logistic);
	E.loglik  = NA;
	if(!is.null(factor.var)) E.loglik  = E.logLikelihood(obs, factor, feature, param, factor.var, factor.cov, subset.info);
	time.used = proc.time() - b.time;

	if(verbose >= 1){
		out = sprintf("%s loglik: CD = %13.10g, E = %13.10g (%.2f sec)  %s\n", prefix, CD.loglik, E.loglik, time.used[3], suffix);
		cat(out);
	}
	if(!is.null(prev.loglik)){
		if(E.loglik < prev.loglik$E){
			cat("\nNOTICE: E[loglik] decreases: ",prev.loglik$E," -> ",E.loglik," (may be due to regularized regression)\n\n",sep="");
			warning(prefix," E[loglik] decreases: ",prev.loglik$E," -> ",E.loglik," (may be due to regularized regression)");
		} 
	}
	
	return(list(CD=CD.loglik, E=E.loglik));
}

###
### Compute E[loglikelihood] (only for the Gaussian model)
###
E.logLikelihood <- function(
	obs, factor.mean, feature, param, factor.var, factor.cov, subset.info=NULL
){
	size = syncheck.multicontext.spec(factor=factor.mean, obs=obs, feature=feature, param=param);
	obs.loglik = Eloglik.lm.random.effect(
		coeff=param$b, var=param$var_y,
		response.mean=obs$y-factor.mean$fScore, feature.mean=feature$x_obs,
		response.var=factor.var$fScore, algo=param$reg.algo
	);
	alpha.loglik = 0;
	if(!is.null(param$var_alpha) && any(param$var_alpha > 0)){
		alpha.loglik = E.loglik.mainEffect(
			coeff=param$g0, slope=param[["q"]], var=param$var_alpha, x=feature$x_src, 
			local.mean=factor.mean$alpha, global.mean=factor.mean$alpha_global, 
			local.var =factor.var$alpha,  global.var =factor.var$alpha_global,
			local.cov =factor.cov$alpha,  subset=subset.info$src.context, algo=param$reg.algo
		);
	}
	beta.loglik = 0;
	if(!is.null(param$var_beta) && any(param$var_beta > 0)){
		beta.loglik = E.loglik.mainEffect(
				coeff=param$d0, slope=param[["r"]], var=param$var_beta, x=feature$x_dst, 
				local.mean=factor.mean$beta, global.mean=factor.mean$beta_global, 
				local.var =factor.var$beta,  global.var =factor.var$beta_global,
				local.cov =factor.cov$beta,  subset=subset.info$dst.context, algo=param$reg.algo
		);
	}
	gamma.loglik = 0;
	if(!is.null(param$var_gamma) && any(param$var_gamma > 0)){
		gamma.loglik = E.loglik.mainEffect(
				coeff=param$h0, slope=NULL, var=param$var_gamma, x=feature$x_ctx, 
				local.mean=factor.mean$gamma, global.mean=NULL, 
				local.var =factor.var$gamma,  global.var =NULL,
				local.cov =NULL,              subset=NULL, algo=param$reg.algo
		);
	}
	u.loglik = 0;
	if(!is.null(param$var_u) && any(param$var_u > 0)){
		if(is.null(param$nLocalFactors)) subset = subset.info$src.id
		else                             subset = subset.info$edge.context;
		u.loglik = E.loglik.interaction(
				coeff=param$G, var=param$var_u, x=feature$x_src, 
				factor.mean=factor.mean$u, factor.var=factor.var$u,
				subset=subset, algo=param$reg.algo
		);
	}
	v.loglik = 0;
	if(!is.null(param$var_v) && any(param$var_v > 0)){
		if(!is.null(factor.mean$u)){
			if(is.null(param$nLocalFactors)) subset = subset.info$dst.id
			else                             subset = subset.info$edge.context;
		}else{
			if(is.null(param$nLocalFactors)) subset = subset.info$any.id
			else                             subset = subset.info$edge.context;
		}
		v.loglik = E.loglik.interaction(
				coeff=param$D, var=param$var_v, x=feature$x_dst, 
				factor.mean=factor.mean$v, factor.var=factor.var$v,
				subset=subset, algo=param$reg.algo
		);
	}
	w.loglik = 0;
	if(!is.null(param$var_w) && any(param$var_w > 0)){
		if(!is.null(param$nLocalFactors)) stop("!is.null(param$nLocalFactors)");
		w.loglik = E.loglik.interaction(
				coeff=param$H, var=param$var_w, x=feature$x_ctx, 
				factor.mean=factor.mean$w, factor.var=factor.var$w, 
				subset=NULL, algo=param$reg.algo
		);
	}
	return(obs.loglik+alpha.loglik+beta.loglik+gamma.loglik+u.loglik+v.loglik+w.loglik);
}
E.loglik.mainEffect <- function(
	coeff, slope, var,
	x, local.mean, global.mean, local.var, global.var, local.cov,
	subset=NULL, algo=NULL
){
	nContexts = if(is.vector(local.mean)) 1 else ncol(local.mean);
	nFeatures = ncol(x);
	if(nContexts == 1){
		if(!is.null(subset)){
			if(is.null(subset[[1]])) stop("is.null(subset[[1]])");
			x = x[subset[[1]],,drop=FALSE];
			local.mean = local.mean[subset[[1]]];
			local.var  = local.var[subset[[1]]];
		}
		
		if(is.list(coeff))        this.coeff = coeff[[1]]
		else if(is.vector(coeff)) this.coeff = coeff
		else                      this.coeff = coeff[,1];
		
		ans = Eloglik.lm.random.effect(
			coeff=this.coeff, var=var,
			response.mean = local.mean, 
			feature.mean  = x,
			response.var  = local.var, algo=algo
		);
		return(ans);
	}
	loglik = rep(NA, nContexts);
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
		z = cbind(global.mean_k, as.matrix(x_k));
		Delta = cbind(global.var_k, matrix(0, nrow=nrow(z), ncol=nFeatures));
		c = cbind(local.cov_k, matrix(0, nrow=nrow(z), ncol=nFeatures));
		eta = c(slope[k], coeff[,k]);
		loglik[k] = Eloglik.lm.random.effect(
				coeff=eta, var=var[k],
				response.mean = local.mean_k, 
				response.var  = local.var_k,
				feature.mean = z,
				feature.var  = Delta, 
				feature.cov  = c
		);
	}
	return(sum(loglik));
}
E.loglik.interaction <- function(
	coeff, var, x, factor.mean, factor.var, subset=NULL, algo=NULL
){
	nFeatures = ncol(x);
	nFactors  = ncol(factor.mean);
	nNodes    = nrow(x);
	if(nrow(factor.mean) != nNodes)   stop("nrow(x) != nrow(factor.mean)");
	if(ncol(factor.var)  != nFactors) stop("ncol(factor.var) != nFactors");

	if(length(var) == 1) var = rep(var, nFactors);
	if(length(var) != nFactors) stop("length(var) != nFactors");
	
	loglik = rep(NA, nFactors);
	
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
		
		if(is.list(coeff)) this.coeff = coeff[[f]]
		else               this.coeff = coeff[,f];
		
		loglik[f] = Eloglik.lm.random.effect(
				coeff=this.coeff, var=var[f],
				response.mean = f.mean, 
				response.var  = f.var, 
				feature.mean  = x.f, algo=algo
		);
	}
	return(sum(loglik));
}

###
### Complete data log likelihood with the constant term removed
###
###	    If is.logistic = TRUE, obs$response will be used to compute likelihood
###                            obs$y        will be used to compute rmse
###
### TODO: Summing over only the observed nodes
###
complete.logLikelihood <- function(
	obs, factor, feature, param,
	is.logistic=FALSE
){
	size = syncheck.multicontext.spec(factor=factor, obs=obs, feature=feature, param=param);
	
	nObs        = size$nObs;
	nSrcNodes   = size$nSrcNodes;
	nDstNodes   = size$nDstNodes;
	nFactors    = size$nFactors;
	nSrcContexts  = size$nSrcContexts;
	nDstContexts  = size$nDstContexts;
	nEdgeContexts = size$nEdgeContexts;
	
	pred.y = predict.y.from.factors(obs=obs, factor=factor, feature=feature, param=param);
	ans = obsLoglik.from.gaussian(pred.y=pred.y, y=obs$response, var_y=param$var_y, is.logistic=is.logistic);
	
	if(!is.null(param$var_alpha) && any(param$var_alpha > 0)){
		pred = reg.predict(model=param$g0, x=feature$x_src, algo=param$reg.algo, ncol=nSrcContexts);
		if(nSrcContexts == 1){
			ans = ans + loglik.gaussian(pred.x=pred, x=factor$alpha, var_x=param$var_alpha);
		}else if(nSrcContexts > 1){
			for(i in 1:nSrcContexts){
				ans = ans + loglik.gaussian(pred.x=pred[,i]+factor$alpha_global*param[["q"]][i], x=factor$alpha[,i], var_x=param$var_alpha[i]);
			}
			ans = ans + loglik.gaussian(pred.x=0, x=factor$alpha_global, var_x=param$var_alpha_global);
		}else stop("nSrcContexts = ", nSrcContexts);
	}
	if(!is.null(param$var_beta) && any(param$var_beta > 0)){
		pred = reg.predict(model=param$d0, x=feature$x_dst, algo=param$reg.algo, ncol=nDstContexts);
		if(nDstContexts == 1){
			ans = ans + loglik.gaussian(pred.x=pred, x=factor$beta, var_x=param$var_beta);
		}else if(nDstContexts > 1){
			for(i in 1:nDstContexts){
				ans = ans + loglik.gaussian(pred.x=pred[,i]+factor$beta_global*param[["r"]][i], x=factor$beta[,i], var_x=param$var_beta[i]);
			}
			ans = ans + loglik.gaussian(pred.x=0, x=factor$beta_global, var_x=param$var_beta_global);
		}else stop("nDstContexts = ", nDstContexts);
	}
	if(!is.null(param$var_gamma) && param$var_gamma > 0){
		pred = reg.predict(model=param$h0, x=feature$x_ctx, algo=param$reg.algo);
		ans = ans + loglik.gaussian(pred.x=pred, x=factor$gamma, var_x=param$var_gamma);
	}
	if(!is.null(param$var_u) && any(param$var_u > 0)){
		pred = reg.predict(model=param$G, x=feature$x_src, algo=param$reg.algo, ncol=nFactors);
		if(length(param$var_u) == 1){
			ans = ans + loglik.gaussian(pred.x=pred, x=factor$u, var_x=param$var_u);
		}else if(length(param$var_u) == ncol(pred)){
			for(k in 1:ncol(pred)) ans = ans + loglik.gaussian(pred.x=pred[,k], x=factor$u[,k], var_x=param$var_u[k]);
		}else stop("length(param$var_u) = ", length(param$var_u));
	}
	if(!is.null(param$var_v) && any(param$var_v > 0)){
		pred = reg.predict(model=param$D, x=feature$x_dst, algo=param$reg.algo, ncol=size$nFactors);
		if(length(param$var_v) == 1){
			ans = ans + loglik.gaussian(pred.x=pred, x=factor$v, var_x=param$var_v);
		}else if(length(param$var_v) == ncol(pred)){
			for(k in 1:ncol(pred)) ans = ans + loglik.gaussian(pred.x=pred[,k], x=factor$v[,k], var_x=param$var_v[k]);
		}else stop("length(param$var_v) = ", length(param$var_v));
	}
	if(!is.null(param$var_w) && any(param$var_w > 0)){
		pred = reg.predict(model=param$H, x=feature$x_ctx, algo=param$reg.algo, ncol=size$nFactors);
		ans = ans + loglik.gaussian(pred.x=pred, x=factor$w, var_x=param$var_w);
	}
	
	return(ans);
}

#####################################################################
### Backup
#####################################################################

###
### Simple random init (NOT UP TO DATE!! Check the ..._EM.R file)
###
###	  alpha, beta, gamma, u, v, w are generated by N(mean=0, sd=factor.sd)
###   b is fitted without any factor
###   g0, d0, h0, G, D, H are all 0 (or NA)
###   q and r are all 1
###   All variances are (factor.sd * sqrt(var_y))^2
###
init.simple_random <- function(
		obs,          # Observations
		feature,      # Features
		interaction.dim,   # Number of interaction factors per entity
		has.gamma, reg.algo=NULL, reg.control=NULL, has.u=TRUE,
		factor.sd=0.2, is.logistic=FALSE
){
	nObs = nrow(obs); nSrcNodes = nrow(feature$x_src); nDstNodes = nrow(feature$x_dst);
	nSrcFeatures = ncol(feature$x_src);  nDstFeatures = ncol(feature$x_dst);
	
	nSrcContexts  = if(!is.null(obs$src.context))  max(obs$src.context)  else 1;
	nDstContexts  = if(!is.null(obs$dst.context))  max(obs$dst.context)  else 1;
	nEdgeContexts = if(!is.null(obs$edge.context)) max(obs$edge.context) else 1;
	
	param = list(reg.algo=reg.algo, reg.control=reg.algo$control, b=NA, 
			g0=as.list(rep(NA,nSrcContexts)), d0=as.list(rep(NA,nDstContexts)));
	if(!is.null(reg.control)){
		for(i in 1:length(reg.control)){
			param$reg.control[[names(reg.control)[i]]] = reg.control[[i]];
		}
	}
	if(is.logistic) stop("Logistic not yet supported");
	if(is.null(reg.algo)){
		temp = lm(obs$y ~ feature$x_obs, model=FALSE);
		param$b = temp$coefficients;
	}else{
		param$b = param$reg.algo$train(x=feature$x_obs, y=obs$y, control=param$reg.control);
	}
	y.pred  = reg.predict(model=param$b, x=feature$x_obs, algo=param$reg.algo);
	var_y   = mean((obs$y - y.pred)^2);
	factor.sd = factor.sd * sqrt(var_y);
	init.var  = factor.sd^2;
	
	factor = list(
			alpha=matrix(rnorm(nSrcNodes*nSrcContexts,sd=factor.sd), nrow=nSrcNodes),
			beta =matrix(rnorm(nDstNodes*nDstContexts,sd=factor.sd), nrow=nDstNodes)
	);
	if(has.gamma && nEdgeContexts >= 2){
		factor$gamma = rnorm(nEdgeContexts,sd=factor.sd);
		param$h0     = NA;
	}
	if(interaction.dim >= 1){
		if(has.u){
			factor$u = matrix(rnorm(nSrcNodes*interaction.dim,sd=factor.sd),nrow=nSrcNodes);
			param$G  = as.list(rep(NA,interaction.dim));
		}
		factor$v = matrix(rnorm(nDstNodes*interaction.dim,sd=factor.sd),nrow=nDstNodes);
		param$D  = as.list(rep(NA,interaction.dim));
		if(nEdgeContexts >= 2){
			factor$w = matrix(rnorm(nEdgeContexts*interaction.dim,sd=factor.sd),nrow=nEdgeContexts);
			param$H = as.list(rep(NA,interaction.dim));
		}
	}
	
	param[["q"]] = rep(1, nSrcContexts);
	param[["r"]] = rep(1, nDstContexts);
	
	param$var_y     = var_y;
	param$var_alpha = rep(init.var, nSrcContexts);
	if(nSrcContexts >= 2) param$var_alpha_global = init.var;
	param$var_beta  = rep(init.var, nDstContexts);
	if(nDstContexts >= 2) param$var_beta_global = init.var;
	if(!is.null(factor$gamma)) param$var_gamma = init.var;
	if(!is.null(factor$u))     param$var_u     = init.var;
	if(!is.null(factor$v))     param$var_v     = init.var;
	if(!is.null(factor$w))     param$var_w     = init.var;
	
	ans = list(factor=factor, param=param);
}
