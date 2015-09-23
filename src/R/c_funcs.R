### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
###
### Author: Bee-Chung Chen

print.address <- function(x){
	if(is.double(x)){
		.Call("print_doublePointer_Call", x);
	}else if(is.integer(x)){
		.Call("print_intPointer_Call", x);
	}else stop("x is not double or integer")
}

###
### out[m] = sum(u[src.id[m], , src.ctx[m]] * v[dst.id[m, , dst.ctx[m]])
###
compute_multicontext_uv <- function(
	u, v, src.id, src.ctx, dst.id, dst.ctx, debug=0
){
	nObs = length(src.id);  nFactors  = dim(u)[2];
	nSrcNodes = dim(u)[1];  nDstNodes = dim(v)[1];
	nSrcCtx   = dim(u)[3];  nDstCtx   = dim(v)[3];

	if(length(dst.id)  != nObs) stop("length(dst.id)  != nObs");
	if(length(src.ctx) != nObs) stop("length(src.ctx) != nObs");
	if(length(dst.ctx) != nObs) stop("length(dst.ctx) != nObs");
	if(dim(v)[2] != nFactors) stop("dim(v)[2] != nFactors");
	if(!is.double(u)) stop("!is.double(u)");
	if(!is.double(v)) stop("!is.double(v)");
	if(!is.integer(src.id))  stop("!is.integer(src.id)");
	if(!is.integer(src.ctx)) stop("!is.integer(src.ctx)");
	if(!is.integer(dst.id))  stop("!is.integer(dst.id)");
	if(!is.integer(dst.ctx)) stop("!is.integer(dst.ctx)");

	out = double(nObs) + 0;

	.Call("computeMultiResponseUV_Call",
		# OUTPUT
		out,
		# INPUT
		u, v, src.id,  dst.id, src.ctx, dst.ctx,
		as.integer(nObs), as.integer(nFactors),
		as.integer(nSrcNodes), as.integer(nSrcCtx),
		as.integer(nDstNodes), as.integer(nDstCtx),
		as.integer(debug)
	);

	return(out);
}

###
### out[k] = t(b) %*% A[k,,] %*% b
###
compute_bAb_3DA <- function(A,b){
	if(length(dim(A)) != 3) stop("length(dim(A)) != 3");
	nCases = as.integer(dim(A)[1]);
	nDim   = as.integer(length(b));
	if(any(dim(A)[2:3] != nDim)) stop("any(dim(A)[2:3] != nDim)");
	if(!is.double(A)) stop("!is.double(A)");
	if(!is.double(b)) stop("!is.double(b)");
	out = double(nCases);
	.Call("compute_bAb_3DA_Call", out, b, A, nCases, nDim);
	return(out + 0);
}

###
### out[k,] = A[k,,] %*% b
###
compute_Ab_3DA <- function(A,b){
	if(length(dim(A)) != 3) stop("length(dim(A)) != 3");
	nCases = as.integer(dim(A)[1]);
	nDim   = as.integer(length(b));
	if(any(dim(A)[2:3] != nDim)) stop("any(dim(A)[2:3] != nDim)");
	if(!is.double(A)) stop("!is.double(A)");
	if(!is.double(b)) stop("!is.double(b)");
	out = matrix(double(1), nrow=nCases, ncol=nDim);
	.Call("compute_Ab_3DA_Call", out, b, A, nCases, nDim);
	return(out + 0);
}

###
### Inversion of a symmetric metrix using Cholesky
###     This function is implemented to make sure the C implementation generates exactly
###     the same output as the R implementation
###
sym_inv.cholesky <- function(x){
    if(!is.double(x)) stop("x is not a double matrix");
    if(!is.matrix(x)) stop("x is not a double matrix");
    if(nrow(x) != ncol(x)) stop("x is not a square matrix");
    n = as.integer(nrow(x));
    ans = .C("sym_inv_byCholesky", x, n, as.integer(1), DUP=TRUE);
    return(ans[[1]]);
}

###
### Eigen decomposition of a symmetric metrix
###     This function is implemented to make sure the C implementation generates exactly
###     the same output as the R implementation
###
sym_eigen <- function(x){
    if(!is.matrix(x)) stop("x is not matrix");
    if(nrow(x) != ncol(x)) stop("nrow(x) != ncol(x)");
    if(!is.double(x)) stop("x is not a double-precision matrix");
    output = list(values  = rep(0.0, nrow(x)),
                  vectors = matrix(0.0, nrow=nrow(x), ncol=ncol(x)));
    .Call("sym_eigen_Call", x, nrow(x), output$values, output$vectors);
    return(output);
}

###
### Sum up each row (or column) of a matrix
###     side=1: Sum up each row and return a vector with length nrow
###     side=2: Sum up each column and return a vector with length ncol
###
sum_margin <- function(A, side){
    if(!is.matrix(A)) stop("'A' is not a matrix");
    if(!is.double(A)) stop("'A' is not double");
    if(side == 1) out = double(nrow(A))
    else if(side == 2) out = double(ncol(A))
    else stop("Unknown side=",side," (side=1 for rows, side=2 for columns)");
    .Call("sum_margin_Call", out, A, as.integer(nrow(A)), as.integer(ncol(A)), as.integer(side));
    return(out);
}

###
### out[k,j] = sum_{x s.t. groupBy[x] = j} mat[k,select[x]] * weight[x]
###
selectColumn_agg_sum <- function(mat, groupBy, select, weight=NULL){
    out = matrix(0.0, nrow=nrow(mat), ncol=max(groupBy));
    if(!is.double(mat)) stop("!is.double(mat)");
    if(!is.integer(groupBy)) stop("!is.integer(groupBy)");
    if(!is.integer(select)) stop("!is.integer(select)");
    if(length(select) != length(groupBy)) stop("length(select) != length(groupBy)");
    nWeightw = 0;
    if(!is.null(weight)){
        nWeights = length(weight);
        if(nWeights != length(select)) stop("length(weight) != length(select)");
        if(!is.double(weight)) stop("!is.double(weight)");
    }
    .Call("selectColumn_agg_sum_Call",
        out, as.integer(nrow(out)), as.integer(ncol(out)),
        mat, as.integer(nrow(mat)), as.integer(ncol(mat)),
        select, groupBy, as.integer(length(select)),
        weight, as.integer(nWeightw)
    );
    return(out);
}

###
### margin = 1: each row sum up to one
### margin = 2: each column sum up to one
###
normalize_sumToOne2D <- function(mat, margin){
    out = matrix(0.0, nrow=nrow(mat), ncol=ncol(mat));
    if(!is.double(mat)) stop("!is.double(mat)");
    .Call("normalize_sumToOne2D_Call",
        out, mat, as.integer(nrow(mat)),as.integer(ncol(mat)), as.integer(margin)
    );
    return(out);
}

###
### out[i] = sum_{p,q} u[u.index[i],p] * B[p,q] * v[v.index,q]
###
compute.uBv <- function(u.index, v.index, u, B, v){

	if(!is.matrix(u)) stop("!is.matrix(u)");
	if(!is.matrix(v)) stop("!is.matrix(v)");
	if(!is.matrix(B)) stop("!is.matrix(B)");
	if(length(v.index) != length(u.index)) stop("length(v.index) != length(u.index)");

	nObs = as.integer(length(u.index));
	nUsers = nrow(u);
	nFactorsPerUser = ncol(u);
	nItems = nrow(v);
	nFactorsPerItem = ncol(v);

	score = double(length=nObs);

	check_type_size(score,"double",nObs);
	check_type_size(u.index,"integer",nObs);
	check_type_size(v.index,"integer",nObs);
	check_type_size(u,"double",c(nUsers, nFactorsPerUser));
	check_type_size(v,"double",c(nItems, nFactorsPerItem));
	check_type_size(B,"double",c(nFactorsPerUser, nFactorsPerItem));

	.Call("compute_uBv_dense_Call",
			# Output
			score,
			# Input
			u.index, v.index, u, B, v,
			nObs, nUsers, nItems, nFactorsPerUser, nFactorsPerItem
	);
	return(score);
}

###
### See condMeanVarSample_singleDim in MCEM_EStep.c
###
condMeanVarSample_singleDim.C <- function(
    option, # 1:Sample, 2:Mean&Var, 3:Sample&Mean&Var
    thisEffIndex, rest, # rest = the o in the paper
    fittedEff, multiplier, var_y, var_eff, debug=0
){
    nObs = as.integer(length(thisEffIndex));
    nThisEff  = as.integer(length(fittedEff));
    nVar_y    = as.integer(length(var_y));
    nVar_eff  = as.integer(length(var_eff));

    if(length(rest) != nObs) stop("length(rest) != nObs");
    if(!is.null(multiplier) && length(multiplier) != nObs) stop("length(multiplier) != nObs");
    if(!(nVar_y == 1 || nVar_y == nObs)) stop("length(var_y) has problem");
    if(!(nVar_eff == 1 || nVar_eff == nThisEff)) stop("length(var_eff) has problem");

    out = list(sample=as.double(NULL), mean=as.double(NULL), var=as.double(NULL));
    if(option == 1 || option == 3)  out$sample = rep(0.0, nThisEff);
    if(option == 2 || option == 3){ out$mean   = rep(0.0, nThisEff);   out$var = rep(0.0, nThisEff);}

    if(!is.double(rest)) stop("!is.double(rest)");
    if(!is.null(multiplier) && !is.double(multiplier)) stop("!is.double(multiplier)");
    if(!is.double(fittedEff)) stop("!is.double(fittedEff)");
    if(!is.double(var_y)) stop("!is.double(var_y)");
    if(!is.double(var_eff)) stop("!is.double(var_eff)");

	check_type_size(out$sample,"double",nThisEff);
	check_type_size(out$mean,"double",nThisEff,isNullOK=(option==1));
	check_type_size(out$var,"double",nThisEff,isNullOK=(option==1));

	ans = .C("condMeanVarSample_singleDim",
        # OUTPUT
        out$sample, out$mean, out$var,
        # INPUT
        as.integer(option), as.integer(thisEffIndex), rest,
        fittedEff, multiplier, var_y, var_eff,
        nObs, nThisEff, nVar_y, nVar_eff,
        as.integer(debug), DUP=FALSE
    );

    return(out);
}

###
### See generateObsIndex in util.c
###
generateObsIndex <- function(
    effIndex, nEff, debug
){
    nEff = as.integer(nEff);
    nObs   = as.integer(length(effIndex));

    out = list(
        obsIndex = integer(nObs),
        start    = integer(nEff),
        num      = integer(nEff)
    );

    if(!is.integer(effIndex)) stop("!is.integer(effIndex)");

    .Call("generateObsIndex_Call",
        out$obsIndex, out$start, out$num,
        # INPUT
        effIndex, nObs, nEff,
        # OTHER
        as.integer(debug)
    );
    return(out);
}

###
### See condMeanVarSample_multiDim in MCEM_EStep.c
###
condMeanVarSample_multiDim.C <- function(
    option, # 1:Sample, 2:Mean&Var, 3:Sample&Mean&Var
    thisEffIndex, otherEffIndex, rest, # rest = the o in the paper
    fittedEff, otherEff, thisCluster, otherCluster, B,
	var_y, var_eff, oi=NULL, debug=0, verbose=0
){
    nObs = as.integer(length(thisEffIndex));
    nThisEff  = as.integer(nrow(fittedEff));
    nOtherEff = as.integer(nrow(otherEff));
	nFactorsPerThis  = as.integer(ncol(fittedEff));
	nFactorsPerOther = as.integer(ncol(otherEff));

    nVar_y = as.integer(length(var_y));
    nVar_eff = as.integer(length(var_eff));
    if(nVar_eff > 1){
        nVar_eff = as.integer(dim(var_eff)[1]);
    }

    if(length(otherEffIndex) != nObs) stop("length(otherEffIndex) != nObs");
    if(length(rest) != nObs) stop("length(rest) != nObs");
    if(length(thisCluster)  != nrow(fittedEff)) stop("length(thisCluster) != nrow(fittedEff)");
	if(length(otherCluster) != nrow(otherEff))  stop("length(otherCluster) != nrow(otherEff)");

	if(length(B) == 1 && B == 1){
		identity_B = as.integer(1);
		nThisClusters = as.integer(0);
		nOtherClusters = as.integer(0);
	}else{
		identity_B = as.integer(0);
		nThisClusters = as.integer(dim(B)[1]);
		nOtherClusters = as.integer(dim(B)[2]);
		if(dim(B)[1] < max(thisCluster))  stop("dim(B)[1] < max(thisCluster)");
		if(dim(B)[2] < max(otherCluster)) stop("dim(B)[2] < max(otherCluster)");
		if(dim(B)[3] != nFactorsPerThis)  stop("dim(B)[3] != nFactorsPerThis");
		if(dim(B)[4] != nFactorsPerOther) stop("dim(B)[4] != nFactorsPerOther");
	}

    if(!(nVar_y == 1 || nVar_y == nObs)) stop("length(var_y) has problem");
    if(!(nVar_eff == 1 || nVar_eff == nThisEff)) stop("length(var_eff) has problem");

    out = list(sample=as.double(NULL), mean=as.double(NULL), var=as.double(NULL));
    if(option == 1 || option == 3)  out$sample = matrix(0.0, nrow=nThisEff, ncol=nFactorsPerThis);
    if(option == 2 || option == 3){ out$mean   = matrix(0.0, nrow=nThisEff, ncol=nFactorsPerThis);  out$var = array(0.0, dim=c(nThisEff, nFactorsPerThis, nFactorsPerThis));}

    if(!is.double(rest)) stop("!is.double(rest)");
    if(!is.double(fittedEff)) stop("!is.double(fittedEff)");
    if(!is.double(otherEff)) stop("!is.double(otherEff)");
    if(!is.double(var_y)) stop("!is.double(var_y)");
    if(!is.double(var_eff)) stop("!is.double(var_eff)");

    if(is.null(oi)){
        cat("obsIndex is null; built it now!\n");
        oi = generateObsIndex(thisEffIndex, nThisEff, debug);
    }

	thisCluster  = thisCluster  - as.integer(1);
	otherCluster = otherCluster - as.integer(1);

    if(length(oi$obsIndex) != nObs) stop("length(oi$obsIndex) != nObs");
    if(length(oi$start) != nThisEff) stop("length(oi$start) != nThisEff");
    if(length(oi$num) != nThisEff) stop("length(oi$num) != nThisEff");
    if(!is.integer(oi$obsIndex)) stop("!is.integer(oi$obsIndex)");
    if(!is.integer(oi$start)) stop("!is.integer(oi$start)");
    if(!is.integer(oi$num)) stop("!is.integer(oi$num)");
	if(!is.integer(thisCluster)) stop("!is.integer(thisCluster)");
	if(!is.integer(otherCluster)) stop("!is.integer(otherCluster)");
	if(!is.double(B)) stop("!is.double(B)");

	check_type_size(out$sample,"double",c(nThisEff,nFactorsPerThis));
	check_type_size(out$mean,"double",c(nThisEff,nFactorsPerThis),isNullOK=(option==1));
	check_type_size(out$var,"double",c(nThisEff,nFactorsPerThis,nFactorsPerThis),isNullOK=(option==1));

	ans = .C("condMeanVarSample_multiDim",
        # OUTPUT
        out$sample, out$mean, out$var,
        # INPUT
        as.integer(option), as.integer(thisEffIndex), as.integer(otherEffIndex), rest,
        fittedEff, otherEff, thisCluster, otherCluster, B,
		var_y, var_eff,
        nObs, nThisEff, nOtherEff, nFactorsPerThis, nFactorsPerOther, nVar_y, nVar_eff,
		nThisClusters, nOtherClusters,
        oi$obsIndex, oi$start, oi$num, identity_B,
        as.integer(debug), as.integer(verbose),
        DUP=FALSE
    );

	return(out);
}

###
### See condMeanVarSample_multiDim in MCEM_EStep.c
###
condMeanVarSample_cluster.C <- function(
	option, # 1:Sample, 2:Probabilities, 3:Sample&Probabilities
	logPrior,
	thisEffIndex, otherEffIndex, rest, # rest = the o in the paper
	thisEff, otherEff, otherCluster, B, c,
	var_y, oi=NULL, debug=0, verbose=0
){
	nObs = as.integer(length(thisEffIndex));
	nThisEff  = as.integer(nrow(thisEff));
	nOtherEff = as.integer(nrow(otherEff));
	nFactorsPerThis  = as.integer(ncol(thisEff));
	nFactorsPerOther = as.integer(ncol(otherEff));

	nThisClusters  = as.integer(dim(logPrior)[2]);
	nOtherClusters = 0;

	if(length(c) == 1 && c==0){
		has_c = as.integer(0);
	}else{
		has_c = as.integer(1);
		nOtherClusters = as.integer(dim(c)[2]);
		check_type_size(c, "double", c(nThisClusters,nOtherClusters));
	}

	if(length(B) == 1 && B == 1){
		identity_B = as.integer(1);
	}else{
		identity_B = as.integer(0);
		nOtherClusters = as.integer(dim(B)[2]);
		check_type_size(B, "double", c(nThisClusters,nOtherClusters,nFactorsPerThis,nFactorsPerOther));
	}

	if(nOtherClusters > 0 && nOtherClusters < max(otherCluster)) stop("nOtherClusters < max(otherCluster)");

	nVar_y = as.integer(length(var_y));
	otherCluster = otherCluster - as.integer(1);

	check_type_size(thisEff, "double", c(nThisEff, nFactorsPerThis));
	check_type_size(otherEff, "double", c(nOtherEff, nFactorsPerOther));
	check_type_size(otherEffIndex, "int", nObs);
	check_type_size(rest, "double", nObs);
	check_type_size(otherCluster, "int", nOtherEff);
	check_type_size(logPrior, "double", c(nThisEff, nThisClusters));

	if(!(nVar_y == 1 || nVar_y == nObs)) stop("length(var_y) has problem");

	out = list();
	if(option == 1 || option == 3) temp = integer(nThisEff);
	if(option == 2 || option == 3) out$prob   = matrix(0.0, nrow=nThisEff, ncol=nThisClusters);

	if(!is.double(var_y)) stop("!is.double(var_y)");

	if(is.null(oi)){
		cat("obsIndex is null; built it now!\n");
		oi = generateObsIndex(thisEffIndex, nThisEff, debug);
	}

	if(length(oi$obsIndex) != nObs) stop("length(oi$obsIndex) != nObs");
	if(length(oi$start) != nThisEff) stop("length(oi$start) != nThisEff");
	if(length(oi$num) != nThisEff) stop("length(oi$num) != nThisEff");
	if(!is.integer(oi$obsIndex)) stop("!is.integer(oi$obsIndex)");
	if(!is.integer(oi$start)) stop("!is.integer(oi$start)");
	if(!is.integer(oi$num)) stop("!is.integer(oi$num)");

	check_type_size(temp,"integer",nThisEff);
	check_type_size(out$prob,"double",c(nThisEff,nThisClusters),isNullOK=(option==1));

	ans = .C("condProbSample_cluster",
			# OUTPUT
			temp, out$prob,
			# INPUT
			as.integer(option), logPrior,
			rest, otherCluster,
			thisEff, otherEff, B, c,
			var_y, as.integer(thisEffIndex), as.integer(otherEffIndex),
			nThisEff, nOtherEff, nObs, nVar_y,
			nThisClusters, nOtherClusters,
			nFactorsPerThis, nFactorsPerOther,
			oi$obsIndex, oi$start, oi$num, identity_B, has_c,
			as.integer(debug), as.integer(verbose),
			DUP=FALSE
	);

	# print(temp);
	out$sample = temp + as.integer(1);

	return(out);
}


###
### MCEM_EStep.C. See MCEM_EStep(...) in MCEM_EStep.cpp
###
### INPUT:  factor  = list(alpha, beta, u, v, s, z); # Initial factor values
###         obs     = data.frame(y, user, item);
###         feature = list(x_dyad, x_user, x_item);
###         param   = list(b, g0, d0, G, D, s_model, z_model,
###                        var_y, var_alpha, var_beta, var_u, var_v, type, B, c);
###
### OUTPUT: mean     = list(alpha, beta, u, v, s, z);
###         sumvar   = list(alpha, beta, u, v, mu);
###         sampvar  = list(alpha, beta, u, v);
###
###   isOldUser[i]: [TRUE/FALSE] Whehter the ith user is an old user; default: all FALSE
###                 For old users, we set g0x_user[i] = alpha[i], Gx_user[i,,] = u[i,,]
###   isOldItem[j]: [TRUE/FALSE] Whehter the jth item is an old item; default: all FALSE
###                 For old items, we set d0x_item[j] = beta[j],  Dx_item[j,,] = v[j,,]
###
### NOTE: factor$alpha and factor$beta cannot be NULL!!
###
MCEM_EStep.C <- function(
	factor, obs, feature, param, nSamples, nBurnIn=1,
	s_algo=NULL, z_algo=NULL,
	userFactorVar=0, itemFactorVar=0,
	isOldUser=NULL, isOldItem=NULL,
	is.cluster_algo_NULL.ok=FALSE,
	debug=0, verbose=0
){
	size = syncheck.cRLFM.spec(factor=factor, obs=obs, feature=feature, param=param);

	if(param$type == "c"){
		stop("type == 'c' is not yet supported");
	}else if(param$type == "q"){
	}else stop("param$type must be either 'c' or 'q'");

	if(size$nUserClusters > 1 && ((!is.cluster_algo_NULL.ok && (is.null(s_algo)||is.null(param$s_model))) || is.null(factor$s))) stop("Please specify s_algo, param$s_model and factor$s!");
	if(size$nItemClusters > 1 && ((!is.cluster_algo_NULL.ok && (is.null(z_algo)||is.null(param$z_model))) || is.null(factor$z))) stop("Please specify z_algo, param$z_model and factor$z!");
	if(is.null(factor$alpha)) stop("alpha cannot be null");
	if(is.null(factor$beta))  stop("beta cannot be null");
	if(length(obs) == 0) stop("No observation data");

    sumvar  = list(alpha=double(1), beta=double(1), u=double(1), v=double(1), mu=double(1));
    sampvar = list();
    if(userFactorVar != 0){
        sampvar$alpha = double(size$nUsers);
        sampvar$u     = array(double(1),dim=c(size$nUsers, size$nFactorsPerUser, size$nFactorsPerUser));
    }
    if(itemFactorVar != 0){
        sampvar$beta = double(size$nItems);
        sampvar$v    = array(double(1),dim=c(size$nItems, size$nFactorsPerItem, size$nFactorsPerItem));
    }

	xb  = feature$x_dyad %*% param$b;
	g0x_user = NULL; Gx_user = NULL;
    if(is.null(isOldUser)){
        if(!is.null(param$g0)) g0x_user = feature$x_user %*% param$g0;
        if(!is.null(param$G))   Gx_user = feature$x_user %*% param$G;
		if(!is.null(s_algo) && !is.null(param$s_model))
			 s_logPrior = log(s_algo$predict(param$s_model, feature$x_user))
		else s_logPrior = matrix(0.0, nrow=size$nUsers, ncol=size$nUserClusters);
	}else{
        if(length(isOldUser) != size$nUsers) stop("length(isOldUser) != nUsers");
        x_user.new = feature$x_user[!isOldUser,,drop=FALSE];
        if(!is.null(param$g0)){ g0x_user = factor$alpha;  g0x_user[!isOldUser] = x_user.new %*% param$g0;}
        if(!is.null(param$G)){   Gx_user = factor$u;      Gx_user[!isOldUser,] = x_user.new %*% param$G;}
		s_logPrior = log(factor$s);
		if(!is.null(s_algo) && !is.null(param$s_model))
			 s_logPrior[!isOldUser,] = log(s_algo$predict(param$s_model, x_user.new))
		else s_logPrior[!isOldUser,] = 0.0;
	}

    d0x_item = NULL;  Dx_item = NULL;
    if(is.null(isOldItem)){
        if(!is.null(param$d0)) d0x_item = feature$x_item %*% param$d0;
        if(!is.null(param$D))   Dx_item = feature$x_item %*% param$D;
		if(!is.null(z_algo) && !is.null(param$z_model))
			 z_logPrior = log(z_algo$predict(param$z_model, feature$x_item))
		else z_logPrior = matrix(0.0, nrow=size$nItems, ncol=size$nItemClusters);
	}else{
        if(length(isOldItem) != size$nItems) stop("length(isOldItem) != nItems");
        x_item.new = feature$x_item[!isOldItem,,drop=FALSE];
        if(!is.null(param$d0)){ d0x_item = factor$beta;   d0x_item[!isOldItem] = x_item.new %*% d0;}
        if(!is.null(param$D)){   Dx_item = factor$v;      Dx_item[!isOldItem,] = x_item.new %*% D;}
		z_logPrior = log(factor$z);
		if(!is.null(z_algo) && !is.null(param$z_model))
			 z_logPrior[!isOldItem,] = log(z_algo$predict(param$z_model, x_item.new))
		else z_logPrior[!isOldItem,] = 0.0;
    }

    # Checks
    if(size$nVar_alpha > 0 && (!is.double(g0x_user) || length(g0x_user) != size$nUsers)) stop("g0x_user!");
	if(size$nVar_beta  > 0 && (!is.double(d0x_item) || length(d0x_item) != size$nItems)) stop("d0x_item!");
	if(size$nVar_u     > 0 && (!is.double(Gx_user)  || any(dim(Gx_user) != c(size$nUsers,size$nFactorsPerUser)))) stop("Gx_user!");
    if(size$nVar_v     > 0 && (!is.double(Dx_item)  || any(dim(Dx_item) != c(size$nItems,size$nFactorsPerItem)))) stop("Dx_item!");
	if(size$nUserClusters > 1 && (!is.double(s_logPrior) || any(dim(s_logPrior) != c(size$nUsers, size$nUserClusters)))) stop("s_logPrior");
	if(size$nItemClusters > 1 && (!is.double(z_logPrior) || any(dim(z_logPrior) != c(size$nItems, size$nItemClusters)))) stop("z_logPrior");

    problem.dim = c(size$nObs, size$nUsers, size$nItems, size$nFactorsPerUser, size$nFactorsPerItem, size$nUserClusters, size$nItemClusters,
			        size$nVar_y, size$nVar_alpha, size$nVar_beta, size$nVar_u, size$nVar_v, length(param$B), length(param$c));
    if(!is.integer(problem.dim)) stop("!is.integer(problem.dim)");

	if(!is.double(xb) || length(xb) != size$nObs) stop("xb");

	check_type_size(factor$alpha, "double", problem.dim[2]);
	check_type_size(factor$beta, "double", problem.dim[3]);
	check_type_size(factor$u, "double", problem.dim[c(2,4)]);
	check_type_size(factor$v, "double", problem.dim[c(3,5)]);
	check_type_size(factor$s, "double", problem.dim[c(2,6)]);
	check_type_size(factor$z, "double", problem.dim[c(3,7)]);
	check_type_size(sumvar$alpha, "double", 1);
	check_type_size(sumvar$beta, "double", 1);
	check_type_size(sumvar$u, "double", 1);
	check_type_size(sumvar$v, "double", 1);
	check_type_size(sumvar$mu, "double", 1);
	check_type_size(sampvar$alpha, "double", problem.dim[2], isNullOK=(userFactorVar==0));
	check_type_size(sampvar$beta, "double", problem.dim[3], isNullOK=(itemFactorVar==0));
	check_type_size(sampvar$u, "double", problem.dim[c(2,4,4)], isNullOK=(userFactorVar==0));
	check_type_size(sampvar$v, "double", problem.dim[c(3,5,5)], isNullOK=(itemFactorVar==0));

	#	//  dim = {1:nObs, 2:nUsers, 3:nItems, 4:nFactorsPerUser, 5:nFactorsPerItem, 6:nUserClusters,
	#   //         7:nItemClusters, 8:nVar_y, 9:nVar_alpha, 10:nVar_beta, 11:nVar_u, 12:nVar_v}

	ans = .C("MCEM_EStep",
		# INPUT (initial factor values) & OUTPUT (Monte Carlo mean of factor values)
		factor$alpha, factor$beta, factor$u, factor$v, factor$s, factor$z,
		# OUTPUT
		sumvar$alpha,  	sampvar$alpha,
		sumvar$beta,   	sampvar$beta,
		sumvar$u, 		sampvar$u,
		sumvar$v, 		sampvar$v,
		sumvar$mu,
		# INPUT
		as.integer(nSamples), as.integer(nBurnIn),
		obs$user, obs$item, obs$y,
		xb, g0x_user, d0x_item, Gx_user, Dx_item, param$B, param$c, s_logPrior, z_logPrior,
		param$var_y, param$var_alpha, param$var_beta, param$var_u, param$var_v,
		problem.dim, as.integer(length(problem.dim)),
	    as.integer(userFactorVar), as.integer(itemFactorVar),
		# OTHER
		as.integer(debug),  as.integer(verbose),
		DUP=FALSE
    );
    output = list(
        mean=factor, sumvar=sumvar
    );
	if(userFactorVar != 0 || userFactorVar != 0) output$sampvar=sampvar;
    return(output);
}

indexWithQuantities <- function(x){
    len = sum(x);
    out = integer(len);
    .Call("indexWithQuantities_Call",
        out, as.integer(x), as.integer(length(x))
    );
    return(out);
}

###
### Batched per-item online gaussian regression & prediction
###
###		obs$y[k] = offset[k] + t(u[obs$user[k],]) * v[obs$item[k],]
###
### Output$prediction = data.frame(y, user, item, predicted.y, batch)
###		   * nrow may be fewer because nBatches * nObsPerBatch sets a threshold
###        * Ordered by item
###
predict.perItem_OnlineBatchGaussianRegression <- function(
	obs, # data.frame(y, user, item)
	u,   # nUsers x nFactors
	offset, # length: nObs
	beta.prior_mean, # length: nItems
	beta.prior_var,  # length: nItems
	v.prior_mean,    # nItems x nFactors
	v.prior_var,     # length: nItems, or nItems x nFactors x nFactors
	discount.factor, # Downweight old batches exponentially by this number
	nBatches, nObsPerBatch,
	output.factors=FALSE, debug=0, verbose=0
){
	obs$temp = 1:nrow(obs);
	obs$offset = offset;
	obs = obs[order(obs$item,obs$temp),];
	obs$temp = NULL;

	nObs = as.integer(nrow(obs));
	nUsers = as.integer(nrow(u));
	nItems = as.integer(nrow(v.prior_mean));
	nFactors = as.integer(ncol(u));
	nBatches = as.integer(nBatches);
	nObsPerBatch = as.integer(nObsPerBatch);

	prediction = double(nrow(obs));
	batch_id   = integer(nrow(obs));
	beta = NULL;
	v    = NULL;
	if(output.factors){
		beta = array(double(1), dim=c(nItems, nBatches));
		v    = array(double(1), dim=c(nItems, nBatches, nFactors));
	}
	check_type_size(obs$y, "double", nObs);
	check_type_size(obs$user, "int", nObs);
	check_type_size(obs$item, "int", nObs);
	check_type_size(u, "double", c(nUsers, nFactors));
	check_type_size(obs$offset, "double", nObs);
	check_type_size(beta.prior_mean, "double", nItems);
	check_type_size(beta.prior_var, "double", nItems);
	check_type_size(v.prior_mean, "double", c(nItems, nFactors));
	if(!is.double(v.prior_var)) stop("v.prior_var should have type 'double'");
	if(!( (is.vector(v.prior_var) && length(v.prior_var)==nItems) ||
		  (length(dim(v.prior_var))==3 && all(dim(v.prior_var)==c(nItems, nFactors, nFactors))) ))
  		stop("Dimensionality of v.prior_var is not correct!!");

	.Call("perItem_online_factor_batch_predict_Call",
		prediction, # nObs x 1 (output)
		batch_id,   # nObs x 1 (output, ID starts from 1)
		beta,       # nItems x nBatches (output)
		v,          # nItems x nBatches x nFactors (output)
		obs$y,      # nObs x 1 (response)
		obs$user,   # nObs x 1 (index start from 1)
		obs$item,   # nObs x 1 (index start from 1)
		u,          # nUsers x nFactors
		obs$offset,     # nObs x 1
		beta.prior_mean, # nItems x 1
		beta.prior_var,  # nItems x 1
		v.prior_mean,    # nItems x nFactors
		v.prior_var,     # nItems x 1 or nItems x nFactors x nFactors
		as.double(discount.factor),
		nObs, nUsers, nItems,
		nFactors, nBatches, nObsPerBatch, as.integer(length(v.prior_var)),
		as.integer(output.factors), as.integer(debug), as.integer(verbose)
	);
	obs$offset = NULL;

	select = batch_id != -1;
	pred = obs[select,];
	pred$predicted.y = prediction[select];
	pred$batch = batch_id[select];
	out = list(
			prediction = pred,
			beta = beta,
			v = v
		  );
	return(out);
}


###
### RRCS_Gaussian_EStep.C. See RRCS_Gaussian_EStep(...) in RRCS_Gaussian.c
###
### INPUT:  factor  = list(alpha, beta, u, v, s, z); # Initial factor values
###         obs     = data.frame(y, user, item);
###         feature = list(x_dyad, x_user, x_item);
###         param   = list(b, g0, d0, G, D, s_model, z_model,
###                        var_y, var_alpha, var_beta, var_u, var_v, type, B, c);
###
### OUTPUT: mean   = list(beta, v) or list(alpha, u);
###         sumvar = list(beta, v) or list(alpha, u);
###         var    = list(beta, v, cov) or list(alpha, u, cov);
###
RRCS_Gaussian_EStep.C <- function(
		factor, obs, feature, param,
		online.side, # 'item' or 'user'
		outputPostVar=TRUE,
		debug=0, verbose=0
){
	size = syncheck.cRLFM.spec(factor=factor, obs=obs, feature=feature, param=param);

	if(param$type != "q") stop("type must be 'q'");
	if(size$nUserClusters != 1) stop("nUserClusters must be 1");
	if(size$nItemClusters != 1) stop("nItemClusters must be 1");

	if(is.null(factor$alpha)) stop("alpha cannot be null");
	if(is.null(factor$beta))  stop("beta cannot be null");
	if(is.null(factor$u))  stop("u cannot be null");
	if(is.null(factor$v))  stop("v cannot be null");
	if(is.null(param$B))  stop("B cannot be null");
	if(length(obs) == 0) stop("No observation data");

	B = param$B; # nOfflineClusters x nOnlineClusters x nOfflineFactors x nOnlineFactors
	if(length(B) == 1) stop("length(B) == 1");
	if(online.side == "item"){
		offline.side = "user";
		online.main = "beta";
		online.factor = "v";
		offline.main = "alpha";
		offline.factor = "u";
		main.priorMean   = drop(feature$x_item %*% param$d0);
		main.priorVar    = rep(param$var_beta, size$nItems);
		factor.priorMean = feature$x_item %*% param$D;
		factor.priorVar  = rep(param$var_v, size$nItems);
		nOnline  = size$nItems;
		nOffline = size$nUsers;
	}else if(online.side == "user"){
		offline.side = "item";
		online.main = "alpha";
		online.factor = "u";
		offline.main = "beta";
		offline.factor = "v";
		B = transpose_B(param$B);
		main.priorMean   = drop(feature$x_user %*% param$g0);
		main.priorVar    = rep(param$var_alpha, size$nUsers);
		factor.priorMean = feature$x_user %*% param$G;
		factor.priorVar  = rep(param$var_u, size$nUsers);
		nOnline  = size$nUsers;
		nOffline = size$nItems;
	}else stop("side has to be either 'item' or 'user'");
	nOnlineFactors   = ncol(factor[[online.factor]]);
	nOfflineFactors  = ncol(factor[[offline.factor]]);

	transFactor = factor[[offline.factor]] %*% matrix(B[1,1,,], nrow=ncol(factor[[offline.factor]]));
	offset =  drop(feature$x_dyad %*% param$b + factor[[offline.main]][obs[[offline.side]]]);
	y = obs$y;
	var_y = rep(param$var_y, nrow(obs));
	online.index  = obs[[online.side]];
	offline.index = obs[[offline.side]];

	problem.dim = c(size$nObs, nOffline, nOnline, nOnlineFactors, length(factor.priorVar));
	if(!is.integer(problem.dim)) stop("!is.integer(problem.dim)");

	# Allocate space for output
	main.postMean = double(nOnline);
	main.sumvar   = double(1);
	factor.postMean = array(double(1), dim=c(nOnline, nOnlineFactors));
	factor.sumvar   = double(1);
	if(outputPostVar){
		main.postVar   = double(nOnline);
		factor.postVar = array(double(1), dim=c(nOnline, nOnlineFactors, nOnlineFactors));
		postCov        = array(double(1), dim=c(nOnline, nOnlineFactors));
	}else{
		main.postVar = NULL;  factor.postVar = NULL;  postCov=NULL;
	}

	check_type_size(main.priorMean,   "double", problem.dim[3]);
	check_type_size(main.priorVar,    "double", problem.dim[3]);
	check_type_size(factor.priorMean, "double", problem.dim[c(3,4)]);
	check_type_size(factor.priorVar,  "double", problem.dim[c(3)]);
	check_type_size(transFactor,      "double", problem.dim[c(2,4)]);
	check_type_size(offline.index,    "int", problem.dim[1]);
	check_type_size(online.index,     "int", problem.dim[1]);
	check_type_size(y,                "double", problem.dim[1]);
	check_type_size(offset,           "double", problem.dim[1]);
	check_type_size(var_y,            "double", problem.dim[1]);

	#  dim = {nObs:1, nOffline:2, nOnline:3, nOnlineFactors:4, factor_priorVar_length:5}

	ans = .C("RRCS_Gaussian_EStep",
			 #OUTPUT
			main.postMean, # nItems x 1
			main.sumvar,   # 1x1
			main.postVar,  # nItems x 1*/
			factor.postMean,  # nItems x nFactors*/
			factor.sumvar,    # 1x1
			factor.postVar,   # nItems x nFactors x nFactors
	 		postCov,          # nItems x nFactors
			# INPUT
			main.priorMean,   # nItems x 1
			main.priorVar,    # nItems x 1
	        factor.priorMean, # nItems x nFactors
			factor.priorVar,  # nItems x 1
			transFactor,      # nUsers x nFactors
	        offline.index,    # nObs x 1
			online.index,     # nObs x 1
	        y, offset, var_y, problem.dim, length(problem.dim),
			# OTHER
			as.integer(outputPostVar), as.integer(debug), as.integer(verbose),
			DUP=FALSE
	);
	output = list();
	output$mean = list(main.postMean, factor.postMean);
	names(output$mean) = c(online.main, online.factor);
	output$sumvar = list(main.sumvar, factor.sumvar);
	names(output$sumvar) = c(online.main, online.factor);
	if(outputPostVar){
		output$var = list(main.postVar, factor.postVar, postCov);
		names(output$var) = c(online.main, online.factor, "cov");
	}

	return(output);
}

###
### Debug code
###
MCEM_EStep_debug.C <- function(
		factor, obs, feature, param, nSamples, nBurnIn=1,
		s_algo=NULL, z_algo=NULL,
		userFactorVar=0, itemFactorVar=0,
		isOldUser=NULL, isOldItem=NULL,
		debug=0, verbose=0
){
	size = syncheck.cRLFM.spec(factor=factor, obs=obs, feature=feature, param=param);

	if(param$type == "c"){
		stop("type == 'c' is not yet supported");
	}else if(param$type == "q"){
	}else stop("param$type must be either 'c' or 'q'");

	if(size$nUserClusters > 1 && (is.null(s_algo) || is.null(factor$s))) stop("Please specify s_algo and factor$s!");
	if(size$nItemClusters > 1 && (is.null(z_algo) || is.null(factor$z))) stop("Please specify z_algo and factor$z!");
	if(is.null(factor$alpha)) stop("alpha cannot be null");
	if(is.null(factor$beta))  stop("beta cannot be null");
	if(length(obs) == 0) stop("No observation data");

	sumvar  = list(alpha=double(1), beta=double(1), u=double(1), v=double(1), mu=double(1));
	sampvar = list();
	if(userFactorVar != 0){
		sampvar$alpha = double(size$nUsers);
		sampvar$u     = array(double(1),dim=c(size$nUsers, size$nFactorsPerUser, size$nFactorsPerUser));
	}
	if(itemFactorVar != 0){
		sampvar$beta = double(size$nItems);
		sampvar$v    = array(double(1),dim=c(size$nItems, size$nFactorsPerItem, size$nFactorsPerItem));
	}

	xb  = feature$x_dyad %*% param$b;
	g0x_user = NULL; Gx_user = NULL;
	if(is.null(isOldUser)){
		if(!is.null(param$g0)) g0x_user = feature$x_user %*% param$g0;
		if(!is.null(param$G))   Gx_user = feature$x_user %*% param$G;
		if(!is.null(s_algo)) s_logPrior = log(s_algo$predict(param$s_model, feature$x_user));
	}else{
		if(length(isOldUser) != size$nUsers) stop("length(isOldUser) != nUsers");
		x_user.new = feature$x_user[!isOldUser,,drop=FALSE];
		if(!is.null(param$g0)){ g0x_user = factor$alpha;  g0x_user[!isOldUser] = x_user.new %*% param$g0;}
		if(!is.null(param$G)){   Gx_user = factor$u;      Gx_user[!isOldUser,] = x_user.new %*% param$G;}
		if(!is.null(s_algo)){ s_logPrior = log(factor$s); s_logPrior[!isOldUser,] = log(s_algo$predict(param$s_model, x_user.new));}
	}

	d0x_item = NULL;  Dx_item = NULL;
	if(is.null(isOldItem)){
		if(!is.null(param$d0)) d0x_item = feature$x_item %*% param$d0;
		if(!is.null(param$D))   Dx_item = feature$x_item %*% param$D;
		if(!is.null(z_algo)) z_logPrior = log(z_algo$predict(param$z_model, feature$x_item));
	}else{
		if(length(isOldItem) != size$nItems) stop("length(isOldItem) != nItems");
		x_item.new = feature$x_item[!isOldItem,,drop=FALSE];
		if(!is.null(param$d0)){ d0x_item = factor$beta;   d0x_item[!isOldItem] = x_item.new %*% d0;}
		if(!is.null(param$D)){   Dx_item = factor$v;      Dx_item[!isOldItem,] = x_item.new %*% D;}
		if(!is.null(z_algo)){ z_logPrior = log(factor$z); z_logPrior[!isOldItem,] = log(z_algo$predict(param$z_model, x_item.new));}
	}

	# Checks
	if(size$nVar_alpha > 0 && (!is.double(g0x_user) || length(g0x_user) != size$nUsers)) stop("g0x_user!");
	if(size$nVar_beta  > 0 && (!is.double(d0x_item) || length(d0x_item) != size$nItems)) stop("d0x_item!");
	if(size$nVar_u     > 0 && (!is.double(Gx_user)  || any(dim(Gx_user) != c(size$nUsers,size$nFactorsPerUser)))) stop("Gx_user!");
	if(size$nVar_v     > 0 && (!is.double(Dx_item)  || any(dim(Dx_item) != c(size$nItems,size$nFactorsPerItem)))) stop("Dx_item!");
	if(size$nUserClusters > 1 && (!is.double(s_logPrior) || any(dim(s_logPrior) != c(size$nUsers, size$nUserClusters)))) stop("s_logPrior");
	if(size$nItemClusters > 1 && (!is.double(z_logPrior) || any(dim(z_logPrior) != c(size$nItems, size$nItemClusters)))) stop("z_logPrior");

	problem.dim = c(size$nObs, size$nUsers, size$nItems, size$nFactorsPerUser, size$nFactorsPerItem, size$nUserClusters, size$nItemClusters,
			size$nVar_y, size$nVar_alpha, size$nVar_beta, size$nVar_u, size$nVar_v);
	if(!is.integer(problem.dim)) stop("!is.integer(problem.dim)");

	if(!is.double(xb) || length(xb) != size$nObs) stop("xb");

	if(!is.double(factor$alpha) || length(factor$alpha) != problem.dim[2]) stop("alpha");
	if(!is.double(factor$beta) || length(factor$beta) != problem.dim[3]) stop("beta");
	if(!is.double(factor$u) || any(dim(factor$u) != problem.dim[c(2,4)])) stop("u");
	if(!is.double(factor$v) || any(dim(factor$v) != problem.dim[c(3,5)])) stop("v");
	if(!is.double(factor$s) || any(dim(factor$s) != problem.dim[c(2,6)])) stop("s");
	if(!is.double(factor$z) || any(dim(factor$z) != problem.dim[c(3,7)])) stop("z");

	#	//  dim = {1:nObs, 2:nUsers, 3:nItems, 4:nFactorsPerUser, 5:nFactorsPerItem, 6:nUserClusters,
	#   //         7:nItemClusters, 8:nVar_y, 9:nVar_alpha, 10:nVar_beta, 11:nVar_u, 12:nVar_v}

	ans = .C("MCEM_EStep",
			# INPUT (initial factor values) & OUTPUT (Monte Carlo mean of factor values)
			# factor$alpha, factor$beta,
			factor$u, factor$v,
			# factor$s, factor$z,
			# OUTPUT
			#sumvar$alpha,  		sampvar$alpha,
			#sumvar$beta,   		sampvar$beta,
			#sumvar$u, 	sampvar$u,
			sumvar$v,
			#sampvar$v,
			#sumvar$mu,
			# INPUT
			as.integer(nSamples), as.integer(nBurnIn),
			obs$user, obs$item, # obs$y,
			#xb, g0x_user, d0x_item, Gx_user,
			Dx_item, param$B,
			#s_logPrior, z_logPrior,
			param$var_y,
			#param$var_alpha, param$var_beta, param$var_u,
			param$var_v,
			problem.dim, as.integer(length(problem.dim)),
			as.integer(userFactorVar), as.integer(itemFactorVar),
			# OTHER
			as.integer(debug),  as.integer(verbose),
			DUP=FALSE
	);
	#void MCEM_EStep(
	#	// INPUT (initial factor values) & OUTPUT (Monte Carlo mean of factor values)
	#	double* alpha_mean/*nUsers x 1*/,           double* beta_mean/*nItems x 1*/,
	#	double* u_mean/*nUsers x nFactorsPerUser*/, double* v_mean/*nItems x nFactorsPerItem*/,
	#	double* s_mean/*nUsers x nUserClusters*/,   double* z_mean/*nItems x nItemClusters*/,
	#	// OUTPUT
	#	double* alpha_sumvar/*1x1*/,  double* alpha_outputVar/*nUsers x 1*/,
	#	double* beta_sumvar/*1x1*/,   double* beta_outputVar/*nItems x 1*/,
	#	double* u_sumvar/*1x1*/,      double* u_outputVar/*nUsers x nFactorsPerUser x nFactorsPerUser*/,
	#	double* v_sumvar/*1x1*/,      double* v_outputVar/*nItems x nFactorsPerItem x nFactorsPerItem*/,
	#	double* err_sumvar/*1x1*/,
	#	// INPUT
	#	const int* nSamples,                    const int* nBurnIn,
	#	const int* user/*nObs x 1*/,            const int* item/*nObs x 1*/,
	#	const double* y/*nObs x 1*/,            const double* xb/*nObs x 1*/,
	#	const double* g0x_user/*nUsers x 1*/,   const double* d0x_item/*nItems x 1*/,
	#	const double* Gx_user/*nUsers x nFactorsPerUser*/,
	#	const double* Dx_item/*nItems x nFactorsPerItem*/,
	#	const double* B/*nUserClusters x nItemClusters x nFactorsPerUser x nFactorsPerItem*/,
	#	const double* s_logPrior/*nUsers x nUserClusters*/, double* z_logPrior/*nItems x nItemClusters*/,
	#	const double* var_y, const double* var_alpha, const double* var_beta,
	#	const double* var_u, const double* var_v,
	#	const int* dim /*12 x 1*/,      const int* nDim /*must be 12*/,
	#	const int* outputUserFactorVar, const int* outputItemFactorVar,
	#	// OTHER
	#	const int* debug,  const int* verbose
	#)
	output = list(
			mean=factor, sumvar=sumvar, sampvar=sampvar
	);
	return(output);
}
