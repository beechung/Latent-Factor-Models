### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

check_type_size <- function(x, type, size, isNullOK=FALSE, check.NA=TRUE){
	if(is.null(x)){
		if(all(size == 0) || isNullOK) return(TRUE)
		else stop("The input is null");
	}
	if(type == "double"){
		if(!is.double(x)) stop("The input should be double");
	}else if(type == "integer" || type == "int"){
		if(!is.integer(x)) stop("The input should be integer");
	}else stop("Unknown type: ",type,sep="");
	d = dim(x);
	if(is.null(d)) d = length(x);
	if(length(d) != length(size) || any(d != size)) stop("Dimensionality mismatch: (",paste(d,collapse=" x "),") vs (",paste(size,collapse=" x "),")");
	if(check.NA && any(is.na(x))) stop("Some elements are NA");
}

###
### Re-index users, items and terms
### INPUT: obs    = data.frame(author_id, voter_id, item_id, context_id, y);
###        x_dyad = data.frame(author_id, voter_id, item_id, feature1, feature2, ...)
###        x_user = data.frame(user_id, feature1, feature2, ...)
###        x_item = data.frame(item_id, feature1, feature2, ...)
### 
reindexData <- function(
    obs, x_dyad=NULL, x_user=NULL, x_item=NULL,
    add.intercept=TRUE,          # whether to add an intercept
    UserIDs=sort(unique(c(obs$author_id, obs$voter_id))), 
							     # selected UserIDs; out$x_user[i,] will correspond to userIDs[i]
    ItemIDs=NULL,  # selected ItemIDs; out$x_item[i,] will correspond to ItemIDs[i]
	ContextIDs=NULL,
	other.columns=NULL
){
    if(is.null(obs$author_id)) stop("author_id must be a column in obs");
	if(is.null(obs$voter_id)) stop("voter_id must be a column in obs");
    if(is.null(obs$y)) stop("y must be a column in obs");
    if(!is.null(x_dyad) && is.null(x_dyad$author_id)) stop("author_id must be a column in x_dyad");
	if(!is.null(x_dyad) && is.null(x_dyad$voter_id)) stop("voter_id must be a column in x_dyad");
	if(!is.null(x_dyad) && is.null(x_dyad$item_id)) stop("item_id must be a column in x_dyad");
    if(!is.null(x_user) && is.null(x_user$user_id)) stop("user_id must be a column in x_user");
    if(!is.null(x_item) && is.null(x_item$item_id)) stop("item_id must be a column in x_item");
	if(!is.null(ItemIDs)&& is.null(obs$item_id)) stop("ItemIDs is not null, but obs$item_id is");
	if(!is.null(ContextIDs)&& is.null(obs$context_id)) stop("ContextIDs is not null, but obs$context_id is");
	
    out = list();
    out$obs = data.frame(author=match(obs$author_id, UserIDs),
						 voter=match(obs$voter_id,  UserIDs),
                         y=as.double(obs$y));
	out$userIDs = UserIDs;
	obs.selected = !is.na(out$obs$author) & !is.na(out$obs$voter);
	if(!is.null(obs$item_id)){
		if(is.null(ItemIDs)) ItemIDs = sort(unique(obs$item_id))
		out$obs$item = match(obs$item_id,  ItemIDs);
		obs.selected = obs.selected & !is.na(out$obs$item);
		out$itemIDs = ItemIDs;
	}
	if(!is.null(obs$context_id)){
		if(is.null(ContextIDs)) ContextIDs = sort(unique(obs$context_id))
		out$obs$context = match(obs$context_id,  ContextIDs);
		obs.selected = obs.selected & !is.na(out$obs$context);
		out$contextIDs = ContextIDs;
	}
	if(!is.null(other.columns)){
		for(k in 1:length(other.columns)){
			out$obs[,other.columns[k]] = obs[,other.columns[k]]
		}
	}
	
	nObs.ignored = sum(!obs.selected);
	if(nObs.ignored > 0) warning("NOTE: ",nObs.ignored," are ignored!!");
	
    out$obs = out$obs[obs.selected,]
    out$feature = list();
    if(add.intercept){
        regFormula = formula(~.);
    }else{
        regFormula = formula(~.-1);
    }
    if(is.null(x_dyad)){
        out$feature$x_dyad = matrix(1.0, nrow=nrow(out$obs), ncol=1);
    }else{
        if(any(obs$author_id != x_dyad$author_id)) stop("obs$author_id != x_dyad$author_id");
		if(any(obs$voter_id != x_dyad$voter_id)) stop("obs$voter_id != x_dyad$voter_id");
		if(!is.null(obs$item_id) && any(obs$item_id != x_dyad$item_id)) stop("obs$item_id != x_dyad$item_id");
        out$feature$x_dyad = model.matrix(regFormula, x_dyad[,-match(c("author_id","voter_id","item_id"), names(x_dyad)),drop=FALSE]);
        out$feature$x_dyad = out$feature$x_dyad[obs.selected,]
    }
    if(is.null(x_user)){
        out$feature$x_user = matrix(1.0, nrow=length(UserIDs), ncol=1);
    }else{
        select = match(UserIDs, x_user$user_id);
        if(any(is.na(select))) stop("Some user IDs are not found in x_user");
        out$feature$x_user = model.matrix(regFormula, x_user[select,-match(c("user_id"), names(x_user)),drop=FALSE]);
    }
	
	if(!is.null(obs$item_id)){
	    if(is.null(x_item)){
	        out$feature$x_item = matrix(1.0, nrow=length(ItemIDs), ncol=1);
	    }else{
	        select = match(ItemIDs, x_item$item_id);
	        if(any(is.na(select))) stop("Some item IDs are not found in x_item");
	        out$feature$x_item = model.matrix(regFormula, x_item[select,-match(c("item_id"), names(x_item)),drop=FALSE]);
	    }
	}
    return(out);
}

indexTestData <- function(
    data.train,
    obs, x_dyad=NULL, x_user=NULL, x_item=NULL,
    add.intercept=TRUE, # whether to add an intercept
	other.columns=NULL
){
    UserIDs = c(data.train$userIDs, setdiff(unique(c(obs$author_id, obs$voter_id)), data.train$userIDs));
	if(!is.null(obs$item_id)) ItemIDs = c(data.train$itemIDs, setdiff(unique(obs$item_id), data.train$itemIDs))
	else                      ItemIDs = NULL;
	if(!is.null(obs$context_id)) ContextIDs = c(data.train$contextIDs, setdiff(unique(obs$context_id), data.train$contextIDs))
	else                         ContextIDs = NULL;
	
    data.test = reindexData(
		obs, x_dyad=x_dyad, x_user=x_user, x_item=x_item,
        UserIDs=UserIDs, ItemIDs=ItemIDs, ContextIDs=ContextIDs,
		add.intercept=add.intercept, other.columns=other.columns);
    return(data.test);
}

###
### Check whehter some users or items do not have any rating or some items
### do not have any terms.
###
check.Obs <- function(obs, nUsers, nItems, error=stop){
    temp = sort(unique(c(obs$author, obs$voter)));
    if(length(temp) != nUsers || any(temp != 1:nUsers)) error("Some users do not have observed ratings.");
    temp = sort(unique(obs$item));
    if(length(temp) != nItems || any(temp != 1:nItems)) error("Some items do not have observed ratings.");
}


###
### Checks for user_id, item_id encoding, features
###
### Looks at data.train${userIDs, itemIDs, feature}
###
check.compatible <- function(data.train, data.test, strict=FALSE){
    if(strict || (!is.null(data.train$userIDs) && !is.null(data.test$userIDs))){
        len = min(length(data.train$userIDs), length(data.test$userIDs));
        if(any(data.train$userIDs[1:len] != data.test$userIDs[1:len])) stop("The two datasets have incompatible user ids");
    }
    if(strict || (!is.null(data.train$itemIDs) && !is.null(data.test$itemIDs))){
        len = min(length(data.train$itemIDs), length(data.test$itemIDs));
        if(any(data.train$itemIDs[1:len] != data.test$itemIDs[1:len])) stop("The two datasets have incompatible item ids");
    }
    check.feature.comp(data.train$feature$x_dyad, data.test$feature$x_dyad, "x_dyad", strict);
    check.feature.comp(data.train$feature$x_user, data.test$feature$x_user, "x_user", strict);
    check.feature.comp(data.train$feature$x_item, data.test$feature$x_item, "x_item", strict);
}
check.feature.comp <- function(x1, x2, name, strict=FALSE){
    if(strict || !is.null(x1) || !is.null(x2)){
        if(is.null(x1) || is.null(x2)) stop("One dataset has ",name,"; one doesn't");
        name1 = dimnames(x1)[[2]];
        name2 = dimnames(x2)[[2]];
        if(strict || !is.null(name1) || !is.null(name2)){
            if(length(name1) != length(name2) || any(name1 != name2))
                stop("Feature ",name," not compatible!");
        }
    }
}

###
### Create a deep copy of the input object
###
deepCopy <- function(x){
    if(is.list(x)){
        out = list();
        for(i in 1:length(x)){
            out[[i]] = deepCopy(x[[i]]);
        }
        if(length(out) > 0){
            names(out) = names(x)[1:length(out)];
        }
        return(out);
    }
    # if(is.array(x)) out = array(rep(x), dim=dim(x))
    # else if(is.vector(x)) out = rep(x);
    if(is.integer(x)){
        out = x + as.integer(0);
    }else if(is.numeric(x)){
        out = x + 0;
	}else if(is.logical(x)){
		out = x & TRUE;
	}else if(is.null(x)){
        out = NULL;
    }else stop("Type not supported");
    
    return(out);
}

###
### Predict the response using the factors
###
predict.y.from.factors <- function(obs, factor, feature, param){
    author = obs$author; voter = obs$voter; item = obs$item;
    ans = feature$x_dyad %*% param$b;
    if(!is.null(factor$alpha)) ans = ans + factor$alpha[author];
    if(!is.null(factor$beta))  ans = ans + factor$beta[voter];
    if(!is.null(factor$v))     ans = ans + sum_margin(factor$v[author,,drop=FALSE] * factor$v[voter,,drop=FALSE], 1);
    return(drop(ans));
}

###
### Complete data log likelihood with the constant term removed
###
logLikelihood <- function(
    obs, factor, feature, param, is.logistic
){
    size = syncheck.factorModel.spec(factor=factor, obs=obs, feature=feature, param=param);
    
    nObs     = size$nObs;
    nUsers   = size$nUsers;
    nItems   = size$nItems;
    nFactors = size$nFactors;
    
	pred.y = predict.y.from.factors(obs, factor, feature, param);
    ans = obsLoglik.from.gaussian(pred.y=pred.y, y=obs$response, var_y=param$var_y, is.logistic=is.logistic);

	if(!is.null(param$var_v)){
        ans = ans + loglik.gaussian(pred.x=feature$x_user %*% param$G, x=factor$v, var_x=param$var_v);
    }
    if(!is.null(param$var_alpha)){
		ans = ans + loglik.gaussian(pred.x=feature$x_user %*% param$g0, x=factor$alpha, var_x=param$var_alpha);
    }
    if(!is.null(param$var_beta)){
		ans = ans + loglik.gaussian(pred.x=feature$x_user %*% param$d0, x=factor$beta, var_x=param$var_beta);
    }
   	
    return(ans);
}

###
### Predict for the Gaussian case
###
###     User 1 to length(fit$factor$alpha) are old users
###
### To disable a FACTOR, set var_FACTOR = NULL
###
pred.gauss <- function(
		fit,          # Fitted model
		obs,          # Observation
		feature,      # Feature values
		factorOnly=FALSE, featureOnly=FALSE,
		is.logistic=FALSE,
		metrics=c(),  # an array of the following "kendall", "pearson"
		output.factor=FALSE
){
	
	userIDs = unique(c(obs$author, obs$voter));
	
	nOldUsers = length(fit$factor$alpha);  nAllUsers = max(userIDs);
	
	oldUserIDs <- userIDs[userIDs <= nOldUsers];
	newUserIDs <- userIDs[userIDs > nOldUsers];
	
	if(!factorOnly){
		if(is.null(feature$x_user) || nrow(feature$x_user) < nAllUsers) stop("some new users in test don't have covariates");
	}
	if(is.null(feature$x_dyad) || nrow(feature$x_dyad) != nrow(obs)) stop("is.null(feature$x_dyad) || nrow(feature$x_dyad) != nrow(obs)");
	
	userIDs = c(oldUserIDs, newUserIDs);
	
	factor = list(
			alpha = getReg1DFactor(fit$factor$alpha, feature$x_user, fit$param$g0, oldUserIDs, newUserIDs, factorOnly, featureOnly, disable=is.null(fit$param$var_alpha)),
			beta  = getReg1DFactor(fit$factor$beta,  feature$x_user, fit$param$d0, oldUserIDs, newUserIDs, factorOnly, featureOnly, disable=is.null(fit$param$var_beta)),
			v     = getReg2DFactor(fit$factor$v,     feature$x_user, fit$param$G,  oldUserIDs, newUserIDs, factorOnly, featureOnly, disable=is.null(fit$param$var_v))
	);
	
	obs$author = match(obs$author, userIDs);
	obs$voter  = match(obs$voter,  userIDs);
	
	# The feature part may not have the right size
	# size = syncheck.cRLFM.spec(factor=factor, obs=obs, feature=feature, param=fit$param);
	
	pred.y = predict.y.from.factors(obs, factor, feature, fit$param);
	
	out = predict.response.from.gaussian(pred.y=pred.y, y=obs$y, param=fit$param, is.logistic=is.logistic);
	
	if("kendall" %in% metrics) out$kendall = cor(obs$y, pred.y, method="kendall");
	if("pearson" %in% metrics) out$pearson = cor(obs$y, pred.y, method="pearson");

	if(output.factor) out$factor = factor;
	
	return(out);
}

###
### Input (factor, feature, param) is one of:
###       (alpha, x_user, g0), (beta, x_item, d0)
###
### output[1:length(oldIDs)] the factor for oldIDs
### output[length(oldIDs)+(1:length(newIDs))] the factor for newIDs
###
getReg1DFactor <- function(factor, feature, param, oldIDs, newIDs, factorOnly, featureOnly, disable, defaultValue=0){
	if(factorOnly && featureOnly) stop("factorOnly && featureOnly");
	if(disable) return(NULL);
	
	out = rep(NA, length(oldIDs)+length(newIDs));
	if(featureOnly){
		out[] = feature[c(oldIDs, newIDs),,drop=FALSE] %*% param;
	} else {
		if(length(oldIDs) > 0) out[1:length(oldIDs)] = factor[oldIDs];
		if(length(newIDs) > 0){
			out[length(oldIDs)+(1:length(newIDs))] <- if(factorOnly) defaultValue else feature[newIDs,,drop=FALSE] %*% param;
		}
	}
	return(out);
}

###
### Input (factor, feature, param) is one of:
###       (u, x_user, G), (v, x_item, D)
###
### output[1:length(oldIDs),] the factor for oldIDs
### output[length(oldIDs)+(1:length(newIDs)),] the factor for newIDs
###
getReg2DFactor <- function(factor, feature, param, oldIDs, newIDs, factorOnly, featureOnly, disable, defaultValue=0){
	if(factorOnly && featureOnly) stop("factorOnly && featureOnly");
	nFactors = 0;
	if(!is.null(factor)) nFactors = ncol(factor);
	if(!is.null(param)){
		if(nFactors != 0 && nFactors != ncol(param)) stop("error");
		nFactors = ncol(param);
	}
	
	if(disable || nFactors == 0) return(matrix(0,nrow=length(oldIDs)+length(newIDs), ncol=0));
	
	out = matrix(NA, nrow=length(oldIDs)+length(newIDs), ncol=nFactors);
	if(featureOnly){
		out[,] = feature[c(oldIDs, newIDs),,drop=FALSE] %*% param;
	} else {
		if(length(oldIDs) > 0) out[1:length(oldIDs),] = factor[oldIDs,];
		if(length(newIDs) > 0){
			out[length(oldIDs)+(1:length(newIDs)),] <- if(factorOnly) defaultValue else feature[newIDs,,drop=FALSE] %*% param;
		}
	}
	return(out);
}

###
### Compute the loglikelihood based on Logistic regression
###
test.loglike.logistic <- function(fit, obs, feature){
	logOdds = pred.gauss(fit=fit, obs=obs, feature=feature)$pred.y;
	label = obs$y;
	if(!all(label %in% c(0,1))) stop("obs$y should be either 0 or 1");
	loglik = -sum(log1p( exp( -(2*label-1)*logOdds ) ));
	return(loglik);
}
test.RMSE <- function(fit, obs, feature){
	return(pred.gauss(fit=fit, obs=obs, feature=feature)$rmse);
}

rep_matrix <- function(matrix, num){
    ans = array(NA, dim=c(num, nrow(matrix), ncol(matrix)));
    for(i in 1:num){
        ans[i,,] = matrix;
    }
    return(ans);
}

###
### Syntactic check of model specification
###
###     factor  = list(alpha, beta, v);
###     obs     = data.frame(author, voter, item, y);
###     feature = list(x_dyad, x_user, x_item);
###     param   = list(b, g0, d0, G, var_y, var_alpha, var_beta, var_v);
###
syncheck.factorModel.spec <- function(
    factor, obs, feature, param, warning=0, print=FALSE,
     factor.name = c("alpha", "beta", "v"),
        obs.name = c("author", "voter", "item", "y"),
    feature.name = c("x_dyad", "x_user", "x_item"),
      param.name = c("b", "g0", "d0", "G", "var_y", "var_alpha", "var_beta", "var_v")
){
    if(warning > 1){
        warning.any.not.in(names(factor),  c(factor.name, "fErr",   "rm.factors.without.obs.in.loglik", "pred.y.square"), "You specified the following unnecessary components in factor: ", print=print);
        warning.any.not.in(names(obs),     c(obs.name, "response"), "You specified the following unnecessary components in obs: ", print=print);
        warning.any.not.in(names(feature), feature.name,"You specified the following unnecessary components in feature: ", print=print);
        warning.any.not.in(names(param),   c(param.name, "xi"), "You specified the following unnecessary components in param: ", print=print);
    }
    if(warning > 2){
        warning.any.not.in(factor.name,  names(factor),  "You did not specify the following components in factor: ", print=print);
        warning.any.not.in(obs.name,     names(obs),     "You did not specify the following components in obs: ", print=print);
        warning.any.not.in(feature.name, names(feature), "You did not specify the following components in feature: ", print=print);
        warning.any.not.in(param.name,   names(param),   "You did not specify the following components in param: ", print=print);
    }
    out = list(
        nObs=get.size(nrow(obs)),
        nUsers=get.size(length(factor[["alpha"]]), length(factor[["beta"]]), nrow(factor[["v"]]), nrow(feature[["x_user"]])), 
        nItems=get.size(max(obs$item), nrow(feature[["x_item"]])),
        nFactors=get.size(ncol(factor[["v"]]),ncol(param[["G"]])),
        nDyadicFeatures=get.size(ncol(feature[["x_dyad"]]),length(param[["b"]])),
        nUserFeatures=get.size(ncol(feature[["x_user"]]), length(param[["g0"]]), length(param[["d0"]]), nrow(param[["G"]])),
        nItemFeatures=get.size(ncol(feature[["x_item"]])),
        nVar_y=as.integer(length(param[["var_y"]])), 
        nVar_alpha=as.integer(length(param[["var_alpha"]])), 
        nVar_beta=as.integer(length(param[["var_beta"]])), 
        nVar_v=as.integer(length(param[["var_v"]]))
    );
    if(out$nVar_v > 1) out$nVar_v = out$nUsers;
    
    if(!is.null(factor[["alpha"]])) check.individual("factor$alpha",is.double(factor[["alpha"]]),length(factor[["alpha"]]),out$nUsers, 
                                    stopIfNull=c("obs$y"=is.null(obs[["y"]]),"feature$x_user"=is.null(feature[["x_user"]]),"param$g0"=is.null(param[["g0"]])));
    if(!is.null(factor[["beta"]]))  check.individual("factor$beta", is.double(factor[["beta"]]), length(factor[["beta"]]), out$nUsers, 
                                    stopIfNull=c("obs$y"=is.null(obs[["y"]]),"feature$x_user"=is.null(feature[["x_user"]]),"param$d0"=is.null(param[["d0"]])));
    if(!is.null(factor[["v"]])) check.individual("factor$v",is.double(factor[["v"]]),dim(factor[["v"]]),c(out$nUsers,out$nFactors), 
                                stopIfNull=c("obs$y"=is.null(obs[["y"]]),"feature$x_user"=is.null(feature[["x_user"]]),"param$G"=is.null(param[["G"]])));
    if(length(obs) > 0){
        warning.any.not.in(obs.name, names(obs), "You did not specify the following components in obs: ", stop=TRUE);
        check.individual("obs$y",is.double(obs[["y"]]),length(obs[["y"]]),out$nObs);
        check.individual("obs$author",is.integer(obs[["author"]]),length(obs[["author"]]),out$nObs);
		check.individual("obs$voter",is.integer(obs[["voter"]]),length(obs[["voter"]]),out$nObs);
		check.individual("obs$item",is.integer(obs[["item"]]),length(obs[["item"]]),out$nObs);
    }
    if(!is.null(feature[["x_dyad"]])) check.individual("feature$x_dyad",is.double(feature[["x_dyad"]]),dim(feature[["x_dyad"]]),c(out$nObs,out$nDyadicFeatures));
    if(!is.null(feature[["x_user"]])) check.individual("feature$x_user",is.double(feature[["x_user"]]),dim(feature[["x_user"]]),c(out$nUsers,out$nUserFeatures));
    if(!is.null(feature[["x_item"]])) check.individual("feature$x_item",is.double(feature[["x_item"]]),dim(feature[["x_item"]]),c(out$nItems,out$nItemFeatures));

    if(!is.null(param[["b"]]))  check.individual("param$b", is.double(param[["b"]]), length(param[["b"]]), out$nDyadicFeatures);
    if(!is.null(param[["g0"]])) check.individual("param$g0",is.double(param[["g0"]]),length(param[["g0"]]),out$nUserFeatures);
    if(!is.null(param[["d0"]])) check.individual("param$d0",is.double(param[["d0"]]),length(param[["d0"]]),out$nUserFeatures);
    if(!is.null(param[["G"]]))  check.individual("param$G", is.double(param[["G"]]), dim(param[["G"]]), c(out$nUserFeatures,out$nFactors));

    if(!is.null(param[["var_y"]]))     check.individual("param$var_y",    is.double(param[["var_y"]]),    length(param[["var_y"]]),    out$nObs,  okDim=1);
    if(!is.null(param[["var_alpha"]])) check.individual("param$var_alpha",is.double(param[["var_alpha"]]),length(param[["var_alpha"]]),out$nUsers,okDim=1);
    if(!is.null(param[["var_beta"]]))  check.individual("param$var_beta", is.double(param[["var_beta"]]), length(param[["var_beta"]]), out$nUsers,okDim=1);
    if(!is.null(param[["var_v"]]))     check.individual("param$var_v",    is.double(param[["var_v"]]),    first.not.null(dim(param[["var_v"]]),length(param[["var_v"]])), c(out$nUsers,out$nFactors,out$nFactors), okDim=1);
    
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
    return(as.integer(s));
}
check.individual <- function(name, typeOK, dim1, dim2, okDim=NULL, stopIfNull=NULL){
    if(!typeOK) stop(name,"'s data type is not correct!");
    if(!is.null(stopIfNull)){
        temp = names(stopIfNull)[stopIfNull];
        if(length(temp) > 0) stop("When ",name," is specified, the following cannot be NULL: ",paste(temp,collapse=", "));
    }
    if(!is.null(okDim) && length(dim1) == length(okDim) && all(dim1 == okDim)) return();
    if(length(dim1) != length(dim2)) stop(name,"'s dim (or length) is not correct: (",paste(dim1,collapse=","),") vs (",paste(dim2,collapse=","),")");
    if(any(dim1 != dim2)) stop(name,"'s dim (or length) is not correct: (",paste(dim1,collapse=","),") vs (",paste(dim2,collapse=","),")");
}
first.not.null <- function(...){
    temp = list(...);
    for(i in 1:length(temp)){
        if(!is.null(temp[[i]])) return(temp[[i]]);
    }
    return(NULL);
}

###
### Check whether two objects x1 and x2 are different, where x1 and x2 can be lists of lists.
###
is.diff <- function(x1, x2, precision=1e-10, prefix=""){
    if(length(x1) != length(x2)){
        cat(prefix,": Different length! (",length(x1)," vs ",length(x2),")\n",sep="");
        return(TRUE);
    }
    if(length(x1) == 0) return(FALSE);
    if(is.list(x1)){
        if(!is.list(x2)){
            cat(prefix,": Different types! (list vs non-list)\n",sep="");
            return(TRUE);
        }
        name1 = sort(names(x1));
        name2 = sort(names(x2));
        if(is.null(name1) || is.null(name2)){
            if(!is.null(name2) || !is.null(name1)){
                cat(prefix,": One has no names; the other has names!\n",sep="");
                return(TRUE);
            }
            ans = FALSE;
            for(i in 1:length(x1)){
                ret = is.diff(x1[[i]], x2[[i]], precision, prefix=paste(prefix,"[[",i,"]]",sep=""));
                ans = ans || ret;
            }
            return(ans);
        }else{
            if(any(name1 != name2)){
                cat(prefix,": Different names!\n",sep="");
                return(TRUE);
            }
            ans = FALSE;
            for(i in 1:length(name1)){
                name = name1[i];
                ret = is.diff(x1[[name]], x2[[name]], precision, prefix=paste(prefix,"$",name,sep=""));
                ans = ans || ret;
            }
            return(ans);
        }
    }else{
        indexDiff = (1:length(x1))[x1 != x2];
        # print(indexDiff);
        if(length(indexDiff) == 0) return(FALSE);
        
        value = cbind(x1[indexDiff], x2[indexDiff]);
        denom = apply(abs(value), 1, max);
        diff  = abs(value[,1]-value[,2]) / denom;
        
        # print(diff);
        
        temp = (1:length(indexDiff))[diff > precision];
        indexDiff = indexDiff[temp];
        diff      = diff[temp];
        if(length(indexDiff) == 0) return(FALSE);
        
        for(i in 1:length(indexDiff)){
            index = indexDiff[i];
            cat(prefix,"[",index,"]: ",x1[index]," vs ",x2[index]," (diff=",diff[i],")\n",sep="");
        }
        return(TRUE);
    }
}

# Get multivariate normal sample
#   mean[k,] is the mean vector for the kth sample point
#   var[k,,] is the variance-covariance matrix for the kth sample point
# output[k,] is the kth sample point
getMVNSample <- function(mean, var=NULL, var.inv=NULL, FUN=my_rmvnorm){
    
    nPoints = nrow(mean);
    nDim    = ncol(mean);
    if(!is.null(var)) temp = dim(var)
    else temp = dim(var.inv) ;
    if(temp[1] != nPoints || temp[2] != nDim || temp[3] != nDim) stop("size mismatch");
    
    if(nDim == 1) return(matrix(rnorm(nPoints, mean, sqrt(var)), nrow=nPoints, ncol=1));
    
    output = matrix(NA, nrow=nPoints, ncol=nDim);
    for(k in 1:nPoints){
        if(!is.null(var.inv)){
            output[k,] = FUN(1, mu=mean[k,], Sigma.inv=var.inv[k,,]);
        }else{
            output[k,] = FUN(1, mu=mean[k,], Sigma=var[k,,]);
        }
    }
    return(output);
}

subsample.ROC <- function(perf, nPoints=1000){
    n = length(perf@x.values[[1]]);
    size = floor(n / nPoints);
    select = ((1:n) %% size) == 1;
    select[n] = TRUE;
    perf@x.values[[1]] = perf@x.values[[1]][select];
    perf@y.values[[1]] = perf@y.values[[1]][select];
    perf@alpha.values[[1]] = perf@alpha.values[[1]][select];
    return(perf);
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
	ans = .C("sum_margin", out, A, as.integer(nrow(A)), as.integer(ncol(A)), as.integer(side), DUP=FALSE);
	return(out);
}

### TODO: IMPORTANT!! This formula needs to be validated!!
###
### Create Gaussian response for the Logistic case
###		label should be either 0 or 1
###
###	Ref: TOMMI S. JAAKKOLA and MICHAEL I. JORDAN, Bayesian parameter estimation via
###		 variational methods. Statistics and Computing (2000) 10, 25–37.
###		 Eq (10) and (11)
###
gaussResponse.forLogistic <- function(label, pred.logOdd){
	if(!all(label %in% c(0,1))) stop("label should be either 0 or 1");
	var_y = 2 * pred.logOdd / tanh(pred.logOdd/2);
	y = (label - 1/2) * var_y;
	return(list(y=y, var_y=var_y));
}
