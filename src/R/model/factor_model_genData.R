### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen


# Generate some normal data for testing
#
#   All feature values (x_dyad, x_user, x_item) follow N(0,1)
#   
#   output$obs is a data.frame with columns (author, voter, item, rating)
#                                               i      j      c      y
#   output$feature$x_dyad[m,] is the feature vector for the mth observation
#   output$feature$x_user[i,] is the feature vector for user i
#   output$feature$x_item[j,] is the feature vector for item j
#   
#       x_user[output$obs[m,"user"],] gives the user feature vector for the mth observation
#       x_item[output$obs[m,"item"],] gives the item feature vector for the mth observation
#       x_dyad[k,] gives the dyadic feature vector fot the kth observation
#
#   output$factor$alpha[i]: The main author effect for user i
#   output$factor$beta[j]:  The main voter  effect for user j
#   output$factor$v[i,k]:   The factor vector for user i
#
#   output$param = a copy of the input parameters
#
genNormalData <- function(
    nUsers, nItems, nObs,
    b, g0, d0, G,
	var_y, var_alpha, var_beta, var_v,
	binary.response=FALSE,
	binary.features=FALSE
){
    if(!is.vector(b))  stop("b should be a vector");
    if(!is.vector(g0)) stop("g0 should be a vector");
    if(!is.vector(d0)) stop("d0 should be a vector");

	nDyadicFeatures = length(b);
	nUserFeatures   = length(g0);
	nItemFeatures   = 1;  # currently not used
	
	if(!is.null(G) && length(dim(G)) != 2)  stop("G should be a matrix");
	if(!is.null(G) && dim(G)[1] != nUserFeatures) stop("dim(G)[1] != nUserFeatures");
	
	nFactors = dim(G)[2]; if(is.null(nFactors)) nFactors = 0;
		
    if(nObs < nUsers || nObs < nItems) stop("nObs < nUsers || nObs < nItems");
    
    x_dyad = matrix(rnorm(nObs*nDyadicFeatures), nrow=nObs,   ncol=nDyadicFeatures,dimnames=list(NULL, sprintf("x_dyad_%03d", 1:nDyadicFeatures)));
    x_user = matrix(rnorm(nUsers*nUserFeatures), nrow=nUsers, ncol=nUserFeatures,  dimnames=list(NULL, sprintf("x_user_%03d", 1:nUserFeatures)));
    x_item = matrix(rnorm(nItems*nItemFeatures), nrow=nItems, ncol=nItemFeatures,  dimnames=list(NULL, sprintf("x_item_%03d", 1:nItemFeatures)));
    
	if(binary.features){
		x_dyad[,] = x_dyad > 0;
		x_user[,] = x_user > 0;
		x_item[,] = x_item > 0;
	}
	x_dyad[,1] = 1;
	x_user[,1] = 1;
	x_item[,1] = 1;
			
	alpha = rnorm(nUsers, mean=x_user %*% g0, sd=sqrt(var_alpha));
	beta  = rnorm(nUsers, mean=x_user %*% d0, sd=sqrt(var_beta));
	
	if(!is.null(G)){
			v = x_user %*% G + rnorm(nUsers*nFactors, mean=0, sd=sqrt(var_v));
	}else{
			v = matrix(0, nrow=nUsers, ncol=nFactors);
	}
	
    author = c(1:nUsers, floor(runif(nObs-nUsers, min=1, max=nUsers+1)));
	voter  = floor(runif(nObs, min=1, max=nUsers));
	voter[voter>=author] = voter[voter>=author]+1;
    item   = c(floor(runif(nObs-nItems, min=1, max=nItems+1)), 1:nItems);

	output=list();
	output$obs = data.frame(author=as.integer(author), voter=as.integer(voter), item=as.integer(item));
	output$feature = list(x_user=x_user, x_item=x_item, x_dyad=x_dyad);
	output$factor  = list(alpha=alpha, beta=beta, v=v);
	output$param   = list(b=b, g0=g0, d0=d0, G=G, var_y=var_y, var_alpha=var_alpha, var_beta=var_beta, var_v=var_v);
	
	pred.y = drop(predict.y.from.factors(output$obs, output$factor, output$feature, output$param));
	if(binary.response){
		output$obs$y = rbinom(n=nObs,size=1,prob=1/(1+exp(-pred.y)));
	}else{
		output$obs$y =  pred.y + rnorm(nObs, mean=0, sd=sqrt(var_y));
	}
    
    if(is.null(G)) output$factor$v = NULL;

    return(output);
}

generate.GaussianData <- function(
	nObs, nUsers, nItems, nFactors,
	nDyadicFeatures, nUserFeatures,
	b.sd=1, g0.sd=1, d0.sd=1, G.sd=1,
	var_y, var_alpha, var_beta, var_v,
	binary.features=FALSE, binary.response=FALSE
){
	b = rnorm(nDyadicFeatures,0,b.sd);
	g0 = rnorm(nUserFeatures, 0, g0.sd);
	d0 = rnorm(nUserFeatures, 0, d0.sd);
	if(nFactors > 0) G = matrix(rnorm(nUserFeatures*nFactors, 0, G.sd), nrow=nUserFeatures)
	else             G = NULL;
	
	ans = genNormalData(
			nUsers=nUsers, nItems=nItems, nObs=nObs,
			b=b, g0=g0, d0=d0, G=G,
			var_y=var_y, var_alpha=var_alpha, var_beta=var_beta, var_v=var_v,
			binary.features=binary.features, binary.response=binary.response
	);
	return(ans);
}

###
### Split a dataset into train and test
###
split.forTesting <- function(data, nAuthors.testOnly, nItems.testOnly, nObs.test){
    size = syncheck.factorModel.spec(data$factor, data$obs, data$feature, data$param);
    user.test = ((1:size$nUsers)[order(runif(size$nUsers))])[1:nAuthors.testOnly];
    item.test = ((1:size$nItems)[order(runif(size$nItems))])[1:nItems.testOnly];
    obs.forTest = (data$obs$author %in% user.test) | (data$obs$item %in% item.test);
    obsID.train = (1:size$nObs)[!obs.forTest];
    obsID.test  = (1:size$nObs)[ obs.forTest];
    need = nObs.test - length(obsID.test);
    if(need < 0){
        warning("nObs.test = ",nObs.test,", but we have already allocated ", length(obsID.test), " for testing because of test-only users and items!");
    }else if(need > 0){
        obs.forTest = ((1:length(obsID.train))[order(runif(length(obsID.train)))])[1:need];
        obsID.test  = c(obsID.test, obsID.train[obs.forTest]);
        obsID.train = obsID.train[-obs.forTest];
    }
    
    userID.train = sort(unique(c(data$obs$author[obsID.train], data$obs$voter[obsID.train])));
    userID.test  = sort(unique(c(data$obs$author[obsID.test ], data$obs$voter[obsID.test])));
    userID.test  = c(userID.train, userID.test[!(userID.test %in% userID.train)]);
    
    itemID.train = sort(unique(data$obs$item[obsID.train]));
    itemID.test  = sort(unique(data$obs$item[obsID.test]));
    itemID.test  = c(itemID.train, itemID.test[!(itemID.test %in% itemID.train)]);
    
    out = list();
    out$train = selectDataSubset(data, obsID.train, userID.train, itemID.train);
    out$test  = selectDataSubset(data, obsID.test,  userID.test,  itemID.test);
    
    size.train = syncheck.factorModel.spec(out$train$factor, out$train$obs, out$train$feature, out$train$param);
    size.test  = syncheck.factorModel.spec(out$test$factor, out$test$obs, out$test$feature, out$test$param);
    
    if(any(sort(unique(c(out$train$obs$author,out$train$obs$voter))) != 1:size.train$nUsers)) stop("error");
    if(any(sort(unique(out$train$obs$item)) != 1:size.train$nItems)) stop("error");
    if(max(out$test$obs$author, out$test$obs$voter) != size.test$nUsers) stop("error");
    if(max(out$test$obs$item) != size.test$nItems) stop("error");
    if(size.test$nUsers != size$nUsers) stop("error");
    if(size.test$nItems != size$nItems) stop("error");
    if(size.train$nItems + nItems.testOnly > size$nItems) stop("error");
    
    return(out);
}

###
### Select a subset of data by obsID (row number), userID and itemID
###
###                 [1]  [2]  [3]  [4]  [5]  <- New ID in the output data
### E.g. userID = c( 1,   2,   4,   8,   9)  <- Original ID in input data
###                                             (ID in data$user, row number in data$factor$alpha, ...)
###
### obsID: the selected observation IDs (row number in the input data$obs)
###
selectDataSubset <- function(data, obsID, userID, itemID){
    # user.newID[oldID] = newID
    # item.newID[oldID] = newID
    user.newID = c(); item.newID = c();
    user.newID[userID] = as.integer(1:length(userID));
    item.newID[itemID] = as.integer(1:length(itemID));
    newObs = data.frame(
        author = user.newID[data$obs$author[obsID]],
		voter  = user.newID[data$obs$voter[obsID]],
		item   = item.newID[data$obs$item[obsID]],
        y      = data$obs$y[obsID]
    );
    newFeature = list(
        x_user=data$feature$x_user[userID,,drop=FALSE], x_item=data$feature$x_item[itemID,,drop=FALSE], x_dyad=data$feature$x_dyad[obsID,,drop=FALSE]
    );
    newFactor = list(
        alpha=data$factor$alpha[userID], beta=data$factor$beta[userID], v=data$factor$v[userID,,drop=FALSE]
	);
    newData = list(
        obs = newObs, feature = newFeature,
        factor = newFactor, param = data$param
    );
    return(newData);
}
