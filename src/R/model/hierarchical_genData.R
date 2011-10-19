### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen


genData.1D.2Levels <- function(
	nItems, nCategories, nRep, nObsFeatures, nItemFeatures,
	q.mean=0, q.sd=1, w_obs.sd=1, w_b.sd=1, offset.sd=1,
	x_b.with.intercept=TRUE, w_obs.with.intercept=TRUE,
	var_a=5, var_b.mean=2, var_b.sdlog=0.2, var_obs.mean=1, var_obs.sdlog=0.2
){
	nObs = nItems * nCategories * nRep;
	offset = rnorm(nObs, mean=0, sd=offset.sd);
	data = data.frame(
			item     = as.integer(rep(1:nItems, nCategories*nRep)),
			category = as.integer(rep(1:nCategories, each=nItems*nRep)),
			obs      = offset,
			var      = rlnorm(nObs, meanlog=log(var_obs.mean), sdlog=var_obs.sdlog)
	);
	q = rnorm(nCategories, mean=q.mean, sd=q.sd);
	var_b = rlnorm(nCategories, meanlog=log(var_b.mean), sdlog=var_b.sdlog);
	a = rnorm(nItems, mean=0, sd=sqrt(var_a));
	b = matrix(0, nrow=nItems, ncol=nCategories);
	ik = cbind(data$item, data$category);
	
	feature = list();
	w_b = NULL; w_obs = NULL;
	if(nItemFeatures > 0){
		feature$x_b = matrix(rnorm(nItems*nItemFeatures), nrow=nItems);
		if(x_b.with.intercept) feature$x_b[,1] = 1;
		w_b = matrix(rnorm(nCategories*nItemFeatures, mean=0, sd=w_b.sd), nrow=nCategories);
		for(k in 1:nCategories) b[,k] = drop(feature$x_b %*% w_b[k,]);
	}
	for(k in 1:nCategories) b[,k] = b[,k] + rnorm(nItems, mean=q[k]*a, sd=sqrt(var_b[k]));
	
	if(nObsFeatures > 0){
		feature$x_obs = matrix(rnorm(nObs*nObsFeatures), nrow=nObs);
		if(w_obs.with.intercept) feature$x_obs[,1] = 1;
		w_obs = rnorm(nObsFeatures, mean=0, sd=w_obs.sd);
		data$obs = data$obs + drop(feature$x_obs %*% w_obs);
	}
	data$obs = data$obs + rnorm(nObs, mean=b[ik], sd=sqrt(data$var));
	out = list(data=data, feature=feature, offset=offset, a=a, b=b, param=list(q=q, w_obs=w_obs, w_b=w_b, var_a=var_a, var_b=var_b, var_obs_adj=1));
	return(out);
}

split.each.item.1D.2Levels <- function(
	problem, train.frac, randomize=TRUE, by.category=TRUE, 
	id.column="item", category.column="category", verbose=0
){
	data   = problem$data;
	offset = problem$offset;
	x_obs  = problem$feature$x_obs;
	
	if(!(id.column %in% names(data))) stop(id.column," is not a column of data");
	
	inTrain = rep(NA, nrow(data));
	obsIndex = split(data.frame(index=1:nrow(data), category=data[[category.column]], stringsAsFactors=FALSE), data[[id.column]]);
	if(verbose > 0) cat("Number of IDs: ",length(obsIndex),"\n",sep="");
	for(i in 1:length(obsIndex)){
		indices    = obsIndex[[i]]$index;
		categories = obsIndex[[i]]$category;
		if(by.category){
			temp = unique(categories);
			if(randomize) temp = temp[order(runif(length(temp)))];
			s = (1:length(temp)) <= ceiling(length(temp)*train.frac);
			cat.train = temp[s];
			select = categories %in% cat.train;
		}else{
			if(randomize) indices = indices[order(runif(length(indices)))];
			select = (1:length(indices)) <= ceiling(length(indices)*train.frac);
		}
		indices.train = indices[select];
		indices.test  = indices[!select];
		inTrain[indices.train] = 1;
		inTrain[indices.test]  = 0;
		if(verbose > 1) cat("  i=",i,":  #train=",length(indices.train),"  #test=",length(indices.test),"\n",sep="");
	}
	if(any(is.na(inTrain))) stop("error");

	train = problem;
	train$data   = data[inTrain == 1,];
	train$offset = offset[inTrain == 1];
	train$feature$x_obs = x_obs[inTrain == 1,,drop=FALSE];
	
	test = problem;
	test$data   = data[inTrain == 0,];
	test$offset = offset[inTrain == 0];
	test$feature$x_obs = x_obs[inTrain == 0,,drop=FALSE];
	
	return(list(train=train, test=test));
}

split.by.item.1D.2Levels <- function(
	problem, train.frac, randomize=TRUE, 
	id.column="item", verbose=0
){
	data   = problem$data;
	offset = problem$offset;
	x_obs  = problem$feature$x_obs;
	
	if(!(id.column %in% names(data))) stop(id.column," is not a column of data");
	
	items.all = unique(data[[id.column]]);
	if(randomize) items.all = items.all[order(runif(length(items.all)))];
	inTrain = data[[id.column]] %in% items.all[1:ceiling(length(items.all)*train.frac)];
	
	if(any(is.na(inTrain))) stop("error");
	
	train = problem;
	train$data   = data[inTrain,];
	train$offset = offset[inTrain];
	train$feature$x_obs = x_obs[inTrain,,drop=FALSE];
	
	test = problem;
	test$data   = data[!inTrain,];
	test$offset = offset[!inTrain];
	test$feature$x_obs = x_obs[!inTrain,,drop=FALSE];
	
	return(list(train=train, test=test));
}
