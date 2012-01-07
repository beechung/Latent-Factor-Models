### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### Find (theta1, theta2) that minimize the following loss
###
###    loss = sum_i E_{y,z}[ (y[i] - theta1' z[i,] - theta2' x[i,])^2 ]
###         = E[sum_i y[i]^2] - 2 theta' a + theta' A theta,
###
### where theta = c( theta1,  theta2 )
###           a = c( E[ sum_i z[i,]*y[i] ],  sum_i x[i,]*E[y[i]] );
###           A = rbind(cbind( E[ sum_i z[i,] z[i,]' ],  sum_i E[z[i,]] x[i,]' ),
###                     cbind( sum_i x[i,] E[z[i,]]',    sum_i x[i,] x[i,]'    ))
###
### Log-likelihood = -(1/2) * ( N log(var_y) + loss / var_y )
###
### INPUT:
###    y[i]     = E[y[i]]
###    sum.y.sq = E[sum_i y[i]^2]
###    z[i,]    = E[z[i,]]
###    sum.zz   = E[ sum_i z[i,] z[i,]' ]
###    x[i,]    = x[i,]
###    sum.zy   = E[ sum_i z[i,]*y[i] ]
###    num      = number of cases involved in sum_i
###
###    subset: A vector of indices for i
###
###	   algo: the fitting algorithm to be used to fit theta2
###          In this case, theta2' x[i,] = theta2(x[i,]); and
###			 theta2 will be fitted first (then treated as a constant when fitting theta1).
###
### OUTPUT: coeff.z=theta1, coeff.x=theta2, loss
###
fit.E.gaussin.loglik <- function(
	y, sum.y.sq, x, z=NULL, sum.zz=NULL, sum.zy=NULL, num, subset=NULL, 
	lambda=0.1, algo=NULL, algo.control=NULL
){
	nCases = length(y);  nFeatures = ncol(x);
	if(!is.null(z)){
		if(is.matrix(z))      nFactors = ncol(z)
		else if(is.vector(z)) nFactors = 1
		else stop("z is not matrix or vector");
		if(nFactors >  1 && nrow(z)   != nCases) stop("nrow(z) != nCases");  
		if(nFactors == 1 && length(z) != nCases) stop("length(z) != nCases");  
		if(length(sum.zy) != nFactors) stop("length(sum.zy) != nFactors");
		if(nFactors >  1 && any(dim(sum.zz) != c(nFactors, nFactors))) stop("any(dim(sum.z.sq) != c(nFactors, nFactors))");
		if(nFactors == 1 && length(sum.zz) != 1) stop("length(sum.zz) != 1");
	}else{
		nFactors = 0;
		if(!is.null(sum.zz)) stop("!is.null(sum.zz)");
		if(!is.null(sum.zy)) stop("!is.null(sum.zy)");
	}
	if(nrow(x) != nCases) stop("nrow(x) != nCases");
	if(length(sum.y.sq) != 1) stop("length(sum.y.sq) != 1")
	if(!is.null(subset)){
		if(is.logical(subset)){
			if(length(subset) != nCases) stop("is.logical(subset) && length(subset) != nCases");
			if(all(subset)) subset = NULL; # no need to select
		}else if(length(subset) == nCases){
			if(all(subset == 1:nCases)) subset = NULL; # no need to select
		}
	}
	if(!is.null(subset)){
		y = y[subset];  x = x[subset,,drop=FALSE];
		if(!is.null(z)) z = z[subset,,drop=FALSE];
	}
	if(length(y) != num) stop("length(y) != num");
	if(is.null(algo)){
		a = c(drop(sum.zy), drop(t(x)%*%y));
		A = matrix(NA, nrow=nFactors+nFeatures, ncol=nFactors+nFeatures);
		if(!is.null(z)){
			A[1:nFactors,1:nFactors] = sum.zz;
			zx = t(z) %*% x;
			A[1:nFactors,nFactors+(1:nFeatures)] = zx;
			A[nFactors+(1:nFeatures),1:nFactors] = t(zx);
		}
		A[nFactors+(1:nFeatures),nFactors+(1:nFeatures)] = t(x) %*% x;
		A.inv = solve(A + lambda*diag(nFactors+nFeatures));
		theta = drop(A.inv %*% a);
		out = list(
				coeff.z = if(nFactors>0) theta[1:nFactors] else NULL,
				coeff.x = theta[nFactors+(1:nFeatures)],
				loss = sum.y.sq - 2 * t(theta) %*% a + t(theta) %*% A %*% theta
		);
	}else{
		theta2 = reg.train(x=x, y=y, algo=algo, control=algo.control);
		theta2.x = reg.predict(model=theta2, x=x, algo=algo);
		if(!is.null(z)){
			E.diff = sum.zy - t(z)%*%theta2.x;
			theta1 = solve(sum.zz + lambda*diag(nFactors)) %*% E.diff;
			loss = sum.y.sq - 2*sum(y*theta2.x) + sum(theta2.x^2) - 2*t(theta1)%*%E.diff + t(theta1)%*%sum.zz%*%theta1;
		}else{
			theta1 = NULL;
			loss = sum.y.sq - 2*sum(y*theta2.x) + sum(theta2.x^2);
		}
		out = list(coeff.z=theta1, coeff.x=theta2, loss=loss);
	}
	return(out);
}
loss.E.gaussin.loglik <- function(
	y, sum.y.sq, x, coeff.x, z=NULL, coeff.z=NULL,
	sum.zz=NULL, sum.zy=NULL, num, subset=NULL, algo=NULL
){
	nCases = length(y);  nFeatures = ncol(x);
	if(!is.null(z)){
		if(is.matrix(z))      nFactors = ncol(z)
		else if(is.vector(z)) nFactors = 1
		else stop("z is not matrix or vector");
		if(nFactors >  1 && nrow(z)   != nCases) stop("nrow(z) != nCases");  
		if(nFactors == 1 && length(z) != nCases) stop("length(z) != nCases");  
		if(nFactors >  1 && any(dim(sum.zz) != c(nFactors, nFactors))) stop("any(dim(sum.z.sq) != c(nFactors, nFactors))");
		if(nFactors == 1 && length(sum.zz) != 1) stop("length(sum.zz) != 1");
		if(length(sum.zy) != nFactors) stop("length(sum.zy) != nFactors");
	}else{
		nFactors = 0;
		if(!is.null(sum.zz)) stop("!is.null(sum.zz)");
		if(!is.null(sum.zy)) stop("!is.null(sum.zy)");
		if(!is.null(coeff.z)) stop("!is.null(coeff.z)");
	}
	if(nrow(x) != nCases) stop("nrow(x) != nCases");
	if(length(sum.y.sq) != 1) stop("length(sum.y.sq) != 1")
	if(!is.null(subset)){
		if(is.logical(subset)){
			if(length(subset) != nCases) stop("is.logical(subset) && length(subset) != nCases");
			if(all(subset)) subset = NULL; # no need to select
		}else if(length(subset) == nCases){
			if(all(subset == 1:nCases)) subset = NULL; # no need to select
		}
	}
	if(!is.null(subset)){
		y = y[subset];  x = x[subset,,drop=FALSE];
		if(!is.null(z)) z = z[subset,,drop=FALSE];
	}
	if(length(y) != num) stop("length(y) != num");
	pred.x = reg.predict(model=coeff.x, x=x, algo=algo);
	ans = sum.y.sq - 2*sum(pred.x * y);
	if(!is.null(z)){
		ans = ans - 2*t(coeff.z)%*%sum.zy;		
		ans = ans + t(coeff.z) %*% sum.zz %*% coeff.z;
		ans = ans + 2 * t(coeff.z) %*% t(z) %*% pred.x;
	}
	ans = ans + sum(pred.x^2);
	return(ans);
}
E.gaussin.loglik <- function(
	y, var_y, sum.y.sq, x, coeff.x, z=NULL, coeff.z=NULL, 
	sum.zz=NULL, sum.zy=NULL, num, subset=NULL, algo=NULL
){
	# sanity check is performed inside loss.E.gaussin.loglik
	loss = loss.E.gaussin.loglik(
			coeff.z=coeff.z, coeff.x=coeff.x, y=y, z=z, x=x, 
			sum.y.sq=sum.y.sq, sum.zz=sum.zz, sum.zy=sum.zy, num=num, subset=subset, algo=algo
	);
	loglik = -(1/2)*(num*log(var_y) + loss / var_y);
	return(loglik);
}

###
### Find the beta and var_y that maximize
###		-(1/2) * sum_i E[ log(var_y) + (w_i/var_y)*(y_i - beta'x_i)^2 ]
### i.e., minimize
###		    sum(w_i > 0) * log(var_y) + (1/var_y) * loss,
###     where loss = sum_i w_i*{ (E[y_i] - beta'E[x_i])^2 + beta'Var[x_i]beta - 2beta'Cov[y_i,x_i] + Var[y_i] }
###
###	Let XX = sum_i { w_i * (Var[x_i] + E[x_i]E[x_i]') }
###     XY = sum_i { w_i * (E[x_i]E[y_i] + Cov[y_i,x_i]) }
### The solution is
###		beta = XX^-1 XY
###    var_y = loss / sum(w_i > 0)
###          = (beta' XX beta - 2 XY beta + sum_i w_i*{E[y_i]^2 + Var[y_i]}) / sum(w_i > 0)
###
### INPUT:
###	   response.mean[i]  = E[y_i]
###    response.var[i]   = Var[y_i]
###     feature.mean[i,] = E[x_i]
###     feature.var[i,,] = Var[x_i] or diag(feature.var[i,]) = Var[x_i]
###     feature.cov[i,]  = Cov[y_i,x_i]
###          weight[i]   = w_i
### NOTE:
###    * Let feature.var[i,,] be a d x d matrix and the length of feature.mean[i,] is n.
###      If n > d, then feature.var correspond to the first d features, and other features
###      have 0 variance.
###    * Let feature.cov[i,] has length d and feature.mean[i,] has length n.
###      If n > d, then feature.cov correspond to the first d features, and other features
###      have 0 covariance.
###
### OUTPUT:
###    out$coeff = beta
###    out$loss  = loss
###    out$num   = sum(w_i > 0)
###    out$var   = var_y = loss / num
###
fit.lm.random.effect <- function(
	response.mean, feature.mean,
	response.var=NULL, feature.var=NULL, feature.cov=NULL, weight=NULL,
	lambda=0 # lambda for ridge regression
){
	nFeatures = ncol(feature.mean);
	if(length(response.mean) != nrow(feature.mean)) stop("length(response.mean) != nrow(feature.mean)");
	if(!is.null(response.var) && length(response.mean) != length(response.var)) stop("length(response.mean) != length(response.var)");
	if(!is.null(feature.cov)  && length(response.mean) != nrow(feature.cov)) stop("length(response.mean) != nrow(feature.cov)");
	if(!is.null(feature.cov)  && ncol(feature.cov) > nFeatures) stop("ncol(feature.cov) > nFeatures");

	if(!is.null(response.var) && any(response.var < 0)) stop("Some response.var < 0");
	if(!is.null(feature.var)  && any(feature.var < 0))  stop("Some feature.var < 0");
	
	if(is.null(lambda) || is.na(lambda)) lambda = 0;
	
	if(!is.null(weight)){
		if(length(weight) != length(response.mean)) stop("length(weight) != nrow(response.mean)");
	}else{
		weight = rep(1, length(response.mean));
	}
	if(any(weight < 0)) stop("any(weight < 0)");

	Var.X = matrix(0, nrow=nFeatures, ncol=nFeatures);
	if(is.null(feature.var)){
		# do nothing
	}else if(is.vector(feature.var)){
		if(length(feature.var) != nrow(feature.mean)) stop("length(feature.var) != nrow(feature.mean)");
		Var.X = sum(feature.var * weight) * diag(nFeatures);
	}else if(length(dim(feature.var)) == 2){
		if(nrow(feature.var) != nrow(feature.mean)) stop("nrow(feature.var) != nrow(feature.mean)");
		if(ncol(feature.var) >  ncol(feature.mean)) stop("ncol(feature.var) > ncol(feature.mean)")
		else if(ncol(feature.var) != ncol(feature.cov)) stop("ncol(feature.var) != ncol(feature.cov)");
		v = apply(feature.var * weight, 2, sum);
		d = ncol(feature.var);
		temp = diag(v, nrow=d, ncol=d);
		Var.X[1:d,1:d] = temp;
	}else if(length(dim(feature.var)) == 3){
		if(dim(feature.var)[1] != nrow(feature.mean)) stop("dim(feature.var)[1] != nrow(feature.mean)");
		if(dim(feature.var)[2] >  ncol(feature.mean)) stop("dim(feature.var)[2] > ncol(feature.mean)")
		else if(dim(feature.var)[2] != ncol(feature.cov)) stop("dim(feature.var)[2] != ncol(feature.cov)");
		if(dim(feature.var)[3] != dim(feature.var)[2]) stop("dim(feature.var)[3] != dim(feature.var)[2]");
		d = dim(feature.var)[2];
		temp = apply(feature.var * weight, c(2,3), sum);
		Var.X[1:d,1:d] = temp;
	}else stop("feature.var error");
	feature.mean.t.W = t(feature.mean * weight);
	E.X.E.X = feature.mean.t.W %*% feature.mean;
	XX = Var.X + E.X.E.X + lambda * diag(nFeatures);
	
	Cov.Y.X = rep(0, nFeatures);
	if(!is.null(feature.cov)){
		d = ncol(feature.cov);
		temp = apply(feature.cov * weight, 2, sum);
		Cov.Y.X[1:d] = temp;
	}
	E.X.E.Y = feature.mean.t.W %*% response.mean;
	XY = drop(E.X.E.Y + Cov.Y.X);
	
	beta = drop(solve(XX) %*% XY);

	if(is.null(response.var)){
		E.y2 = sum(response.mean^2 * weight);
	}else{
		E.y2 = sum((response.mean^2 + response.var) * weight);
	}
	loss = drop(t(beta) %*% (Var.X + E.X.E.X) %*% beta - 2 * t(beta) %*% XY + E.y2);
	num  = sum(weight > 0);
	var_y = loss / num;
	if(var_y < 0){
		loss.2 = loss.lm.random.effect(
			coeff=beta, response.mean=response.mean, feature.mean=feature.mean, 
			response.var=response.var, feature.var=feature.var, feature.cov=feature.cov, weight=weight
		);
		stop("Variance = ",var_y," < 0 (loss: ",loss," vs ",loss.2,")");
	}
	out = list(coeff=beta, loss=loss, num=num, var=var_y);
	return(out);
}
### loss = sum_i w_i*{ (E[y_i] - beta'E[x_i])^2 + beta'Var[x_i]beta - 2beta'Cov[y_i,x_i] + Var[y_i] }
loss.lm.random.effect <- function(
	coeff, response.mean, feature.mean,
	response.var=NULL, feature.var=NULL, feature.cov=NULL, weight=NULL
){
	nFeatures = ncol(feature.mean);
	nObs      = length(response.mean);
	if(nObs != nrow(feature.mean)) stop("nObs != nrow(feature.mean)");
	if(!is.null(response.var) && nObs != length(response.var)) stop("nObs != length(response.var)");
	if(!is.null(feature.cov)  && nObs != nrow(feature.cov)) stop("nObs != nrow(feature.cov)");
	if(!is.null(feature.cov)  && ncol(feature.cov) > nFeatures) stop("ncol(feature.cov) > nFeatures (",ncol(feature.cov)," vs ",nFeatures,")");
	if(length(coeff) != nFeatures) stop("length(ceoff) != nFeatures");
	
	if(!is.null(weight)){
		if(length(weight) != nObs) stop("length(weight) != nObs");
	}else weight = rep(1, nObs);
	if(any(weight < 0)) stop("any(weight < 0)");
	if(is.null(response.var)) response.var = 0;
	
	err.square = drop((response.mean - feature.mean %*% coeff)^2);
	
	if(is.null(feature.var)){
		beta.V.beta = 0;
	}else if(is.vector(feature.var)){
		if(length(feature.var) != nObs) stop("length(feature.var) != nObs");
		beta.V.beta = feature.var * sum(coeff^2);
	}else if(length(dim(feature.var)) == 2){
		if(nrow(feature.var) != nObs) stop("nrow(feature.var) != nObs");
		if(ncol(feature.var) > nFeatures) stop("ncol(feature.var) > nFeatures");
		d = ncol(feature.var);
		beta.V.beta = drop(coeff[1:d] %*% (t(feature.var) * coeff[1:d]));
	}else if(length(dim(feature.var)) == 3){
		if(dim(feature.var)[1] != nObs) stop("dim(feature.var)[1] != nObs");
		if(dim(feature.var)[2] > nFeatures) stop("dim(feature.var)[2] > nFeatures");
		if(dim(feature.var)[3] != dim(feature.var)[2]) stop("dim(feature.var)[3] != dim(feature.var)[2]");
		d = dim(feature.var)[2];
		temp = rep(coeff[1:d], each=nObs*d);
		V.beta = apply(feature.var*temp, c(1,2), sum);
		beta.V.beta = drop(V.beta %*% coeff[1:d]);
	}else stop("feature.var error");
	
	beta.Cov = 0;
	if(!is.null(feature.cov)) beta.Cov = drop(feature.cov %*% coeff[1:ncol(feature.cov)]);
	loss = sum( weight * (err.square + beta.V.beta - 2*beta.Cov + response.var) );
	return(loss);
}

### E[loglikelihood] 
###	  = -(1/2) * sum_i E[ log(var_y) + (w_i/var_y)*(y_i - beta'x_i)^2 ]
###	  = -(1/2) * ( sum(w_i > 0) * log(var_y) + (1/var_y) * loss ),
###     where loss = sum_i w_i*{ (E[y_i] - beta'E[x_i])^2 + beta'Var[x_i]beta - 2beta'Cov[y_i,x_i] + Var[y_i] }
Eloglik.lm.random.effect <- function(
	coeff, var, response.mean, feature.mean,
	response.var=NULL, feature.var=NULL, feature.cov=NULL, weight=NULL,
	algo=NULL
){
	nObs = length(response.mean);
	if(length(var) != 1){
		if(!is.null(weight)) stop("length(var) != 1");
		sum_log_var = sum(log(var));
		weight = 1/var;
		var = 1;
	}else{
		if(is.null(weight)) weight = rep(1, nObs);
		sum_log_var = sum(weight > 0) * log(var);
	}
	if(any(weight < 0)) stop("any(weight < 0)");
	if(length(weight) != nObs) stop("length(weight) != nObs");
	if(is.null(algo)){
		loss = loss.lm.random.effect(
			coeff=coeff, response.mean=response.mean, feature.mean=feature.mean, 
			response.var=response.var, feature.var=feature.var, feature.cov=feature.cov, weight=weight
		);
	}else{
		if(!is.null(feature.var) || !is.null(feature.cov)) stop("When algo is not NULL, you cannot specify feature.var or feature.cov");
		loss = reg.loss(model=coeff, x=feature.mean, y=response.mean, algo=algo, y.var=response.var, weight=weight);
	}
	loglik = -(1/2)*(sum_log_var + loss / var);
	return(loglik);
}

###
### model = NA would predict 0 for any input
###
reg.predict <- function(model, x, algo, ncol=NULL){
	if(is.null(nrow(x))) stop("is.null(nrow(x))");
	if(is.null(ncol)){
		if(length(model) == 1 && is.na(model)){
			return(rep(0.0,nrow(x)));
		}else if(is.list(model) && !is.null(model$name)){
			return(drop(algo$predict(model=model, x=x)));
		}else if(is.vector(model)){
			if(length(model)==0 && ncol(x)==0) return(rep(0.0,nrow(x)));
			return(drop(x %*% model));
		}else{
			cat("str(model) = \n");
			str(model);
			stop("Unknown input model");	
		}
	}else if(is.matrix(model)){
		if(ncol(model) != ncol) stop("ncol(model) != ncol");
		if(ncol(x)==0 && nrow(model)==0) return(array(0.0, dim=c(nrow(x),ncol)));
		return(as.matrix(x %*% model));
	}else{
		if(ncol != length(model)) stop("ncol != length(model)");
		out = matrix(NA, nrow=nrow(x), ncol=ncol);
		for(i in 1:ncol){
			out[,i] = reg.predict(model=model[[i]], x=x, algo=algo, ncol=NULL);
		}
		return(out);
	}
}

###
### Model:
###	   y[i] ~ N(mean = h(x[i,]),  var = var_y/weight[i])
### OUTPUT:
###    out$coeff = h (the regression function)
###    out$loss  = sum_i weight[i] * { (E[y[i]] - h(x[i,]))^2 + Var[y[i]] }
###    out$num   = sum(weight > 0)
###    out$var   = var_y = loss / num
###
reg.train <- function(x, y, algo, y.var=0, weight=NULL, control=NULL, max.nObs=NULL){
	if(length(y) != nrow(x)) stop("length(y) != nrow(x)");
	if(length(y.var) != 1 && length(y.var) != length(y)) stop("length(y.var) != length(y)");
	if(!is.null(weight) && any(weight == 0)){
		select = weight != 0;
		y = y[select];
		x = x[select,,drop=FALSE];
	}
	if(!is.null(max.nObs) && length(y) > max.nObs){
		select = sample(length(y), max.nObs);
		y = y[select];
		x = x[select,,drop=FALSE];
	}
	if(is.null(control)) control = algo$control;
	if(is.null(weight))  weight=rep(1,length(y));
	if(any(weight <= 0)) stop("any(weight <= 0)");
	model = algo$train(x=x, y=y, weight=weight, control=control);
	pred  = drop(algo$predict(model=model, x=x));
	loss  = sum( weight * ((y - pred)^2 + y.var) );
	out   = list(coeff=model, loss=loss, num=length(y), var=loss/length(y));
	return(out);
}

reg.loss <- function(model, x, y, algo, y.var=0, weight=NULL){
	pred  = reg.predict(model=model, x=x, algo=algo);
	if(is.null(weight)) weight = 1;
	loss  = sum( weight * ((y - pred)^2 + y.var) );
	return(loss);
}

###
### Gaussian likelihood function
###
loglik.gaussian <- function(pred.x, x, var_x){
	if(length(var_x) != 1) stop("length(var_x) != 1");
	loglik = -(1/2) * ( sum((x - pred.x)^2 / var_x) + length(x) * log(var_x) );
	return(loglik);
}

###
### Get feature matrix
###	The input data x is a data frame in one of the two formats:
### (1) Dense format:
###     x = data.frame(id, feature1, feature2, ...)
### (2) Sparse format:
###     x = data.frame(id, index, value)
###           id is the case ID
###           index is the feature index
###           value is the feature value
### Output a matrix (may be in the sparseMatrix format), in which
### the nth row correspond to selected.id[n]
get.feature.matrix <- function(
	x, # input data
	id.colname,  # name of the ID column in x
	selected.id, # a list of ID to select
	add.intercept=FALSE, # whether to add a column of all ones
	err.prefix="", err.x.name=NULL, err.select.name=NULL
){
	if(add.intercept){ regFormula = formula(~.);   default = 1.0; }
	else{              regFormula = formula(~.-1); default = 0.0; }
	
	if(is.null(x)) return(matrix(default, nrow=length(selected.id), ncol=1));
	
	x.name = if(is.null(err.x.name)) "input table x" else err.x.name;
	if(!(id.colname %in% names(x))) stop(err.prefix,"Cannot find ID column '", id.colname,"' in ",x.name);
	if(ncol(x) == 3 && all(c(id.colname, "index", "value") %in% names(x))){
		# Sparse format
		nCases    = length(selected.id);
		nFeatures = max(x$index);
		x$row = match(x[[id.colname]], selected.id);
		x = x[!is.na(x$row),];
		if(add.intercept){
			nFeatures = nFeatures + 1;
			intercept = data.frame(row=1:nCases, index=nFeatures, value=1.0);
			x = rbind(x[,c("row","index","value")], intercept);
		}
		out = sparseMatrix(i=x$row, j=x$index, x=x$value, dims=c(nCases, nFeatures));
	}else{
		# Dense format
		select = match(selected.id, x[[id.colname]]);
		if(any(is.na(select))){
			temp = if(is.null(err.select.name)) "" else paste(" in ",err.select.name,sep="");
			stop(err.prefix,"Some IDs",temp," cannot be found in ",x.name,"$",id.colname);
		} 
		out = model.matrix(regFormula, x[select,-match(id.colname, names(x)),drop=FALSE]);
	}
	return(out);
}

is.sparse.feature <- function(x){
	if(is.data.frame(x) && ncol(x) == 3 && all(c("index", "value") %in% names(x))) return(TRUE);
	return(FALSE);
}

###
### For non-Gaussian response
### 
### Initialize the observation table
###   Observation table: obs = data.frame(user, item, y)
###	  Add the following to the code before the EM procedure.
###      obs = init.obs(obs, is.logistic);
###   After the call, obs$response will contain the binary response.
###   Note obs$y will be used to store the Gaussian response values
###   for variational approximation
init.obs <- function(obs, is.logistic){
	if(is.null(obs)) return(NULL);
	obs$response = obs$y; # Make a copy of the original response
	                      # since obs$y may be changed for variational approx.
	if(is.logistic){
		labels = unique(obs$response);
		if(length(labels) != 2) stop("The response is not binary: ",paste(labels[1:min(10,length(labels))],collapse=", "));
		labels = sort(labels);
		if(any(labels != c(0,1)) && any(labels != c(-1,1))) stop("Binary response must be {0,1} or {-1,1}");
	}
	return(obs);
}
### Generate Gaussain response
### Add the following code in the EM loop before the E-step
###    response = generate.response(obs=obs, param=param, is.logistic=is.logistic, verbose=verbose);
###    obs$y = response$y;
###    param$var_y = response$var_y;
generate.response <- function(obs, param, is.logistic, verbose=0){
	if(is.logistic){
		# variational approximation
		if(verbose >= 1) cat("generate gaussian response for logistic\n");
		if(length(param$xi) != nrow(obs)) stop("length(param$xi) != nObs");
		if(length(obs$response) != nrow(obs)) stop("length(param$xi) != nObs");
		response = obs$response;
		if(all(response %in% c(0,1))) response[response == 0] = -1;
		if(!all(response %in% c(-1,1))) stop("Binary response must be {-1,1} at this point");
		var_y = 1/(2 * logistic.lambda(param$xi));
		y = response * var_y / 2;
	}else{
		var_y = param$var_y;
		y     = obs$y;
	}
	return(list(y=y, var_y=var_y));
}
### Update the parameters
### Add the following code in the EM loop after the M-step
###    param = update.param(factor.mean=mc_e$mean, param=param.new, obs=obs, is.logistic=is.logistic);
update.param <- function(factor.mean, param, obs, is.logistic, factor.var=NULL, x_obs=NULL){
	if(is.logistic){
		if(!is.null(factor.mean$pred.y.square)){
			if(length(factor.mean$pred.y.square) != nrow(obs)) stop("length(factor.mean$pred.y.square) != nrow(obs)");
			param$xi = sqrt(factor.mean$pred.y.square);
		}else if(!is.null(factor.mean$fScore) && !is.null(factor.var$fScore)){
			if(length(factor.mean$fScore) != nrow(obs)) stop("length(factor.mean$fScore) != nrow(obs)");
			if(length(factor.var$fScore) != nrow(obs)) stop("length(factor.var$fScore) != nrow(obs)");
			xb = reg.predict(model=param$b, x=x_obs, algo=param$reg.algo);
			pred = xb + factor.mean$fScore;
			param$xi = sqrt(pred^2 + factor.var$fScore);
		}else stop("Cannot compute xi for variational logistic");
	}
	return(param);
}
# input pred.y is based on the Gaussian model
obsLoglik.from.gaussian <- function(pred.y, y, var_y, is.logistic){
	if(length(y) != length(pred.y)) stop("length(y) != length(pred.y)");
	if(is.logistic){
		if(all(y %in% c(0,1))) y[y == 0] = -1;
		if(!all(y %in% c(-1,1))) stop("Binary response must be {-1,1} at this point");
		loglik = sum( -log1p(exp(-y * pred.y) ) );
		attr(loglik, "loss") = -loglik / length(y);
	}else{
		loglik = loglik.gaussian(pred.x=pred.y, x=y, var_x=var_y);
		attr(loglik, "loss") = sqrt(mean( (y - pred.y)^2 ));
	}
	return(loglik);
}
# input pred.y is based on the Gaussian model
# output$pred.y is the input for the Gaussian model
#               is the predicted probability for the Logistic model
predict.response.from.gaussian <- function(pred.y, y, param, is.logistic){
	if(!is.null(y)) loglik = obsLoglik.from.gaussian(pred.y=pred.y, y=y, var_y=param$var_y, is.logistic=is.logistic)
	else            loglik = NULL;
	if(is.logistic){
		pred.y = 1/(1+exp(-pred.y));
		if(!is.null(y)){
			if(min(y) == -1) temp.y = 2*pred.y - 1
			else             temp.y = pred.y;
		}
	}else{
		temp.y = pred.y;
	}
	if(!is.null(y)){
		rmse = sqrt(mean( (y - temp.y)^2 ));
		mae  = mean( abs(y - temp.y) );
	}else{
		rmse = NULL;
		mae  = NULL;
	}
	return(list(pred.y=pred.y, true.y=y, rmse=rmse, mae=mae, loglik=loglik, test.loss=attr(loglik, "loss")));
}
logistic.lambda <- function(xi){
	return(tanh(xi/2) / (4*xi));
}

###
### Output to the output directory
###	IDs: a list of userIDs, itemIDs, etc (can be null)
### prediction$test.loss must exist if prediction is not null
###	loglik is the complete data training log likelihood
###
output.to.dir <- function(
	out.dir, factor, param, IDs, prediction, loglik, 
	minTestLoss, nSamples, iter, out.level, out.overwrite, 
	TimeEStep, TimeMStep, TimeTest, verbose,
	other=NULL, name="est", data.train=NULL
){
	if(!is.null(prediction) && is.null(prediction$test.loss)) stop("prediction$test.loss does not exist");
	if(!is.null(prediction) && is.null(prediction$rmse))      stop("prediction$rmse does not exist");
	if(is.null(attr(loglik, "loss"))) stop("is.null(attr(loglik, 'loss'))");
	
	if(out.level <= 0) return(NULL);
	b.time.write = proc.time();
	if(is.null(out.dir)) stop("Please specify out.dir");
	if(iter == 0){
		if(file.exists(paste(out.dir,"/",name,".last",sep="")) && !out.overwrite){
			stop("Output File '",out.dir,"' EXISTS!!");
		}else if(!file.exists(out.dir)){
			dir.create(out.dir, recursive=TRUE);
		}
		smry.file = paste(out.dir,"/summary",sep="");
		if(file.exists(smry.file)) file.remove(smry.file);
	}
	
	thisTestLoss = if(is.null(prediction)) -1 else prediction$test.loss;
	TestRMSE     = if(is.null(prediction)) -1 else prediction$rmse;
	
	if(iter == 0){
		if(out.level >= 2) save(file=paste(out.dir,"/",name,".0",sep=""), list=c("factor", "param"));
	}else{
		file = paste(out.dir,"/",name,".last",sep="");
		if(out.level >= 2 && iter >= 2){
			file.prev = paste(out.dir,"/",name,".",(iter-1),sep="");
			if(file.exists(file.prev)) file.remove(file.prev);
			file.rename(file, file.prev);
		}
		save(file=file, list=c("factor", "param", "prediction", "IDs", "other", "data.train"));
		if(!is.null(prediction)){
			if(thisTestLoss == minTestLoss) file.copy(file, paste(out.dir,"/",name,".minTestLoss",sep=""), overwrite=TRUE);
		}
		save(file=paste(out.dir,"/param.",iter,sep=""), list=c("param"));
	}
	if(length(nSamples) != 1){
		nSamples = if(iter > 0) nSamples[iter] else 0;
	}
	
	summary = data.frame(Method="MCEM", Iter=iter, nSteps=nSamples, CDlogL=loglik, TestLoss=thisTestLoss, LossInTrain=attr(loglik,"loss"), TestRMSE=TestRMSE, TimeEStep=TimeEStep, TimeMStep=TimeMStep, TimeTest=TimeTest);
	file = paste(out.dir,"/summary",sep="");
	if(file.exists(file)) write.table(summary, file=file, append=TRUE,  quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	else                  write.table(summary, file=file, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE);
	if(verbose > 0){
		time.used.write = proc.time() - b.time.write;
		cat("write a model & summary info on to disk (used ",time.used.write[3]," sec)\n",sep="");
	}
}
