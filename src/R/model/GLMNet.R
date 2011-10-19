### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### GLM net
###
###		 model$glmnet:  The random forest model (fitted by randomForest)
###      model$control: The control parameters
###      model$name:    The name of the model
###		 model$type:    "cv" or "regular"
###      model$lambda:  if type="cv",      "lambda.1se" or "lambda.min"
###                     if type="regular", the value of lambda

library(glmnet);

GLMNet.generateRndModel <- function(
	featureNames, frac.zero=0.5
){
	nFeatures = length(featureNames);
	X = matrix(rnorm(10*nFeatures*nFeatures), nrow=10*nFeatures, dimnames=list(NULL, featureNames));
	beta = rnorm(nFeatures);
	beta[ceiling((1-frac.zero)*nFeatures):nFeatures] = 0;
	y = drop(X %*% beta);
	model = list();
	model$glmnet = glmnet(x=X,y=y, nlambda=5);
	model$control = list();
	model$name = "GLMNet";
	model$type = "regular";
	model$lambda = model$glmnet$lambda[length(model$glmnet$lambda)];
	
	return(model);
}

GLMNet.train <- function(
	x, y, weight, control=list(type="cv", lambda="lambda.1se", nfolds=4, family="gaussian", alpha=0.5, nlambda=20)
){
	if(any(!(names(control) %in% c("type", "lambda", "nfolds", "family", "alpha", "nlambda")))) stop("Unknown control parameter!");
	if(length(y) != nrow(x)) stop("length(y) != nrow(x)");
	
	model = list(control=control);
	model$name = "GLMNet";
	model$lambda = control$lambda;
	if(control$type == "regular"){
		model$glmnet = glmnet(x=x, y=y, family=control$family, weights=weight, alpha=control$alpha, nlambda=control$nlambda);
		model$type = "regular";
	}else if(control$type == "cv"){
		model$glmnet = cv.glmnet(x=x, y=y, family=control$family, weights=weight, alpha=control$alpha, nlambda=control$nlambda, nfolds=control$nfolds);
		model$type = "cv";
	}else stop("Unknown type: ", control$type);
	
	return(model);
}

GLMNet.predict <- function(
	model, x
){	
	if(model$name != "GLMNet") stop("The input model is not a GLMNet");
	ans = drop(predict(object=model$glmnet, newx=x, s=model$lambda, type="response"));
	return(ans);
}

GLMNet = list(
		name="GLMNet",
		generateRndModel=GLMNet.generateRndModel,
		train=GLMNet.train,
		predict=GLMNet.predict,
		control=list(type="cv", lambda="lambda.1se", nfolds=4, family="gaussian", alpha=0.5, nlambda=20)
);
