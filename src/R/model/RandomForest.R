### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### Random Forest
###
###		 model$rf:      The random forest model (fitted by randomForest)
###      model$control: The control parameters
###      model$name:    The name of the model

library(randomForest);

RandomForest.generateRndModel <- function(
	featureNames, nLeaves=20
){
	nFeatures = length(featureNames);
	X = matrix(rnorm(nLeaves*nFeatures), nrow=nLeaves, dimnames=list(NULL, featureNames));
	y = rnorm(nLeaves);
	model = list();
	model$rf = randomForest(x=X,y=y,ntree=1,mtry=nFeatures,replace=FALSE,sampsize=nLeaves,nodesize=1);
	model$control = list();
	model$name = "randomForest";
	
	return(model);
}

RandomForest.train <- function(
	x, y, weight, control=list(min.sup=0.002, ntree=100)
){
	if(any(!(names(control) %in% c("min.sup", "ntree")))) stop("Unknown control parameter!");
	if(length(y) != nrow(x)) stop("length(y) != nrow(x)");
	
	model = list(control=control);
	model$rf = randomForest(x=x, y=y, ntree=control$ntree, nodesize=max(1,length(y)*control$min.sup));
	model$name = "randomForest";
	
	return(model);
}

RandomForest.predict <- function(
	model, x
){	
	if(model$name != "randomForest") stop("The input model is not a randomForest");
	ans = predict(object=model$rf, newdata=x);
	return(ans);
}

RandomForest = list(
		name="RandomForest",
		generateRndModel=RandomForest.generateRndModel,
		train=RandomForest.train,
		predict=RandomForest.predict,
		control=list(min.sup=0.002, ntree=100)
);
