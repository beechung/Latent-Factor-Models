### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### EM hierarchical smoothing: 1-dimensional, 2 levels
###
###	INPUT:
###	 (1) data = data.frame(item, category, obs, var);
###
### 		(item[n], category[n], obs[n], var[n])
###               i            k
###      Each tuple specifies that we have observation obs[n] of item i in category k
###      with observation variance proportional to var[n].
###      An item can be in multiple categories.
### 
###	     Consider different user segments interacting with different items.
###		 category[n] would represent the user segment of the nth observation.
###
###		 IMPORTANT NOTE: Item and category indices start from 1 (NOT 0)
###
###  (2) feature = list(x_obs, x_b); (optional)
###      x_obs[n,j]: The jth feature of the nth observation.
###                  x_obs[n,] correspond to data[n,].
###                  This is a matrix of size: #Observations x #ObsFeatures
###      x_b[i,j]:   The jth feature of item i
###                  This is a matrix of size: #Items x #ItemFeatures
###
###      To include an intercept, add of a column of all 1s to x_obs and x_b.
###
###  (3) offset (optional)
###
### MODEL:
### 	 obs[n] ~ N(mean = b[i,k] + sum(w_obs * x_obs[n,]) + offset[n],  var=var[n]*var_obs_adj)
###      b[i,k] ~ N(mean = q[k]*a[i] + sum(w_b[k,] * x_b[i,]),           var=var_b[k])
###      a[i]   ~ N(mean = 0,                                            var=var_a)
###
###    where b[i,k] is the factor for item i in category k;
###          a[i]   is the global factor for item i
###          w_obs  is the regression weight vector for the observation features.
###          w_b    is the regression weight vector for the item features.
###
### OUTPUT:
###    Factors:  factor  = list(a=list(mean,var), b=list(mean,var,cov));
###	Parameters:  param   = list(var_obs_adj, var_b, var_a, q, w_obs, w_b);
###
###	Observations can be pre-aggregated.  If there are multiple observations
###     of (item=i, category=k), then it is equivalent to replace all these
###     observations by a single observation.
###     For simplicity, assume there are no observation features and no offsets, and
###           all those observations have the same variance sigma^2.
###		Then, the new observation is the sample mean of all those observations, and
###           the new variance is sigma^2 / sample_size
###
###     This is based on the following equivalence:
###		          obs[n] ~ Normal(mu, sigma^2), for n=1,...,N
###     (sum_n obs[n])/N ~ Normal(mu, sigma^2/N)
###
###

###
### Simulate some data
###
nItems=1000; nRep=3;  # 1000 items, each has 11 observations
nCategories=13;        # 13 user segments
nObsFeatures=0; nItemFeatures=2;  # no observation features, two item features
source("src/R/model/hierarchical_genData.R");
set.seed(1);
d = genData.1D.2Levels(
		nItems=nItems, nCategories=nCategories, nRep=nRep, nObsFeatures, nItemFeatures, offset.sd=0,
		q.mean=1, q.sd=3, w_obs.sd=1, w_b.sd=1, var_a=5, var_b.mean=1, var_b.sdlog=0.2, var_obs.mean=0.2, var_obs.sdlog=0.1
);
# Check the data format
str(d); # d$data, d$feature and d$offset are input to the fitting algorithm.
        # Others are the ground truth used to generate the data.
# d$feature$x_obs = matrix(1,nrow=nrow(d$data),ncol=1); # (uncomment this you want to add a global intercept)
# Randomly split data into training and test if desired
split = split.each.item.1D.2Levels(problem=d, train.frac=0.75, by.category=FALSE);

###
### Use EM to fit the model based on the training data
###
dyn.load("lib/c_funcs.so");
source("src/R/util.R")
source("src/R/model/hierarchical_utils.R");
set.seed(0);
model = fit.hierModel.1Dim2Levels(
		data=split$train$data, lambda=1e-3, nEM.iter=10, feature=split$train$feature, offset=split$train$offset,
		data.test=split$test$data, feature.test=split$test$feature, offset.test=split$test$offset, # optional
		init.var=1, init.q=1, init.factor.a.sd=1, init.factor.b.sd=1, # optional
		verbose=1, debug=0
);
# check the output
str(model);

###
### Smooth the observations based on the fitted parameters.
###
factor.test = smooth.hierModel.1Dim2Levels(
	data=split$test$data, param=model$param, feature=split$test$feature, offset=split$test$offset
);
# make prediction using the smoothed factors
pred = predict.hierModel.1Dim2Levels(
	data=split$test$data, factor=factor.test, param=model$param, feature=split$test$feature, offset=split$test$offset
);
pred$rmse;


###
### Some baseline methods
###
baseline = fit.a.only.1Dim2Levels(data=split$train$data, offset=split$train$offset);
predict.hierModel.1Dim2Levels(data=split$test$data, factor=baseline, offset=split$test$offset)$rmse;

baseline = fit.b.only.1Dim2Levels(data=split$train$data, offset=split$train$offset);
predict.hierModel.1Dim2Levels(data=split$test$data, factor=baseline, offset=split$test$offset)$rmse;

