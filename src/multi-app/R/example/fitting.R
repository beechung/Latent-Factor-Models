### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### Run the EM algorithm
###
###		Working directory: src/multi-app
###     To compile, just type make
###
dyn.load("c_funcs.so");
source("R/utils.R");
source("R/model/util.R");

# Simulate some data
source("R/model/generate-data.R");
set.seed(0);
data = generate.GaussianData(
		nUsers=200, nApps=3,
		nGlobalFactors = 2,
		nFeatures     = c(2,   4,   3),
		nItems        = c(5,   2,   3),
		nLocalFactors = c(2,   3,   4),
		frac.missing  = c(0.1, 0.3, 0.2),
		var_x         = c(0.05,0.1, 0.15), 
		var_y         = c(0.05,0.2, 0.4),
		var_z         = c(2,   1,   3),
		var_u=1, A.sd=1, B.sd=1, b.sd=1, alpha.sd=2, beta.sd=2,
		w_x.mean=1, w_x.var=0.01, w_y.mean=1, w_y.var=0.01
);
size = check.syntax.all(feature=data$feature, response=data$response, param=data$param, factor=data$factor, check.indices=TRUE)

# Create random train/test split
select.train = runif(nrow(data$feature)) <= 0.75
feature.train = data$feature[ select.train,];
feature.test  = data$feature[!select.train,];
select.train = runif(nrow(data$response)) <= 0.75
response.train = data$response[ select.train,];
response.test  = data$response[!select.train,];

# Fit the model using EM
dyn.load("c_funcs.so");
source("R/utils.R");
source("R/model/util.R");
source("R/model/fit-EM.R");
set.seed(1);
ans = fit.EM(
		# Input data
		feature=feature.train,   # data.frame(user, app, index, x, w)
		response=response.train, # data.frame(user, app, item,  y, w)
		# Model setup
		nGlobalFactors=2,       # num of global factors per user
		nLocalFactors=c(2,3,4), # num of local  factors per user (may be a vector of length: #app)
		identity.A=FALSE, identity.B=FALSE, # whether A/B is a identity matrix
		# Test data (optional)
		test.feature=feature.test,
		test.response=response.test,
		# Model-fitting parameters
		nIter=20, # num of EM iterations
		ridge.lambda=c(A=0, B=0, beta=0), # lambda values for ridge regression for A, B, beta
		keep.users.without.obs.in.E.step=FALSE,
		keep.users.without.obs.in.M.step=FALSE,
		# Initialization parameters (optional)
		param=NULL, # directly set all the parameters (if not NULL, the following will not be used)
		var_x=1, var_y=1, var_z=1, # initial var (may be vectors of length: #app)
		var_u=1,                   # initial var (length: 1)
		A.sd=1, B.sd=1, beta.sd=1, # A_k ~ N(mean=0, sd=A.sd), and so on.
		# Output options
		out.level=1,  # out.level=1: Save the model in out.dir/model.last and out.dir/model.minTestLoss
		out.dir="/tmp/test-abcde", # out.level=2: Save the model after each iteration i in out.dir/model.i
		out.overwrite=TRUE, # whether to allow overwriting existing files
		# Debugging options (optional)
		debug=0, verbose=2, use.C=TRUE,
		show.marginal.loglik=TRUE  # Very computationally expensive if TRUE (for debugging only)
);
# Now, ans$model is the fitted model
#      The fitted model is also stored in /tmp/test-abcde/model.last (as specified using out.dir)
#      A summary of the EM iterations is in /tmp/test-abcde/summary.txt
load("/tmp/test-abcde/model.last");
str(model); # This is the fitted model
read.table("/tmp/test-abcde/summary.txt", header=TRUE, sep="\t");

# Check whether the marginal log likelihood keeps increasing
#  ans$marginal.loglik[k+1]: the log likelihood after the kth iteration
ans$marginal.loglik

# Marginal likelihood of the ground-truth model
marginal.loglik.R(feature=data$feature, response=data$response, param=data$param);

# Check the test-set loss (RMSE for the gaussian model)
# Test-set loss of the response
ans$test.loss.y
# Test-set loss of the features
ans$test.loss.x

# Make predictions
model = ans$model
pred = predict.x.and.y(feature=feature.test, response=response.test, param=model$param, factor=model$factor);
# Now pred$pred.y contains the predicted response values (pred$pred.y[m] corresponds to response.test[m,])
#     pred$pred.x contains the predicted feature  values (pred$pred.x[m] corresponds to  feature.test[m,])


###
### Fit a model with no features
###
feature.train.empty = feature.train[0,]
feature.test.empty  = feature.test[0,]

dyn.load("c_funcs.so");
source("R/utils.R");
source("R/model/util.R");
source("R/model/fit-EM.R");
set.seed(1);
ans = fit.EM(
		# Input data
		feature=feature.train.empty, # data.frame(user, app, index, x, w)
		response=response.train,     # data.frame(user, app, item,  y, w)
		# Model setup
		nGlobalFactors=2,       # num of global factors per user
		nLocalFactors=c(2,3,4), # num of local  factors per user (may be a vector of length: #app)
		identity.A=FALSE, identity.B=FALSE, # whether A/B is a identity matrix
		# Test data (optional)
		test.feature=feature.test.empty,
		test.response=response.test,
		# Model-fitting parameters
		nIter=20, # num of EM iterations
		ridge.lambda=c(A=0, B=0, beta=0), # lambda values for ridge regression for A, B, beta
		keep.users.without.obs.in.E.step=FALSE,
		keep.users.without.obs.in.M.step=FALSE,
		# Initialization parameters (optional)
		param=NULL, # directly set all the parameters (if not NULL, the following will not be used)
		var_x=1, var_y=1, var_z=1, # initial var (may be vectors of length: #app)
		var_u=1,                   # initial var (length: 1)
		A.sd=1, B.sd=1, beta.sd=1, # A_k ~ N(mean=0, sd=A.sd), and so on.
		# Output options
		out.level=1,  # out.level=1: Save the model in out.dir/model.last and out.dir/model.minTestLoss
		out.dir="/tmp/test-abcde", # out.level=2: Save the model after each iteration i in out.dir/model.i
		out.overwrite=TRUE, # whether to allow overwriting existing files
		# Debugging options (optional)
		debug=0, verbose=2, use.C=TRUE,
		show.marginal.loglik=TRUE  # Very computationally expensive if TRUE (for debugging only)
);

# Make predictions
model = ans$model
pred = predict.x.and.y(feature=NULL, response=response.test, param=model$param, factor=model$factor);
# Now pred$pred.y contains the predicted response values (pred$pred.y[m] corresponds to response.test[m,])
#     pred$pred.x contains the predicted feature  values (pred$pred.x[m] corresponds to  feature.test[m,])


###
### Fit a model with no features and a single context
###
feature.train.empty = feature.train[0,]
feature.test.empty  = feature.test[0,]
response.train$app = as.integer(1);
response.test$app  = as.integer(1);

dyn.load("c_funcs.so");
source("R/utils.R");
source("R/model/util.R");
source("R/model/fit-EM.R");
set.seed(1);
ans = fit.EM(
		# Input data
		feature=feature.train.empty, # data.frame(user, app, index, x, w)
		response=response.train,     # data.frame(user, app, item,  y, w)
		# Model setup
		nGlobalFactors=5,       # num of global factors per user
		nLocalFactors=5, # num of local  factors per user (may be a vector of length: #app)
		identity.A=TRUE, identity.B=FALSE, # whether A/B is a identity matrix
		# Test data (optional)
		test.feature=feature.test.empty,
		test.response=response.test,
		# Model-fitting parameters
		nIter=20, # num of EM iterations
		ridge.lambda=c(A=0, B=0, beta=0), # lambda values for ridge regression for A, B, beta
		keep.users.without.obs.in.E.step=FALSE,
		keep.users.without.obs.in.M.step=FALSE,
		fix.var_u=1e-6,
		# Initialization parameters (optional)
		param=NULL, # directly set all the parameters (if not NULL, the following will not be used)
		var_x=1, var_y=1, var_z=1, # initial var (may be vectors of length: #app)
		var_u=1,                   # initial var (length: 1)
		A.sd=1, B.sd=1, beta.sd=1, # A_k ~ N(mean=0, sd=A.sd), and so on.
		# Output options
		out.level=1,  # out.level=1: Save the model in out.dir/model.last and out.dir/model.minTestLoss
		out.dir="/tmp/test-abcde", # out.level=2: Save the model after each iteration i in out.dir/model.i
		out.overwrite=TRUE, # whether to allow overwriting existing files
		# Debugging options (optional)
		debug=0, verbose=2, use.C=TRUE,
		show.marginal.loglik=TRUE  # Very computationally expensive if TRUE (for debugging only)
);

# Make predictions
model = ans$model
pred = predict.x.and.y(feature=NULL, response=response.test, param=model$param, factor=model$factor);
# Now pred$pred.y contains the predicted response values (pred$pred.y[m] corresponds to response.test[m,])
#     pred$pred.x contains the predicted feature  values (pred$pred.x[m] corresponds to  feature.test[m,])

#####################################################################
### Fit Logistic Model
#####################################################################
is.y.logistic=TRUE; # whether to use a logistic model for y
is.x.logistic=FALSE; # whether to use a logistic model for x

dyn.load("c_funcs.so");
source("R/utils.R");
source("R/model/util.R");

# Simulate some data
source("R/model/generate-data.R");
set.seed(0);
data = generate.GaussianData(
		nUsers=4009, nApps=3,
		nGlobalFactors = 2,
		nFeatures     = c(2,   4,   3),
		nItems        = c(5,   2,   3),
		nLocalFactors = c(2,   3,   4),
		frac.missing  = c(0.1, 0.3, 0.2),
		var_x         = c(0.05,0.1, 0.15), 
		var_y         = c(0.05,0.2, 0.4),
		var_z         = c(0.2, 0.1, 0.3),
		is.y.logistic=is.y.logistic,
		is.x.logistic=is.x.logistic,
		y.bias=-3.5,
		var_u=1, A.sd=1, B.sd=1, b.sd=1, alpha.sd=2, beta.sd=1,
		w_x.mean=1, w_x.var=0.01, w_y.mean=1, w_y.var=0.01
);
mean(data$response$y);
# hist(data$y.prob);
size = check.syntax.all(feature=data$feature, response=data$response, param=data$param, factor=data$factor, check.indices=TRUE)

# Create random train/test split
select.train = runif(nrow(data$feature)) <= 0.75
feature.train = data$feature[ select.train,];
feature.test  = data$feature[!select.train,];
select.train = runif(nrow(data$response)) <= 0.75
response.train = data$response[ select.train,];
response.test  = data$response[!select.train,];
y.prob.test    = data$y.prob[!select.train]

# Fit the logistic model using EM
dyn.load("c_funcs.so");
source("R/utils.R");
source("R/model/util.R");
source("R/model/fit-EM.R");
set.seed(1);
ans = fit.EM(
		# Input data
		feature=feature.train,   # data.frame(user, app, index, x, w)
		response=response.train, # data.frame(user, app, item,  y, w)
		# Model setup
		nGlobalFactors=2,       # num of global factors per user
		nLocalFactors=c(2,3,4), # num of local  factors per user (may be a vector of length: #app)
		identity.A=FALSE, identity.B=FALSE, # whether A/B is a identity matrix
		# Test data (optional)
		test.feature=feature.test,
		test.response=response.test,
		# Model-fitting parameters
		nIter=20, # num of EM iterations
		ridge.lambda=c(A=0, B=0, beta=0), # lambda values for ridge regression for A, B, beta
		keep.users.without.obs.in.E.step=FALSE,
		keep.users.without.obs.in.M.step=FALSE,
		is.y.logistic=is.y.logistic,
		is.x.logistic=is.x.logistic,
		# Initialization parameters (optional)
		param=NULL, # directly set all the parameters (if not NULL, the following will not be used)
		var_x=1, var_y=1, var_z=1, # initial var (may be vectors of length: #app)
		var_u=1,                   # initial var (length: 1)
		A.sd=1, B.sd=1, beta.sd=1, # A_k ~ N(mean=0, sd=A.sd), and so on.
		# Output options
		out.level=1,  # out.level=1: Save the model in out.dir/model.last and out.dir/model.minTestLoss
		out.dir="/tmp/test-abcde", # out.level=2: Save the model after each iteration i in out.dir/model.i
		out.overwrite=TRUE, # whether to allow overwriting existing files
		# Debugging options (optional)
		debug=0, verbose=2, use.C=TRUE,
		show.marginal.loglik=FALSE  # Very computationally expensive if TRUE (for debugging only)
);
# Now, ans$model is the fitted model
#      The fitted model is also stored in /tmp/test-abcde/model.last (as specified using out.dir)
#      A summary of the EM iterations is in /tmp/test-abcde/summary.txt
load("/tmp/test-abcde/model.last");
str(model); # This is the fitted model
read.table("/tmp/test-abcde/summary.txt", header=TRUE, sep="\t");

# Check the test-set loss (RMSE for the gaussian model, average -log(likelihood) for logistic)
# Test-set loss of the response
ans$test.loss.y
# Test-set loss of the features
ans$test.loss.x

# Make predictions
model = ans$model
pred = predict.x.and.y(feature=feature.test, response=response.test, param=model$param, factor=model$factor);
# Now pred$pred.y contains the predicted response values (pred$pred.y[m] corresponds to response.test[m,])
#     pred$pred.x contains the predicted feature  values (pred$pred.x[m] corresponds to  feature.test[m,])

# Fit the Gaussian model
dyn.load("c_funcs.so");
source("R/utils.R");
source("R/model/util.R");
source("R/model/fit-EM.R");
set.seed(1);
ans.G = fit.EM(
		# Input data
		feature=feature.train,   # data.frame(user, app, index, x, w)
		response=response.train, # data.frame(user, app, item,  y, w)
		# Model setup
		nGlobalFactors=2,       # num of global factors per user
		nLocalFactors=c(2,3,4), # num of local  factors per user (may be a vector of length: #app)
		identity.A=FALSE, identity.B=FALSE, # whether A/B is a identity matrix
		# Test data (optional)
		test.feature=feature.test,
		test.response=response.test,
		# Model-fitting parameters
		nIter=20, # num of EM iterations
		ridge.lambda=c(A=0, B=0, beta=0), # lambda values for ridge regression for A, B, beta
		keep.users.without.obs.in.E.step=FALSE,
		keep.users.without.obs.in.M.step=FALSE,
		# Initialization parameters (optional)
		param=NULL, # directly set all the parameters (if not NULL, the following will not be used)
		var_x=1, var_y=1, var_z=1, # initial var (may be vectors of length: #app)
		var_u=1,                   # initial var (length: 1)
		A.sd=1, B.sd=1, beta.sd=1, # A_k ~ N(mean=0, sd=A.sd), and so on.
		# Output options
		out.level=1,  # out.level=1: Save the model in out.dir/model.last and out.dir/model.minTestLoss
		out.dir="/tmp/test-abcde", # out.level=2: Save the model after each iteration i in out.dir/model.i
		out.overwrite=TRUE, # whether to allow overwriting existing files
		# Debugging options (optional)
		debug=0, verbose=2, use.C=TRUE,
		show.marginal.loglik=FALSE  # Very computationally expensive if TRUE (for debugging only)
);
pred.G = predict.x.and.y(feature=feature.test, response=response.test, param=ans.G$model$param, factor=ans.G$model$factor);


# Plot ROC curves
data = data.frame(truth=y.prob.test, logistic=pred$pred.y, gaussian=pred.G$pred.y);
col = c( 1, 2, 3, 4, 6, 8);
lty = c( 1, 2, 1, 5, 3, 4);
library(ROCR);
set.seed(0);
for(i in 1:length(data)){
	name = names(data)[i];
	if(i > 1) par(new=TRUE);
	pred.roc = prediction(data[,name],response.test$y);
	perf = performance(pred.roc,"tpr","fpr");
	pts = subsample.ROC(perf);
	plot(pts, type="l", col=col[i], lty=lty[i]);
}
legend("bottomright", legend=names(data), lty=lty[1:length(data)], col=col[1:length(data)], box.lty=0, cex=0.9);
abline(0,1);

