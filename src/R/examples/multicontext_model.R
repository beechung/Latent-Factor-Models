### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### Preparation:
###    (1) Set your path/alias to run the right version of R
###    (2) make  (in public-factor-models/, not in any subdirectory)
###    (3) Take a look at public-factor-models/src/R/model/Notation-multicontext.txt
###        for the specification of the model.
###    (4) Run R (in public-factor-models/, not in any subdirectory)
###

###
### Example 1: Run the fitting code with synthetic data
###
# (1) Generate some data
#     See src/R/model/multicontext_model_genData.R for details
library(Matrix);
dyn.load("lib/c_funcs.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_genData.R");
source("src/R/model/multicontext_model_utils.R");
set.seed(0);
d = generate.GaussianData(
		nSrcNodes=203, nDstNodes=203, nObs=10003, 
		nSrcContexts=4, nDstContexts=5, nEdgeContexts=3, nFactors=2, has.gamma=FALSE, has.u=TRUE,
		nObsFeatures=2, nSrcFeatures=3, nDstFeatures=3, nCtxFeatures=1,
		b.sd=1, g0.sd=1, d0.sd=1, h0.sd=0, G.sd=1, D.sd=1, H.sd=0, q.sd=1, r.sd=1,
		q.mean=5, r.mean=5,
		var_y=0.1, var_alpha=0.5, var_beta=0.5, var_gamma=1, var_v=1, var_u=1, var_w=1,
		var_alpha_global=0.2, var_beta_global=0.2,
		has.intercept=FALSE,
		sparse.matrices=TRUE, frac.zeroFeatures=0.2
);
# (2) Create training/test split
select.train = runif(nrow(d$obs),min=0,max=1) < 0.75;
obs = d$obs;  names(obs) = c("src_id", "dst_id", "src_context", "dst_context", "ctx_id", "y");
obs.train = obs[ select.train,];  x_obs.train = data.frame(as.matrix(d$feature$x_obs)[ select.train,,drop=FALSE]);
obs.test  = obs[!select.train,];  x_obs.test  = data.frame(as.matrix(d$feature$x_obs)[!select.train,,drop=FALSE]);
x_src = data.frame(src_id=1:nrow(d$feature$x_src), as.matrix(d$feature$x_src));
x_dst = data.frame(dst_id=1:nrow(d$feature$x_dst), as.matrix(d$feature$x_dst));
# The following are input data tables:
#     obs.train, obs.test, x_src, x_dst
# obs.train and obs.test contain the training and test rating data
#      The columns of these two tables are:
#      1. src_id: e.g., user_id
#      2. dst_id: e.g., item_id
#      3. src_context: (optional) This is the context in which the source node gives the rating
#      4. dst_context: (optional) This is the context in which the destination node receives the rating
#      5. ctx_id:      (optional) This is the context of this (src_id, dst_id) pair
#      6. y: This is the rating that the source node gives the destination node
# Note: You may set all/any of src_context, dst_context, ctx_id to NULL if there is no context info
#            or set all of them to the same vector
#       The number of contexts cannot be too many; otherwise, the program will be very slow
str(obs.train); # to see the data structure
str(obs.test);  # to see the data structure
# x_src is the source node (e.g., user) feature table
#       The first column src_id specifies the source node ID
str(x_src); # to see the data structure
# x_dst is the destination node (e.g., item) feature table
#       The first column dst_id specifies the destination node ID
str(x_dst); # to see the data structure

# (3) Index training data
#     See src/R/model/multicontext_model_utils.R: indexData() for details
data.train = indexData(
		obs=obs.train, src.dst.same=TRUE, rm.self.link=TRUE,
		x_obs=x_obs.train, x_src=x_src, x_dst=x_dst,
		add.intercept=FALSE,
);
# (4) Index test data
#     See src/R/model/multicontext_model_utils.R: indexTestData() for details
data.test = indexTestData(
		data.train=data.train, obs=obs.test,
		x_obs=x_obs.test, x_src=x_src, x_dst=x_dst,
);
# (5) Setup the model(s) to be fitted
#     See src/R/model/multicontext_model_EM.R: run.multicontext(), fit.multicontext()
#     Note run.multicontext() is a wrapper to fit multiple models using fit.multicontext().
setting = data.frame(
		name          = c("wuv", "wvv"),
		nFactors      = c(    2,     2), # number of interaction factors
		has.u         = c(    T,     F), # whether to use u_i' v_j or v_i' v_j
		has.gamma     = c(    F,     F), # just set to F
		nLocalFactors = c(    0,     0), # just set to 0
		is.logistic   = c(    F,     F)  # whether to use the logistic model for binary rating
);
dyn.load("lib/c_funcs.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_genData.R");
source("src/R/model/multicontext_model_utils.R");
source("src/R/model/multicontext_model_MStep.R");
source("src/R/model/multicontext_model_EM.R");
set.seed(2);
# (6) Run the fitting code
#     See src/R/model/multicontext_model_EM.R: run.multicontext(), fit.multicontext()
#     Note run.multicontext() is a wrapper to fit multiple models using fit.multicontext().
ans = run.multicontext(
		obs=data.train$obs,         # Observation table
		feature=data.train$feature, # Features
		setting=setting,    # Model setting
		nSamples=200,   # Number of samples drawn in each E-step: could be a vector of size nIter.
		nBurnIn=20,     # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
		nIter=20,       # Number of EM iterations
		test.obs=data.test$obs,         # Test data: Observations for testing (optional)
		test.feature=data.test$feature, #            Features for testing     (optional)
		ridge.lambda=1,
		IDs=data.test$IDs,
		out.level=1,         # out.level=1: Save the factor & parameter values to out.dir/model.last and out.dir/model.minTestLoss
		out.dir="/tmp/test", # out.level=2: Save the factor & parameter values of each iteration i to out.dir/model.i
		out.overwrite=TRUE,     # whether to overwrite the output directory if it exists
		debug=0,      # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
		verbose=1,    # Set to 0 to disable console output; Set to 100 to print everything to the console
		verbose.M=2
);
# There may be some warning messages, which are mostly debugging messages and do not mean real problems.

# (7) Checking the model summary
ans$summary[,c("name", "nFactors", "has.u", "has.gamma", "nLocalFactors", "is.logistic", "best.test.loss", "last.test.loss")];

# (8) Load the fitted model(s)
#     Here, I only use the "wuv" model as an example
# (8.1) Check the summary file
read.table("/tmp/test_wuv/summary", header=TRUE, sep="\t", as.is=TRUE);
# (8.2) Load the model
load("/tmp/test_wuv/model.last");
#       Now, factor and param contain the fitted model
str(factor);
str(param);
# (8.3) Make prediction
prediction = predict.multicontext(
	model=list(factor=factor, param=param), 
	obs=data.test$obs, feature=data.test$feature, is.logistic=FALSE
);
# Now, prediction$pred.y contains the predicted rating for data.test$obs
str(prediction);


###
### Example 2: Run the fitting code with synthetic data using SPARSE feature matrix
###
# (1) Generate some data
#     See src/R/model/multicontext_model_genData.R for details
library(Matrix);
dyn.load("lib/c_funcs.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_genData.R");
source("src/R/model/multicontext_model_utils.R");
set.seed(0);
d = generate.GaussianData(
		nSrcNodes=1003, nDstNodes=1003, nObs=100003, 
		nSrcContexts=3, nDstContexts=3, nEdgeContexts=1, nFactors=3, has.gamma=FALSE, has.u=FALSE,
		nObsFeatures=13, nSrcFeatures=19, nDstFeatures=23, nCtxFeatures=1,
		b.sd=1, g0.sd=1, d0.sd=1, h0.sd=0, G.sd=1, D.sd=1, H.sd=0, q.sd=1, r.sd=1,
		q.mean=5, r.mean=5,
		var_y=0.1, var_alpha=0.5, var_beta=0.5, var_gamma=1, var_v=1, var_u=1, var_w=1,
		var_alpha_global=0.2, var_beta_global=0.2,
		has.intercept=FALSE,
		sparse.matrices=TRUE, index.value.format=TRUE, frac.zeroFeatures=0.5
);
names(d$obs) = c("src_id", "dst_id", "src_context", "dst_context", "ctx_id", "y");
d$obs$ctx_id = NULL;
rating.data = d$obs;
x_obs=d$feature$x_obs[order(d$feature$x_obs$row,d$feature$x_obs$col),];  names(x_obs) = c("obs_id", "index", "value");
x_src=d$feature$x_src[order(d$feature$x_src$row,d$feature$x_src$col),];  names(x_src) = c("src_id", "index", "value");
x_dst=d$feature$x_dst[order(d$feature$x_dst$row,d$feature$x_dst$col),];  names(x_dst) = c("dst_id", "index", "value");

#
# Input data: rating.data, x_obs, x_src, x_dst (you need to prepare these four tables for your data)
# Note: All ID numbers start from 1 (not 0)
#
str(rating.data); # see the data structure
# rating.data is the rating data table with the following columns:
#      1. src_id: e.g., user_id or voter_id
#      2. dst_id: e.g., item_id or author_id
#      3. src_context: (optional) This is the context in which the source node gives the rating
#      4. dst_context: (optional) This is the context in which the destination node receives the rating
#      5. y: This is the rating that the source node gives the destination node
#      6. ctx_id: (optional) This is the context of this (src_id, dst_id) pair
# Note: You may set all/any of src_context, dst_context, ctx_id to NULL if there is no context info
#       The number of contexts cannot be too many; otherwise, the program will be very slow
str(x_obs);
# x_obs is the feature table for observations with the following columns
#      1. obs_id: observation ID (obs_id=n corresponds to the nth row of rating.data)
#      2. index:  feature index
#      3. value:  feature value
str(x_src);
# x_src is the feature table for source nodes with the following columns
#      1. src_id: source node ID (this correspond to the src_id column in rating.data)
#      2. index:  feature index
#      3. value:  feature value
str(x_dst);
# x_dst is the feature table for destination nodes with the following columns
#      1. dst_id: destination node ID (this correspond to the dst_id column in rating.data)
#      2. index:  feature index
#      3. value:  feature value

# (2) Create training/test split
set.seed(1);
select.train = sample(nrow(rating.data), floor(nrow(rating.data)*0.75));
select.test  = setdiff(1:nrow(rating.data), select.train);
obs.train = rating.data[select.train,];  x_obs.train = x_obs[x_obs$obs_id %in% select.train,];  x_obs.train$obs_id = match(x_obs.train$obs_id, select.train);
obs.test  = rating.data[select.test, ];  x_obs.test  = x_obs[x_obs$obs_id %in% select.test, ];  x_obs.test$obs_id  = match(x_obs.test$obs_id,  select.test);

# (3) Index training data
#     See src/R/model/multicontext_model_utils.R: indexData() for details
data.train = indexData(
		obs=obs.train, src.dst.same=TRUE, rm.self.link=TRUE,
		x_obs=x_obs.train, x_src=x_src, x_dst=x_dst,
		add.intercept=FALSE,
);
# (4) Index test data
#     See src/R/model/multicontext_model_utils.R: indexTestData() for details
data.test = indexTestData(
		data.train=data.train, obs=obs.test,
		x_obs=x_obs.test, x_src=x_src, x_dst=x_dst,
);
# (5) Setup the model(s) to be fitted
#     See src/R/model/multicontext_model_EM.R: run.multicontext(), fit.multicontext()
#     Note run.multicontext() is a wrapper to fit multiple models using fit.multicontext().
setting = data.frame(
		name          = c( "uv",  "vv"),
		nFactors      = c(    3,     3), # number of interaction factors
		has.u         = c(    T,     F), # whether to use u_i' v_j or v_i' v_j
		has.gamma     = c(    F,     F), # just set to F
		nLocalFactors = c(    0,     0), # just set to 0
		is.logistic   = c(    F,     F)  # whether to use the logistic model for binary rating
);
# (6) Run the fitting code
#     See src/R/model/multicontext_model_EM.R: run.multicontext(), fit.multicontext()
#     Note run.multicontext() is a wrapper to fit multiple models using fit.multicontext().
dyn.load("lib/c_funcs.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_genData.R");
source("src/R/model/multicontext_model_utils.R");
source("src/R/model/multicontext_model_MStep.R");
source("src/R/model/multicontext_model_EM.R");
source("src/R/model/GLMNet.R");
rnd.seed=1;
ans = run.multicontext(
		obs=data.train$obs,         # Observation table
		feature=data.train$feature, # Features
		setting=setting,    # Model setting
		nSamples=200,   # Number of samples drawn in each E-step: could be a vector of size nIter.
		nBurnIn=20,     # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
		nIter=10,       # Number of EM iterations
		test.obs=data.test$obs,         # Test data: Observations for testing (optional)
		test.feature=data.test$feature, #            Features for testing     (optional)
		reg.algo=GLMNet,
		IDs=data.test$IDs,
		rnd.seed.init=rnd.seed, rnd.seed.fit=rnd.seed,
		out.level=1,         # out.level=1: Save the factor & parameter values to out.dir/model.last and out.dir/model.minTestLoss
		out.dir="/tmp/test", # out.level=2: Save the factor & parameter values of each iteration i to out.dir/model.i
		out.overwrite=TRUE,     # whether to overwrite the output directory if it exists
		debug=0,      # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
		verbose=1,    # Set to 0 to disable console output; Set to 100 to print everything to the console
		verbose.M=2
);
ans$summary[,c("name", "nFactors", "has.u", "has.gamma", "nLocalFactors", "is.logistic", "best.test.loss", "last.test.loss")];

