### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen and Liang Zhang

###
### Quick Start
###

# (1) Read input data
input.dir = "test-data/multicontext_model/simulated-mtx-uvw-10K"
# (1.1) Training observations and observation features
obs.train = read.table(paste(input.dir,"/obs-train.txt",sep=""), 
		sep="\t", header=FALSE, as.is=TRUE);
names(obs.train) = c("src_id", "dst_id", "src_context", 
		"dst_context", "ctx_id", "y");
x_obs.train = read.table(paste(input.dir,"/dense-feature-obs-train.txt",
				sep=""), sep="\t", header=FALSE, as.is=TRUE);
# (1.2) Test observations and observation features
obs.test = read.table(paste(input.dir,"/obs-test.txt",sep=""), 
		sep="\t", header=FALSE, as.is=TRUE);
names(obs.test) = c("src_id", "dst_id", "src_context", 
		"dst_context", "ctx_id", "y");
x_obs.test = read.table(paste(input.dir,"/dense-feature-obs-test.txt",
				sep=""), sep="\t", header=FALSE, as.is=TRUE);
# (1.3) User/item/context features
x_src = read.table(paste(input.dir,"/dense-feature-user.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_src)[1] = "src_id";
x_dst = read.table(paste(input.dir,"/dense-feature-item.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_dst)[1] = "dst_id";
x_ctx = read.table(paste(input.dir,"/dense-feature-ctxt.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_ctx)[1] = "ctx_id";

# (2) Fit Models
source("src/R/BST.R");
# (2.1) Fit a model without features
ans = fit.bst(obs.train=obs.train, obs.test=obs.test, 
		out.dir="/tmp/bst/quick-start", model.name="uvw3", 
		nFactors=3, nIter=10);
# (2.2) Fit a model with features
ans = fit.bst(obs.train=obs.train, obs.test=obs.test, x_obs.train=x_obs.train, 
		x_obs.test=x_obs.test, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx,
		out.dir="/tmp/bst/quick-start", 
		model.name="uvw3-F", nFactors=3, nIter=10);

# (3) Check the Output
# (3.1) Check the summary of EM iterations
read.table("/tmp/bst/quick-start_uvw3-F/summary", header=TRUE);
# (3.2) Check the fitted model
load("/tmp/bst/quick-start_uvw3-F/model.last");
str(param);
str(factor);
str(data.train);

# (4) Make Predictions
pred = predict.bst(
		model.file="/tmp/bst/quick-start_uvw3-F/model.last",
		obs.test=obs.test, x_obs.test=x_obs.test,
		x_src=x_src, x_dst=x_dst, x_ctx=x_ctx);

# Fit Multiple Models in One Call
ans = fit.bst(obs.train=obs.train, obs.test=obs.test, x_obs.train=x_obs.train, 
		x_obs.test=x_obs.test, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx,
		out.dir = "/tmp/bst/quick-start", 
		model.name=c("uvw1", "uvw2"), nFactors=c(1,2), nIter=10);

# Fit the Original BST Model
obs.train$src_context = obs.train$dst_context = obs.train$ctx_id
obs.test$src_context  = obs.test$dst_context  = obs.test$ctx_id
ans = fit.bst(obs.train=obs.train, obs.test=obs.test,
		x_obs.train=x_obs.train, x_obs.test=x_obs.test, 
		x_src=x_src, x_dst=x_dst, x_ctx=x_ctx,
		out.dir="/tmp/bst/quick-start", model.name="original-bst", 
		nFactors=3, nIter=10, src.dst.same=TRUE,
		control=fit.bst.control(has.gamma=FALSE, rm.self.link=TRUE));

# Fit RLFM
obs.train$src_context = obs.train$dst_context = obs.train$ctx_id = NULL;
obs.test$src_context  = obs.test$dst_context  = obs.test$ctx_id  = NULL;
ans = fit.bst(obs.train=obs.train, obs.test=obs.test, 
		x_obs.train=x_obs.train, x_obs.test=x_obs.test, 
		x_src=x_src, x_dst=x_dst,
		out.dir="/tmp/bst/quick-start", model.name="uvw3-F", 
		nFactors=3, nIter=10);

###
### Appendix Example 1: Fit the BST model with dense features
###
library(Matrix);
dyn.load("lib/c_funcs.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_utils.R");
set.seed(0);

# (1) Read input data
input.dir = "test-data/multicontext_model/simulated-mtx-uvw-10K"
# (1.1) Training observations and observation features
obs.train = read.table(paste(input.dir,"/obs-train.txt",sep=""), 
            sep="\t", header=FALSE, as.is=TRUE);
names(obs.train) = c("src_id", "dst_id", "src_context", 
                     "dst_context", "ctx_id", "y");
x_obs.train = read.table(paste(input.dir,"/dense-feature-obs-train.txt",
              sep=""), sep="\t", header=FALSE, as.is=TRUE);
# (1.2) Test observations and observation features
obs.test = read.table(paste(input.dir,"/obs-test.txt",sep=""), 
           sep="\t", header=FALSE, as.is=TRUE);
names(obs.test) = c("src_id", "dst_id", "src_context", 
                    "dst_context", "ctx_id", "y");
x_obs.test = read.table(paste(input.dir,"/dense-feature-obs-test.txt",
             sep=""), sep="\t", header=FALSE, as.is=TRUE);
# (1.3) User/item/context features
x_src = read.table(paste(input.dir,"/dense-feature-user.txt",sep=""),
        sep="\t", header=FALSE, as.is=TRUE);
names(x_src)[1] = "src_id";
x_dst = read.table(paste(input.dir,"/dense-feature-item.txt",sep=""),
        sep="\t", header=FALSE, as.is=TRUE);
names(x_dst)[1] = "dst_id";
x_ctx = read.table(paste(input.dir,"/dense-feature-ctxt.txt",sep=""),
        sep="\t", header=FALSE, as.is=TRUE);
names(x_ctx)[1] = "ctx_id";

# (2) Index data: Put the input data into the right form
#     Convert IDs into numeric indices and 
#     Convert some data frames into matrices
# (2.1) Index training data
#     See src/R/model/multicontext_model_utils.R: indexData() for details
data.train = indexData(
		obs=obs.train, src.dst.same=FALSE, rm.self.link=FALSE,
		x_obs=x_obs.train, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx,
		add.intercept=TRUE
);
# (2.2) Index test data
#     See src/R/model/multicontext_model_utils.R: indexTestData() for details
data.test = indexTestData(
		data.train=data.train, obs=obs.test,
		x_obs=x_obs.test, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx
);

# (3) Setup the model(s) to be fitted
setting = data.frame(
		name          = c("uvw1", "uvw2"),
		nFactors      = c(     1,      2), # number of interaction factors
		has.u         = c(  TRUE,   TRUE), # whether to use u_i' v_j or v_i' v_j
		has.gamma     = c( FALSE,  FALSE), # whether to include gamma_k in the model
		nLocalFactors = c(     0,      0), # just set to 0
		is.logistic   = c( FALSE,  FALSE)  # whether to use the logistic response model
);

# (4) Run the fitting code
#     See src/R/model/multicontext_model_EM.R: run.multicontext() for details
dyn.load("lib/c_funcs.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_genData.R");
source("src/R/model/multicontext_model_utils.R");
source("src/R/model/multicontext_model_MStep.R");
source("src/R/model/multicontext_model_EM.R");
set.seed(2);
out.dir = "/tmp/unit-test/simulated-mtx-uvw-10K";
ans = run.multicontext(
		data.train=data.train, # training data
		data.test=data.test,   # test data (optional)
		setting=setting, # Model setting
		nSamples=200,    # Number of samples drawn in each E-step: could be a vector of size nIter.
		nBurnIn=20,      # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
		nIter=10,        # Number of EM iterations
		approx.interaction=TRUE, # predict E[uv] as E[u]E[v].
		reg.algo=NULL,     # The regression algorithm to be used in the M-step (NULL => linear regression)
		reg.control=NULL,  # The control paramter for reg.algo
		# initialization parameters
		var_alpha=1, var_beta=1, var_gamma=1, 
		var_v=1, var_u=1, var_w=1, var_y=NULL,
		relative.to.var_y=FALSE, var_alpha_global=1, var_beta_global=1,
		# others
		out.level=1,      # out.level=1: Save the factor & parameter values to out.dir/model.last and out.dir/model.minTestLoss
		out.dir=out.dir,  # out.level=2: Save the factor & parameter values of each iteration i to out.dir/model.i
		out.overwrite=TRUE,  # whether to overwrite the output directory if it exists
		debug=0,      # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
		verbose=1,    # Set to 0 to disable console output; Set to 100 to print everything to the console
		verbose.M=2,
		ridge.lambda=1, # Add diag(lambda) to X'X in linear regression
		rnd.seed.init=0, rnd.seed.fit=1
);

# Check the output
read.table(paste(out.dir,"_uvw2/summary",sep=""), header=TRUE, sep="\t", as.is=TRUE);

# Load the model
load(paste(out.dir,"_uvw2/model.last",sep=""));
# It loads param, factor, IDs, prediction
str(param);
str(factor);

# Make prediction
pred = predict.multicontext(
		model=list(factor=factor, param=param), 
		obs=data.test$obs, feature=data.test$feature, is.logistic=FALSE
);
# Now, pred$pred.y contains the predicted rating for data.test$obs
str(pred);


###
### Appendix Example 2:
###    Fit the BST model with sparse features
###    glmnet is used to fit prior regression parameters
###
library(Matrix);
dyn.load("lib/c_funcs.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_utils.R");
set.seed(0);

# (1) Read input data
input.dir = "test-data/multicontext_model/simulated-mtx-uvw-10K"
# (1.1) Training observations and observation features
obs.train = read.table(paste(input.dir,"/obs-train.txt",sep=""), 
		sep="\t", header=FALSE, as.is=TRUE);
names(obs.train) = c("src_id", "dst_id", "src_context", 
		"dst_context", "ctx_id", "y");
x_obs.train = read.table(paste(input.dir,"/sparse-feature-obs-train.txt",
				sep=""), sep="\t", header=FALSE, as.is=TRUE);
names(x_obs.train) = c("obs_id", "index", "value");
# (1.2) Test observations and observation features
obs.test = read.table(paste(input.dir,"/obs-test.txt",sep=""), 
		sep="\t", header=FALSE, as.is=TRUE);
names(obs.test) = c("src_id", "dst_id", "src_context", 
		"dst_context", "ctx_id", "y");
x_obs.test = read.table(paste(input.dir,"/sparse-feature-obs-test.txt",
				sep=""), sep="\t", header=FALSE, as.is=TRUE);
names(x_obs.test) = c("obs_id", "index", "value");
# (1.3) User/item/context features
x_src = read.table(paste(input.dir,"/sparse-feature-user.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_src) = c("src_id", "index", "value");
x_dst = read.table(paste(input.dir,"/sparse-feature-item.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_dst) = c("dst_id", "index", "value");
x_ctx = read.table(paste(input.dir,"/sparse-feature-ctxt.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_ctx) = c("ctx_id", "index", "value");

# (2) Index data: Put the input data into the right form
#     Convert IDs into numeric indices and 
#     Convert some data frames into matrices
# (2.1) Index training data
#     See src/R/model/multicontext_model_utils.R: indexData() for details
data.train = indexData(
		obs=obs.train, src.dst.same=FALSE, rm.self.link=FALSE,
		x_obs=x_obs.train, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx,
		add.intercept=TRUE
);
# (2.2) Index test data
#     See src/R/model/multicontext_model_utils.R: indexTestData() for details
data.test = indexTestData(
		data.train=data.train, obs=obs.test,
		x_obs=x_obs.test, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx
);

# (3) Setup the model(s) to be fitted
setting = data.frame(
		name          = c("uvw1", "uvw2"),
		nFactors      = c(     1,      2), # number of interaction factors
		has.u         = c(  TRUE,   TRUE), # whether to use u_i' v_j or v_i' v_j
		has.gamma     = c( FALSE,  FALSE), # whether to include gamma_k in the model
		nLocalFactors = c(     0,      0), # just set to 0
		is.logistic   = c( FALSE,  FALSE)  # whether to use the logistic response model
);

# (4) Run the fitting code
#     See src/R/model/multicontext_model_EM.R: run.multicontext() for details
dyn.load("lib/c_funcs.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_genData.R");
source("src/R/model/multicontext_model_utils.R");
source("src/R/model/multicontext_model_MStep.R");
source("src/R/model/multicontext_model_EM.R");
source("src/R/model/GLMNet.R");
set.seed(2);
out.dir = "/tmp/tutorial-BST/example-2";
ans = run.multicontext(
		data.train=data.train, # training data
		data.test=data.test,   # test data (optional)
		setting=setting,    # Model setting
		nSamples=200,   # Number of samples drawn in each E-step: could be a vector of size nIter.
		nBurnIn=20,     # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
		nIter=10,       # Number of EM iterations
		approx.interaction=TRUE, # predict E[uv] as E[u]E[v].
		reg.algo=GLMNet,   # The regression algorithm to be used in the M-step (NULL => linear regression)
		# initialization parameters
		var_alpha=1, var_beta=1, var_gamma=1, 
		var_v=1, var_u=1, var_w=1, var_y=NULL,
		relative.to.var_y=FALSE, var_alpha_global=1, var_beta_global=1,
		# others
		out.level=1,      # out.level=1: Save the factor & parameter values to out.dir/model.last and out.dir/model.minTestLoss
		out.dir=out.dir,  # out.level=2: Save the factor & parameter values of each iteration i to out.dir/model.i
		out.overwrite=TRUE,  # whether to overwrite the output directory if it exists
		debug=0,      # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
		verbose=1,    # Set to 0 to disable console output; Set to 100 to print everything to the console
		verbose.M=2,
		ridge.lambda=1, # Add diag(lambda) to X'X in linear regression
		rnd.seed.init=0, rnd.seed.fit=1
);


###
### Appendix Example 3: Add more EM iterations to an already fitted model
###
###   Example scenario: After running Example 2 with 10 EM iterations,
###   you feel that the model has not yet converged and want to add
###   5 more EM iterations to the "uvw2" model specified in the
###   setting.
###
###   To run Example 3, you must first run example 2.
###
###   Note: In the following, Step 1 and Step 2 are exactly the same
###   as those in Example 2.
###
library(Matrix);
dyn.load("lib/c_funcs.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_utils.R");
set.seed(0);

# (1) Read input data
input.dir = "test-data/multicontext_model/simulated-mtx-uvw-10K"
# (1.1) Training observations and observation features
obs.train = read.table(paste(input.dir,"/obs-train.txt",sep=""), 
		sep="\t", header=FALSE, as.is=TRUE);
names(obs.train) = c("src_id", "dst_id", "src_context", 
		"dst_context", "ctx_id", "y");
x_obs.train = read.table(paste(input.dir,"/sparse-feature-obs-train.txt",
				sep=""), sep="\t", header=FALSE, as.is=TRUE);
names(x_obs.train) = c("obs_id", "index", "value");
# (1.2) Test observations and observation features
obs.test = read.table(paste(input.dir,"/obs-test.txt",sep=""), 
		sep="\t", header=FALSE, as.is=TRUE);
names(obs.test) = c("src_id", "dst_id", "src_context", 
		"dst_context", "ctx_id", "y");
x_obs.test = read.table(paste(input.dir,"/sparse-feature-obs-test.txt",
				sep=""), sep="\t", header=FALSE, as.is=TRUE);
names(x_obs.test) = c("obs_id", "index", "value");
# (1.3) User/item/context features
x_src = read.table(paste(input.dir,"/sparse-feature-user.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_src) = c("src_id", "index", "value");
x_dst = read.table(paste(input.dir,"/sparse-feature-item.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_dst) = c("dst_id", "index", "value");
x_ctx = read.table(paste(input.dir,"/sparse-feature-ctxt.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_ctx) = c("ctx_id", "index", "value");

# (2) Index data: Put the input data into the right form
#     Convert IDs into numeric indices and 
#     Convert some data frames into matrices
data.train = indexData(
		obs=obs.train, src.dst.same=FALSE, rm.self.link=FALSE,
		x_obs=x_obs.train, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx,
		add.intercept=TRUE
);
data.test = indexTestData(
		data.train=data.train, obs=obs.test,
		x_obs=x_obs.test, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx
);

# (3) Load the "uvw2" model from Example 2.
#     If the follwoing file does not exist, run Example 2.
load("/tmp/tutorial-BST/example-2_uvw2/model.last");
model = list(factor=factor, param=param);

# (4) Run 5 additional EM iterations
dyn.load("lib/c_funcs.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_genData.R");
source("src/R/model/multicontext_model_utils.R");
source("src/R/model/multicontext_model_MStep.R");
source("src/R/model/multicontext_model_EM.R");
source("src/R/model/GLMNet.R");
out.dir = "/tmp/tutorial-BST/example-3_uvw2";
set.seed(2);
ans = fit.multicontext(
		data.train=data.train, # training data
		data.test=data.test,   # test data (optional)
		init.model=model, # Initial model = list(factor, param)
		nSamples=200, # Number of samples drawn in each E-step: could be a vector of size nIter.
		nBurnIn=20,   # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
		nIter=5,      # Number of EM iterations
		is.logistic=FALSE,
		out.level=1,     # out.level=1: Save the factor & parameter values to out.dir/model.last and out.dir/model.minTestLoss
		out.dir=out.dir, # out.level=2: Save the factor & parameter values of each iteration i to out.dir/model.i
		out.overwrite=TRUE,
		verbose=1,     # Set to 0 to disable console output; Set to 100 to print everything to the console
		verbose.M=2,
		ridge.lambda=1 # Add diag(lambda) to X'X in linear regression
);

# Check the output
read.table(paste(out.dir,"/summary",sep=""), header=TRUE, sep="\t", as.is=TRUE);

# Load the model
load(paste(out.dir,"/model.last",sep=""));
# It loads param, factor, IDs, prediction
str(param, max.level=2);
str(factor);

###
### Appendix Example 4:
###    Fit the RLFM model with sparse features
###    glmnet is used to fit prior regression parameters
###
library(Matrix);
dyn.load("lib/c_funcs.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_utils.R");
set.seed(0);

# (1) Read input data
input.dir = "test-data/multicontext_model/simulated-mtx-uvw-10K"
# (1.1) Training observations and observation features
obs.train = read.table(paste(input.dir,"/obs-train.txt",sep=""), 
		sep="\t", header=FALSE, as.is=TRUE);
names(obs.train) = c("src_id", "dst_id", "src_context", 
		"dst_context", "ctx_id", "y");
x_obs.train = read.table(paste(input.dir,"/sparse-feature-obs-train.txt",
				sep=""), sep="\t", header=FALSE, as.is=TRUE);
names(x_obs.train) = c("obs_id", "index", "value");
# (1.2) Test observations and observation features
obs.test = read.table(paste(input.dir,"/obs-test.txt",sep=""), 
		sep="\t", header=FALSE, as.is=TRUE);
names(obs.test) = c("src_id", "dst_id", "src_context", 
		"dst_context", "ctx_id", "y");
x_obs.test = read.table(paste(input.dir,"/sparse-feature-obs-test.txt",
				sep=""), sep="\t", header=FALSE, as.is=TRUE);
names(x_obs.test) = c("obs_id", "index", "value");
# (1.3) User/item/context features
x_src = read.table(paste(input.dir,"/sparse-feature-user.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_src) = c("src_id", "index", "value");
x_dst = read.table(paste(input.dir,"/sparse-feature-item.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_dst) = c("dst_id", "index", "value");
# (1.4) Ignore the context information
obs.train$src_context = obs.train$dst_context = obs.train$ctx_id = NULL;
obs.test$src_context  = obs.test$dst_context  = obs.test$ctx_id = NULL;
x_ctx = NULL;

# (2) Index data: Put the input data into the right form
#     Convert IDs into numeric indices and 
#     Convert some data frames into matrices
# (2.1) Index training data
#     See src/R/model/multicontext_model_utils.R: indexData() for details
data.train = indexData(
		obs=obs.train, src.dst.same=FALSE, rm.self.link=FALSE,
		x_obs=x_obs.train, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx,
		add.intercept=TRUE
);
# (2.2) Index test data
#     See src/R/model/multicontext_model_utils.R: indexTestData() for details
data.test = indexTestData(
		data.train=data.train, obs=obs.test,
		x_obs=x_obs.test, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx
);

# (3) Setup the model(s) to be fitted
setting = data.frame(
		name          = c( "uv1",  "uv2"),
		nFactors      = c(     1,      2), # number of interaction factors
		has.u         = c(  TRUE,   TRUE), # whether to use u_i' v_j or v_i' v_j
		has.gamma     = c( FALSE,  FALSE), # whether to include gamma_k in the model
		nLocalFactors = c(     0,      0), # just set to 0
		is.logistic   = c( FALSE,  FALSE)  # whether to use the logistic response model
);

# (4) Run the fitting code
#     See src/R/model/multicontext_model_EM.R: run.multicontext() for details
dyn.load("lib/c_funcs.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_genData.R");
source("src/R/model/multicontext_model_utils.R");
source("src/R/model/multicontext_model_MStep.R");
source("src/R/model/multicontext_model_EM.R");
source("src/R/model/GLMNet.R");
set.seed(2);
out.dir = "/tmp/tutorial-BST/example-4";
ans = run.multicontext(
		data.train=data.train, # training data
		data.test=data.test,   # test data (optional)
		setting=setting,    # Model setting
		nSamples=200,   # Number of samples drawn in each E-step: could be a vector of size nIter.
		nBurnIn=20,     # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
		nIter=10,       # Number of EM iterations
		approx.interaction=TRUE, # predict E[uv] as E[u]E[v].
		reg.algo=GLMNet,   # The regression algorithm to be used in the M-step (NULL => linear regression)
		# initialization parameters
		var_alpha=1, var_beta=1, var_gamma=1, 
		var_v=1, var_u=1, var_w=1, var_y=NULL,
		relative.to.var_y=FALSE, var_alpha_global=1, var_beta_global=1,
		# others
		out.level=1,      # out.level=1: Save the factor & parameter values to out.dir/model.last and out.dir/model.minTestLoss
		out.dir=out.dir,  # out.level=2: Save the factor & parameter values of each iteration i to out.dir/model.i
		out.overwrite=TRUE,  # whether to overwrite the output directory if it exists
		debug=0,      # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
		verbose=1,    # Set to 0 to disable console output; Set to 100 to print everything to the console
		verbose.M=2,
		ridge.lambda=1, # Add diag(lambda) to X'X in linear regression
		rnd.seed.init=0, rnd.seed.fit=1
);


###
### Appendix Example 5:
###    Fit the BST model with dense features
###    Do not give the test data to the fitting procedure,
###    and then later load the test data, index it (in the correct way)
###    and predict the response in the test data.
###
library(Matrix);
dyn.load("lib/c_funcs.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_utils.R");
set.seed(0);

# (1) Read only the training data (NOT the test data)
input.dir = "test-data/multicontext_model/simulated-mtx-uvw-10K"
# (1.1) Training observations and observation features
obs.train = read.table(paste(input.dir,"/obs-train.txt",sep=""), 
		sep="\t", header=FALSE, as.is=TRUE);
names(obs.train) = c("src_id", "dst_id", "src_context", 
		"dst_context", "ctx_id", "y");
x_obs.train = read.table(paste(input.dir,"/dense-feature-obs-train.txt",
				sep=""), sep="\t", header=FALSE, as.is=TRUE);
# (1.2) User/item/context features
x_src = read.table(paste(input.dir,"/dense-feature-user.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_src)[1] = "src_id";
x_dst = read.table(paste(input.dir,"/dense-feature-item.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_dst)[1] = "dst_id";
x_ctx = read.table(paste(input.dir,"/dense-feature-ctxt.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_ctx)[1] = "ctx_id";

# (2) Index the training data: Put the input data into the right form
#     Convert IDs into numeric indices and 
#     Convert some data frames into matrices
data.train = indexData(
		obs=obs.train, src.dst.same=FALSE, rm.self.link=FALSE,
		x_obs=x_obs.train, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx,
		add.intercept=TRUE
);

# (3) Setup the model(s) to be fitted
setting = data.frame(
		name          = c("uvw1", "uvw2"),
		nFactors      = c(     1,      2), # number of interaction factors
		has.u         = c(  TRUE,   TRUE), # whether to use u_i' v_j or v_i' v_j
		has.gamma     = c( FALSE,  FALSE), # whether to include gamma_k in the model
		nLocalFactors = c(     0,      0), # just set to 0
		is.logistic   = c( FALSE,  FALSE)  # whether to use the logistic response model
);

# (4) Run the fitting code without supplying the test data
#     See src/R/model/multicontext_model_EM.R: run.multicontext() for details
dyn.load("lib/c_funcs.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_genData.R");
source("src/R/model/multicontext_model_utils.R");
source("src/R/model/multicontext_model_MStep.R");
source("src/R/model/multicontext_model_EM.R");
set.seed(2);
out.dir = "/tmp/tutorial-BST/example-5";
ans = run.multicontext(
		data.train=data.train, # training data
		setting=setting,    # Model setting
		nSamples=200,   # Number of samples drawn in each E-step: could be a vector of size nIter.
		nBurnIn=20,     # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
		nIter=10,       # Number of EM iterations
		approx.interaction=TRUE, # predict E[uv] as E[u]E[v].
		reg.algo=NULL,     # The regression algorithm to be used in the M-step (NULL => linear regression)
		reg.control=NULL,  # The control paramter for reg.algo
		# initialization parameters
		var_alpha=1, var_beta=1, var_gamma=1, 
		var_v=1, var_u=1, var_w=1, var_y=NULL,
		relative.to.var_y=FALSE, var_alpha_global=1, var_beta_global=1,
		# others
		out.level=1,      # out.level=1: Save the factor & parameter values to out.dir/model.last and out.dir/model.minTestLoss
		out.dir=out.dir,  # out.level=2: Save the factor & parameter values of each iteration i to out.dir/model.i
		out.overwrite=TRUE,  # whether to overwrite the output directory if it exists
		debug=0,      # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
		verbose=1,    # Set to 0 to disable console output; Set to 100 to print everything to the console
		verbose.M=2,
		ridge.lambda=1, # Add diag(lambda) to X'X in linear regression
		rnd.seed.init=0, rnd.seed.fit=1
);

# Quit from R and enter R again.

out.dir = "/tmp/tutorial-BST/example-5";
# Check the output
read.table(paste(out.dir,"_uvw2/summary",sep=""), header=TRUE, sep="\t", as.is=TRUE);

# Load the model
load(paste(out.dir,"_uvw2/model.last",sep=""));
# It loads param, factor, data.train
str(param);
str(factor);
str(data.train); # this does not include the actual data!!

# Read the test data
input.dir = "test-data/multicontext_model/simulated-mtx-uvw-10K"
obs.test = read.table(paste(input.dir,"/obs-test.txt",sep=""), 
		sep="\t", header=FALSE, as.is=TRUE);
names(obs.test) = c("src_id", "dst_id", "src_context", 
		"dst_context", "ctx_id", "y");
x_obs.test = read.table(paste(input.dir,"/dense-feature-obs-test.txt",
				sep=""), sep="\t", header=FALSE, as.is=TRUE);
x_src = read.table(paste(input.dir,"/dense-feature-user.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_src)[1] = "src_id";
x_dst = read.table(paste(input.dir,"/dense-feature-item.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_dst)[1] = "dst_id";
x_ctx = read.table(paste(input.dir,"/dense-feature-ctxt.txt",sep=""),
		sep="\t", header=FALSE, as.is=TRUE);
names(x_ctx)[1] = "ctx_id";

# Index the test data
dyn.load("lib/c_funcs.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_utils.R");
data.test = indexTestData(
		data.train=data.train, obs=obs.test,
		x_obs=x_obs.test, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx
);

# Make prediction
pred = predict.multicontext(
		model=list(factor=factor, param=param), 
		obs=data.test$obs, feature=data.test$feature, is.logistic=FALSE
);
# Now, pred$pred.y contains the predicted rating for data.test$obs
str(pred);
