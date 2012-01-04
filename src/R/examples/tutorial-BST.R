
###
### Example 1: Fit the BST model with dense features
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

# (3) Fit models
# (3.1) Setup the model(s) to be fitted
setting = data.frame(
		name          = c("uvw1", "uvw2"),
		nFactors      = c(     1,      2), # number of interaction factors
		has.u         = c(  TRUE,   TRUE), # whether to use u_i' v_j or v_i' v_j
		has.gamma     = c( FALSE,  FALSE), # whether to include gamma_k in the model
		nLocalFactors = c(     0,      0), # just set to 0
		is.logistic   = c( FALSE,  FALSE)  # whether to use the logistic response model
);
# (3.2) Run the fitting code
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
		obs=data.train$obs,         # Observation table
		feature=data.train$feature, # Features
		setting=setting,    # Model setting
		nSamples=200,   # Number of samples drawn in each E-step: could be a vector of size nIter.
		nBurnIn=20,     # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
		nIter=10,       # Number of EM iterations
		test.obs=data.test$obs,         # Test data: Observations for testing (optional)
		test.feature=data.test$feature, #            Features for testing     (optional)
		approx.interaction=TRUE, # predict E[uv] as E[u]E[v].
		reg.algo=NULL,     # The regression algorithm to be used in the M-step (NULL => linear regression)
		reg.control=NULL,  # The control paramter for reg.algo
		# initialization parameters
		var_alpha=1, var_beta=1, var_gamma=1, 
		var_v=1, var_u=1, var_w=1, var_y=NULL,
		relative.to.var_y=FALSE, var_alpha_global=1, var_beta_global=1,
		# others
		IDs=data.test$IDs,
		out.level=1,      # out.level=1: Save the factor & parameter values to out.dir/model.last and out.dir/model.minTestLoss
		out.dir=out.dir,  # out.level=2: Save the factor & parameter values of each iteration i to out.dir/model.i
		out.overwrite=TRUE,  # whether to overwrite the output directory if it exists
		debug=0,      # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
		verbose=1,    # Set to 0 to disable console output; Set to 100 to print everything to the console
		verbose.M=2,
		ridge.lambda=1, # Add diag(lambda) to X'X in linear regression
		rnd.seed.init=0, rnd.seed.fit=1
);
