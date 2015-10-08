
###
### Example: Run the fitting code with synthetic data using SPARSE feature matrix
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
		nSrcContexts=1, nDstContexts=1, nEdgeContexts=0, nFactors=3, has.gamma=FALSE, has.u=TRUE,
		nObsFeatures=13, nSrcFeatures=19, nDstFeatures=23,
		b.sd=1, g0.sd=1, d0.sd=1, G.sd=1, D.sd=1,
		var_y=0.1, var_alpha=0.5, var_beta=0.5, var_v=1, var_u=1,
		has.intercept=TRUE, binary.response=TRUE,
		sparse.matrices=TRUE, index.value.format=TRUE, frac.zeroFeatures=0.5
);
names(d$obs) = c("src_id", "dst_id", "y");
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
#      1. src_id: e.g., user_id
#      2. dst_id: e.g., item_id
#      3. y: This is the rating that the source node gives the destination node
str(x_obs);
# x_obs is the feature table for observations with the following columns
#      1. obs_id: observation ID (obs_id=n corresponds to the nth row of rating.data)
#      2. index:  feature index
#      3. value:  feature value
#   NOTE: The first column must be the intercept (i.e., a column with value 1)
str(x_src);
# x_src is the feature table for source nodes with the following columns
#      1. src_id: source node ID (this correspond to the src_id column in rating.data)
#      2. index:  feature index
#      3. value:  feature value
#   NOTE: The first column must be the intercept (i.e., a column with value 1)
str(x_dst);
# x_dst is the feature table for destination nodes with the following columns
#      1. dst_id: destination node ID (this correspond to the dst_id column in rating.data)
#      2. index:  feature index
#      3. value:  feature value
#   NOTE: The first column must be the intercept (i.e., a column with value 1)

# (2) Create training/test split
set.seed(1);
select.train = sample(nrow(rating.data), floor(nrow(rating.data)*0.75));
select.test  = setdiff(1:nrow(rating.data), select.train);
obs.train = rating.data[select.train,];  x_obs.train = x_obs[x_obs$obs_id %in% select.train,];  x_obs.train$obs_id = match(x_obs.train$obs_id, select.train);
obs.test  = rating.data[select.test, ];  x_obs.test  = x_obs[x_obs$obs_id %in% select.test, ];  x_obs.test$obs_id  = match(x_obs.test$obs_id,  select.test);

# (3) Index training data
#     See src/R/model/multicontext_model_utils.R: indexData() for details
data.train = indexData(
		obs=obs.train, src.dst.same=FALSE,
		x_obs=x_obs.train, x_src=x_src, x_dst=x_dst,
		add.intercept=FALSE,
);
# (4) Index test data
#     See src/R/model/multicontext_model_utils.R: indexTestData() for details
data.test = indexTestData(
		data.train=data.train, obs=obs.test,
		x_obs=x_obs.test, x_src=x_src, x_dst=x_dst,
);
# (5) Run the fitting code
#     See src/RLFM-ars-logistic/R/fit.MCEM.logistic.R: fit.ARS.logistic()
#     NOTE: The order of the following source(...) calls might be important
#           since we have not yet clean up the functions. Some function may
#           need to be overwritten in order for the code to work.
library(glmnet);
dyn.load("lib/c_funcs.so");
dyn.load("lib/arslogistic.so");
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/multicontext_model_EM.R");
source("src/R/model/multicontext_model_utils.R");
source("src/RLFM-ars-logistic/R/c_funcs.R");
source("src/RLFM-ars-logistic/R/util.R");
source("src/RLFM-ars-logistic/R/fit.MCEM.logistic.R");
source("src/RLFM-ars-logistic/R/regression.R");
set.seed(1); # NOTE: set.seed doesn't work because the ARS code uses its own random number generator
ans =  fit.ARS.logistic(
        nIter=10,      # Number of EM iterations
        nSamples=10, nBurnin=20,  # Number of samples and burnin drawn in each E-step: could be a vector of size nIter.
        data.train=data.train, # Training data = list(obs, feature)
        nFactors=3,   # Number of factors (i.e., number of dimensions of u)
        init.model=NULL, # Initial model = list(factor, param). Set to NULL to use the default.
        # initialization parameters
        var_alpha=1, var_beta=1, var_v=1, var_u=1,
        # others
        out.level=2,  # out.level=1: Save the parameter values out.dir/est.highestCDL and out.dir/est.last
                      # out.level=2: Save the parameter values of each iteration i to out.dir/est.i
        out.dir="/tmp/fit-ARS-logistic-example",
        out.append=TRUE,
        debug=10,      # Set to 0 to disable internal sanity checking; Set to 10 for most detailed sanity checking
        verbose=10,    # Set to 0 to disable console output; Set to 10 to print everything to the console
        use.glmnet=TRUE,
        # ARS parameters
        ars_ninit=3, ars_qcent=c(5.0,50.0,95.0), # number of initial points and the quantiles of the initial points
        ars_xl=-5, ars_xu=5, # lower bound and upper bound of ARS samples
        ars_alpha=0.5,
        center=FALSE # center the random effects at every iteration of the ARS?
);
