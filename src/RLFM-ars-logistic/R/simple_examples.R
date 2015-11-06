library(Matrix);
source("src/R/c_funcs.R");
source("src/R/util.R");
source("src/R/model/util.R");
source("src/R/model/multicontext_model_genData.R");
source("src/R/model/multicontext_model_utils.R");
set.seed(0);

d = genMainEffectData(
  nSrcNodes=3, nDstNodes=5, nObs=1000000, intercept=-1,
  var_y=1, var_alpha=0.5, var_beta=0.5, binary.response=TRUE
)

dyn.load("lib/c_funcs.so");
dyn.load("lib/arslogistic.so");
source("src/R/model/multicontext_model_EM.R");
source("src/R/model/multicontext_model_utils.R");
source("src/RLFM-ars-logistic/R/c_funcs.R");
source("src/RLFM-ars-logistic/R/util.R");
source("src/RLFM-ars-logistic/R/fit.MCEM.logistic.R");
source("src/RLFM-ars-logistic/R/regression.R");
set.seed(1); # NOTE: set.seed doesn't work because the ARS code uses its own random number generator
library(arm)

ans =  fit.ARS.logistic(
        nIter=3,      # Number of EM iterations
        nSamples=100, nBurnin=20,  # Number of samples and burnin drawn in each E-step: could be a vector of size nIter.
        data.train=list(obs=d$obs,feature=d$feature), # Training data = list(obs, feature)
        nFactors=3,   # Number of factors (i.e., number of dimensions of u)
        init.model=NULL, # Initial model = list(factor, param). Set to NULL to use the default.
        # initialization parameters
        var_alpha=0.1, var_beta=0.1, var_v=0.1, var_u=0.1,
        # others
        out.level=2,  # out.level=1: Save the parameter values out.dir/est.highestCDL and out.dir/est.last
                      # out.level=2: Save the parameter values of each iteration i to out.dir/est.i
        out.dir="./fit-ARS-logistic-example",
        out.append=TRUE,
        debug=10,      # Set to 0 to disable internal sanity checking; Set to 10 for most detailed sanity checking
        verbose=10,    # Set to 0 to disable console output; Set to 10 to print everything to the console
        use.glmnet=FALSE,
        # ARS parameters
        ars_ninit=3, ars_qcent=c(5.0,50.0,95.0), # number of initial points and the quantiles of the initial points
        ars_xl=-3, ars_xu=3, # lower bound and upper bound of ARS samples
        ars_alpha=0.5,
        main.effects=TRUE,
        fit.regression=FALSE,
        identifiable=FALSE
);
