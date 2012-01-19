### Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### (1) Take a look at Notation.txt to understand the naming convention and
###     the model specification.
###
### (2) Compile the C code:
###     R CMD SHLIB C/util.c C/MCEM_EStep.c C/MCEM_EStep2.c C/MC_predict.c C/MC_predict2.c -o C/c_funcs.so
###

###
### (A) Generate sample data
###
source("R/model/genSimpleData.R");
source("R/utils.R");
dyn.load("c_funcs.so");
source("R/c_funcs.R");

# Generate ground truth
set.seed(0);
data.truth = genNormalData(
    nUsers=421, nItems=219, nObs=30011, nTerms=1009, corpusSize=40011,
    b=c(1, -1), var_y=0.02,
    g0=c(-2, -1, 2), d0=c(1, 2, -3, -1), c0=c(.5, .3, -.9), var_alpha=.05, var_beta=.03, var_gamma=.02,
    G=cbind(c(.9, .3, -.8), c(.4, -.7, .5)), 
    D=cbind(c(-.8, .4, .6, -.3), c(1, .3, -2, .4)),  
    H=cbind(c(2, .3, -3), c(-2, 1, 1.5), c(-.3, .1, .5)),  
    var_u=.08, var_v=.03, var_s=.09,
    eta=1, lambda=2
);

# Check whether the dataset (data.truth) has consistent dimensions and data types
# and return the dimensions (e.g., nUsers, nItems, nObs, etc).
# No error or warning means correct!! (warnings do not mean errors; I just want
# you to make sure you know what you are doing.)
# If you read real data, you should initialize your factor and param 
size = syncheck.LDA_RLFM.spec(data.truth$factor, data.truth$obs, data.truth$corpus, data.truth$feature, data.truth$param, warning=100);

# Generate training-test split
set.seed(1);
data = split.forTesting(
    data.truth,         # Input data
    nUsers.testOnly=20, # Randomly pick 20 users to be in the test set only (to test cold-start)
    nItems.testOnly=20, # Randomly pick 20 items to be in the test set only (to test cold-start)
    nObs.test=10000     # Randomly pick 10000 observations to be in the test set
);                      #      The rest 20011 observations will be the training data
# Now, data$train is the training data
#      data$test  is the test data
# Add some noise to the training data; otherwise, the initial factor and param will be the
# ground truth (because I start with the ground truth).
data.train = add.noise(data$train, factor.sd=1, factor.prob=0.3, param.coef.sd=1, param.var.sd=0.01);

###
### (B) Fit a model
###
dyn.load("c_funcs.so");
source("R/c_funcs.R");
source("R/utils.R");
source("R/model/MCEM_MStep.R");
source("R/model/fit.MCEM.R");

set.seed(2);
ans = fit.MCEM(
    nIter=5,         # Number of EM iterations
    nSamples=100,    # Number of samples drawn in each E-step: could be a vector of size nIter.
    nBurnIn=10,      # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
    factor=data.train$factor,   # Initial factor values
    obs=data.train$obs,         # Observed rating
    feature=data.train$feature, # Feature values
    param=data.train$param,     # Initial parameter values
    corpus=data.train$corpus,   # The text corpus
    try=list(lambda=c(0.5,1,2,4,8), eta=c(0.5, 1, 2, 4)), # Values of lambda and eta that you want to try
    out.level=1,                # out.level=1: Save the factor & parameter values to out.dir/est.highestCDL and out.dir/est.last
    out.dir="/tmp/test/lda-rlfm", # out.level=2: Save the factor & parameter values of each iteration i to out.dir/est.i
    out.append=FALSE,
    debug=1,     # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
    verbose=1,   # Set to 0 to disable console output; Set to 100 to print everything to the console
    verbose.M=1, # Verbose setting for the M-step
    use.C=TRUE,  # Whether to use the C implementation (R implementation does not have full functionalities)
    lm=T         # Whether to use lm to fit linear regression (otherwise bayesglm will be used, which will be slow)
);
# Now, ans$est.last is the fitted model at the last iteration
#      ans$est.last$factor: Monte-Carlo mean of the factor at the last E-step
#      ans$est.last$param:  The fitted parameter values at the last M-step
#      ans$trainingCDL:     Complete data log likelihood at each iteration
# Take a look at /tmp/test/lda-rlfm/  (see the out.dir input argument to fit.MCEM)
#      load("/tmp/test/lda-rlfm/est.last"); # You can do this even when the program is running
# Now, factor and param are the estimates from the last finished iteration.


###
### (C) Test the model
###
source("R/model/prediction.R");
eval = pred.gauss(
    fit=ans$est.last,          # Fitted model
    obs=data$test$obs,         # Observation in the test set
    feature=data$test$feature, # Feature values in the test set
    corpus=data$test$corpus,   # The text corpus in the test set
    factorOnly=F, featureOnly=F
);
str(eval);


###
### (D) Find a model having the min test-set RMSE
###
dyn.load("c_funcs.so");
source("R/c_funcs.R");
source("R/utils.R");
source("R/model/MCEM_MStep.R");
source("R/model/fit.MCEM.R");
source("R/model/prediction.R");

set.seed(2);
ans = fit.MCEM(
    nIter=5,         # Number of EM iterations
    nSamples=100,    # Number of samples drawn in each E-step: could be a vector of size nIter.
    nBurnIn=10,      # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
    factor=data.train$factor,   # Initial factor values
    obs=data.train$obs,         # Observed rating
    feature=data.train$feature, # Feature values
    param=data.train$param,     # Initial parameter values
    corpus=data.train$corpus,   # The text corpus
    try=list(lambda=c(0.5,1,2,4,8), eta=c(0.5, 1, 2, 4)), # Values of lambda and eta that you want to try
    test.obs=data$test$obs,         # TEST SET: Dyadic observations for testing
    test.feature=data$test$feature, #           Feature values for testing
    test.corpus=data$test$corpus,   #           Text corpus for testing
    out.level=1,                  # out.level=1: Save the factor & parameter values to out.dir/est.highestCDL and out.dir/est.last
    out.dir="/tmp/test/lda-rlfm", # out.level=2: Save the factor & parameter values of each iteration i to out.dir/est.i
    out.append=FALSE,
    debug=1,     # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
    verbose=1,   # Set to 0 to disable console output; Set to 100 to print everything to the console
    verbose.E=1, # Verbose setting for the E-step
    verbose.M=1, # Verbose setting for the M-step
    use.C=TRUE,  # Whether to use the C implementation (R implementation does not have full functionalities)
    lm=T         # Whether to use lm to fit linear regression (otherwise bayesglm will be used, which will be slow)
);
# Now, ans$est.minRMSE is the fitted model at the last iteration with min RMSE
#      ans$est.minRMSE$factor: Monte-Carlo mean of the factor at the E-step
#      ans$est.minRMSE$param:  The fitted parameter values at the M-step
#      ans$testRMSE:           Test set RMSE at each iteration
#      ans$est.last$factor: Monte-Carlo mean of the factor at the last E-step
#      ans$est.last$param:  The fitted parameter values at the last M-step
#      ans$trainingCDL:     Complete data log likelihood at each iteration
# Take a look at /tmp/test/lda-rlfm/  (see the out.dir input argument to fit.MCEM)
#      load("/tmp/test/lda-rlfm/est.minRMSE"); # You can do this even when the program is running
# Now, factor and param are the estimates with min RMSE in the finished iterations.


###
### (E) Find a REDUCED model with the min test-set RMSE
###
dyn.load("c_funcs.so");
source("R/c_funcs.R");
source("R/utils.R");
source("R/model/MCEM_MStep.R");
source("R/model/fit.MCEM.R");
source("R/model/prediction.R");

factor1 = data.train$factor;
param1  = data.train$param;
# Disable gamma (i.e., set gamma = 1)
factor1$gamma = NULL; param1$c0 = NULL; param1$var_gamma = NULL;
# Disable factor
factor1$u = NULL; factor1$v = NULL; param1$G = NULL; param1$D = NULL; param1$var_u = NULL; param1$var_v = NULL;
# Disable topic
factor1$s = NULL; factor1$corpus_topic = NULL; param1$H = NULL; param1$var_s = NULL; param1$eta = NULL; param1$lambda = NULL;  factor1$z_avg = NULL;

set.seed(2);
ans = fit.MCEM(
    nIter=5,         # Number of EM iterations
    nSamples=100,    # Number of samples drawn in each E-step: could be a vector of size nIter.
    nBurnIn=10,      # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
    factor=factor1,   # Initial factor values
    obs=data.train$obs,         # Observed rating
    feature=data.train$feature, # Feature values
    param=param1,               # Initial parameter values
    corpus=data.train$corpus,   # The text corpus
    try=list(lambda=c(0.5,1,2,4,8), eta=c(0.5, 1, 2, 4)), # Values of lambda and eta that you want to try
    test.obs=data$test$obs,         # TEST SET: Dyadic observations for testing
    test.feature=data$test$feature, #           Feature values for testing
    test.corpus=data$test$corpus,   #           Text corpus for testing
    out.level=1,                  # out.level=1: Save the factor & parameter values to out.dir/est.highestCDL and out.dir/est.last
    out.dir="/tmp/test/lda-rlfm", # out.level=2: Save the factor & parameter values of each iteration i to out.dir/est.i
    out.append=FALSE,
    debug=1,     # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
    verbose=1,   # Set to 0 to disable console output; Set to 100 to print everything to the console
    verbose.M=1, # Verbose setting for the M-step
    use.C=TRUE,  # Whether to use the C implementation (R implementation does not have full functionalities)
    lm=T         # Whether to use lm to fit linear regression (otherwise bayesglm will be used, which will be slow)
);
# You will see some warning message because you did not specify some parameters,
# but that is find.


###
### (F) Run plain LDA and fix (without sampling) the LDA topics
###
dyn.load("c_funcs.so");
source("R/c_funcs.R");
source("R/utils.R");
source("R/model/MCEM_MStep.R");
source("R/model/fit.MCEM.R");
source("R/model/Gibbs_LDA.R");
source("R/model/prediction.R");

# Do plain LDA
lda.out = Gibbs_LDA(
    corpus=data.train$corpus, param=data.train$param, 
    try=list(lambda=c(0.1,0.5,1,2,4,8), eta=c(0.1,0.5,1,2,4,8)), 
    nTopics=3, nIter=5, nSamples=100, nBurnIn=10, verbose=1
);

# Setup the topic distribution for items
data.train$factor$corpus_topic = lda.out$corpus_topic;
data.train$factor$z_avg = lda.out$mean$z_avg;
data.train$factor$phi = lda.out$mean$phi;
data.train$param$eta = lda.out$param$eta;
data.train$param$lambda = lda.out$param$lambda;

drawTopicSample=FALSE;

set.seed(2);
ans = fit.MCEM(
    nIter=5,         # Number of EM iterations
    nSamples=100,    # Number of samples drawn in each E-step: could be a vector of size nIter.
    nBurnIn=10,      # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
    factor=data.train$factor,   # Initial factor values
    obs=data.train$obs,         # Observed rating
    feature=data.train$feature, # Feature values
    param=data.train$param,     # Initial parameter values
    corpus=data.train$corpus,   # The text corpus
    drawTopicSample=drawTopicSample,
    test.obs=data$test$obs,         # TEST SET: Dyadic observations for testing
    test.feature=data$test$feature, #           Feature values for testing
    test.corpus=data$test$corpus,   #           Text corpus for testing
    out.level=1,                  # out.level=1: Save the factor & parameter values to out.dir/est.highestCDL and out.dir/est.last
    out.dir="/tmp/test/lda-rlfm", # out.level=2: Save the factor & parameter values of each iteration i to out.dir/est.i
    out.append=TRUE,
    debug=1,     # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
    verbose=1,   # Set to 0 to disable console output; Set to 100 to print everything to the console
    verbose.E=0, # Verbose setting for the E-step
    verbose.M=1, # Verbose setting for the M-step
    use.C=TRUE,  # Whether to use the C implementation (R implementation does not have full functionalities)
    lm=T         # Whether to use lm to fit linear regression (otherwise bayesglm will be used, which will be slow)
);


###
### (G) Use multinomial probabilities (instead of samples) for the topics by
###     setting: drawTopicSample=FALSE, useTopicProb=TRUE in fit.MCEM(...)
###
###     z_avg[j,k] =  avg_w Pr(term w in item j belongs to topic k)
###     corpus_topic[w,k] = Pr(term w in corpus$item[w] belongs to topic k), (No actual multinomial sampling)
###              instead of  I{term w in corpus$item[w] belongs to topic k} from samples
###
###     IMPORTANT NOTE: When you set drawTopicSample=FALSE and useTopicProb=TRUE,
###         factor$corpus_topic has to be a matrix of corpusSize x nTopics (instead of a vector of length corpusSize)    
###
###     Also, to use a per-item lambda, set lambda.exponentOnNumRatings = k (default 0) in fit.MCEM(...):
###         lambda_j = lambda * (#ratingsForItem_j + 1)^k
###
library(arm);
dyn.load("c_funcs.so");
source("R/c_funcs.R");
source("R/utils.R");
source("R/model/MCEM_MStep.R");
source("R/model/fit.MCEM.R");
source("R/model/Gibbs_LDA.R");
source("R/model/prediction.R");

nTopics = 3;

# Do plain LDA
lda.out = Gibbs_LDA(
    corpus=data.train$corpus, param=data.train$param, 
    try=list(lambda=c(0.001,0.005,0.01,0.05,0.1,0.5,1,2), eta=c(0.001,0.005,0.01,0.05,0.1,0.5,1,2)), 
    nTopics=nTopics, nIter=5, nSamples=100, nBurnIn=10, verbose=1
);

# Setup the topic distribution for items using unsupervised LDA as initialization
# 
data.train$factor$corpus_topic = matrix(0,nrow=nrow(data.train$corpus),ncol=nTopics);
data.train$factor$corpus_topic[cbind(1:nrow(data.train$factor$corpus_topic), lda.out$corpus_topic)] = 1;
data.train$param$eta = lda.out$param$eta;
data.train$param$lambda = lda.out$param$lambda;

drawTopicSample=FALSE;
useTopicProb=TRUE;

set.seed(2);
ans = fit.MCEM(
    nIter=5,          # Number of EM iterations
    nSamples=100,    # Number of samples drawn in each E-step: could be a vector of size nIter.
    nBurnIn=10,      # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
    factor=data.train$factor,   # Initial factor values
    obs=data.train$obs,         # Observed rating
    feature=data.train$feature, # Feature values
    param=data.train$param,     # Initial parameter values
    corpus=data.train$corpus,   # The text corpus
    try=list(lambda=c(0.001,0.005,0.01,0.05,0.1,0.5,1,2), eta=c(0.001,0.005,0.01,0.05,0.1,0.5,1,2)), # Values of lambda and eta that you want to try
    test.obs=data$test$obs,         # TEST SET: Dyadic observations for testing
    test.feature=data$test$feature, #           Feature values for testing
    test.corpus=data$test$corpus,   #           Text corpus for testing
    drawTopicSample=drawTopicSample,
    useTopicProb=useTopicProb,
    lambda.exponentOnNumRatings=1,  # lambda * (#ratingsForItem_j + 1)^lambda.exponentOnNumRatings (default: 0) 
    out.level=1,     # out.level=1: Save the factor & parameter values to out.dir/est.highestCDL and out.dir/est.last
    out.dir="/tmp/test/lda-rlfm", # out.level=2: Save the factor & parameter values of each iteration i to out.dir/est.i
    out.append=TRUE,
    debug=0,     # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
    verbose=1,   # Set to 0 to disable console output; Set to 100 to print everything to the console
    verbose.M=1, # Verbose setting for the M-step
    verbose.E=1, # Verbose setting for the E-step
    use.C=TRUE,  # Whether to use the C implementation (R implementation does not have full functionalities)
    lm=FALSE     # Whether to use lm to fit linear regression (set to FALSE to use bayesglm will be used, which will be slow)
);

