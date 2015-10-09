library(Matrix);
dyn.load("lib/c_funcs.so");
dyn.load("lib/arslogistic.so");
source("src/R/c_funcs.R");
source("src/arslogistic/Rfuncs.R");
nFactors = 5;
nObs = 10000;
nSamples = 1000;
data = sim.logistic(nFactors,nObs);
ars_ninit=10;
ars_qcent=seq(5,95,length=ars_ninit); # number of initial points and the quantiles of the initial points
ars_xl=-5;
ars_xu=5;
nFactors = as.integer(nFactors);
nObs = as.integer(nObs);
nSamples = as.integer(nSamples);
ars_ninit=as.integer(ars_ninit);
offset = rep(0,nObs);
samples = matrix(0,nrow=nSamples,ncol=nFactors);
beta_mean = rep(0, nFactors);
beta_var = rep(10, nFactors);

.Call("logistic_ars_Call",offset, beta_mean, beta_var, data$X, data$y, nFactors,
                        nObs, samples, ars_qcent, ars_ninit, ars_ninit, ars_xl, ars_xu, nSamples);



