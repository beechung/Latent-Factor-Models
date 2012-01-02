
###
### Create a small test data for the multi-context model
###
out.dir = "test-data/multicontext_model/simulated-mtx-uvw-10K"
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
		nSrcContexts=4, nDstContexts=3, nEdgeContexts=11, nFactors=2, has.gamma=FALSE, has.u=TRUE,
		nObsFeatures=2, nSrcFeatures=3, nDstFeatures=3, nCtxFeatures=1,
		b.sd=1, g0.sd=1, d0.sd=1, h0.sd=0, G.sd=1, D.sd=1, H.sd=0, q.sd=1, r.sd=1,
		q.mean=5, r.mean=5,
		var_y=0.1, var_alpha=0.5, var_beta=0.5, var_gamma=1, var_v=1, var_u=1, var_w=1,
		var_alpha_global=0.2, var_beta_global=0.2,
		has.intercept=FALSE, sparse.matrices=TRUE, frac.zeroFeatures=0.1
);
# (2) Create training/test split
obs = d$obs;  names(obs) = c("src_id", "dst_id", "src_context", "dst_context", "ctx_id", "y");
set.seed(1);
select.train = runif(nrow(d$obs),min=0,max=1) < 0.75;
obs.train = obs[ select.train,];  x_obs.train = d$feature$x_obs[ select.train,,drop=FALSE];
obs.test  = obs[!select.train,];  x_obs.test  = d$feature$x_obs[!select.train,,drop=FALSE];
# (3) Write tables
# (3.1) Observations
write.table(obs.train, file=paste(out.dir,"/obs-train.txt",sep=""), row.names=FALSE, col.names=FALSE, sep="\t");
write.table(obs.test,  file=paste(out.dir,"/obs-test.txt", sep=""), row.names=FALSE, col.names=FALSE, sep="\t");
# (3.2) Features in dense format
x_src.d = data.frame(src_id=1:nrow(d$feature$x_src), as.matrix(d$feature$x_src));
x_dst.d = data.frame(dst_id=1:nrow(d$feature$x_dst), as.matrix(d$feature$x_dst));
x_ctx.d = data.frame(ctx_id=1:nrow(d$feature$x_ctx), as.matrix(d$feature$x_ctx));
x_obs.train.d = data.frame(as.matrix(x_obs.train));
x_obs.test.d  = data.frame(as.matrix(x_obs.test));
write.table(x_src.d, file=paste(out.dir,"/dense-feature-user.txt", sep=""), row.names=FALSE, col.names=FALSE, sep="\t");
write.table(x_dst.d, file=paste(out.dir,"/dense-feature-item.txt", sep=""), row.names=FALSE, col.names=FALSE, sep="\t");
write.table(x_ctx.d, file=paste(out.dir,"/dense-feature-ctxt.txt", sep=""), row.names=FALSE, col.names=FALSE, sep="\t");
write.table(x_obs.train.d, file=paste(out.dir,"/dense-feature-obs-train.txt", sep=""), row.names=FALSE, col.names=FALSE, sep="\t");
write.table(x_obs.test.d,  file=paste(out.dir,"/dense-feature-obs-test.txt",  sep=""), row.names=FALSE, col.names=FALSE, sep="\t");
# (3.2) Features in sparse format
x_src.s = matrix.to.index.value(as.matrix(d$feature$x_src), order.by="row"); names(x_src.s) = c("src_id", "index", "value");
x_dst.s = matrix.to.index.value(as.matrix(d$feature$x_dst), order.by="row"); names(x_dst.s) = c("dst_id", "index", "value");
x_ctx.s = matrix.to.index.value(as.matrix(d$feature$x_ctx), order.by="row"); names(x_ctx.s) = c("ctx_id", "index", "value");
x_obs.train.s = matrix.to.index.value(as.matrix(x_obs.train), order.by="row"); names(x_obs.train.s) = c("obs_id", "index", "value");
x_obs.test.s  = matrix.to.index.value(as.matrix(x_obs.test),  order.by="row"); names(x_obs.test.s)  = c("obs_id", "index", "value");
write.table(x_src.s, file=paste(out.dir,"/sparse-feature-user.txt", sep=""), row.names=FALSE, col.names=FALSE, sep="\t");
write.table(x_dst.s, file=paste(out.dir,"/sparse-feature-item.txt", sep=""), row.names=FALSE, col.names=FALSE, sep="\t");
write.table(x_ctx.s, file=paste(out.dir,"/sparse-feature-ctxt.txt", sep=""), row.names=FALSE, col.names=FALSE, sep="\t");
write.table(x_obs.train.s, file=paste(out.dir,"/sparse-feature-obs-train.txt", sep=""), row.names=FALSE, col.names=FALSE, sep="\t");
write.table(x_obs.test.s,  file=paste(out.dir,"/sparse-feature-obs-test.txt",  sep=""), row.names=FALSE, col.names=FALSE, sep="\t");

param  = d$param;
factor = d$factor;
save(param,factor,file=paste(out.dir,"/ground-truth.RData",sep=""));

