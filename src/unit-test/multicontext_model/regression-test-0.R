### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Liang Zhang

source("src/R/BST.R");

# (1) Read input data
input.dir = "test-data/multicontext_model/simulated-mtx-uvw-10K"
# (1.1) Training observations and observation features
obs.train = read.table(paste(input.dir,"/obs-train.txt",sep=""), sep="\t", header=FALSE, as.is=TRUE);
names(obs.train) = c("src_id", "dst_id", "src_context", "dst_context", "ctx_id", "y");
x_obs.train = read.table(paste(input.dir,"/dense-feature-obs-train.txt",sep=""), sep="\t", header=FALSE, as.is=TRUE);
# (1.2) Test observations and observation features
obs.test = read.table(paste(input.dir,"/obs-test.txt",sep=""), sep="\t", header=FALSE, as.is=TRUE);
names(obs.test) = c("src_id", "dst_id", "src_context", "dst_context", "ctx_id", "y");
x_obs.test = read.table(paste(input.dir,"/dense-feature-obs-test.txt",sep=""), sep="\t", header=FALSE, as.is=TRUE);
# (1.3) User/item/context features
x_src = read.table(paste(input.dir,"/dense-feature-user.txt",sep=""), sep="\t", header=FALSE, as.is=TRUE);
names(x_src)[1] = "src_id";
x_dst = read.table(paste(input.dir,"/dense-feature-item.txt",sep=""), sep="\t", header=FALSE, as.is=TRUE);
names(x_dst)[1] = "dst_id";
x_ctx = read.table(paste(input.dir,"/dense-feature-ctxt.txt",sep=""), sep="\t", header=FALSE, as.is=TRUE);
names(x_ctx)[1] = "ctx_id";

# (2) Call BST
ans = fit.bst(obs.train=obs.train, obs.test=obs.test, out.dir = "/tmp/unit-test/simulated-mtx-uvw-10K", model.name=c("uvw1", "uvw2"), nFactors=c(1,2), nIter=10);
#ans = fit.bst(obs.train=obs.train, obs.test=obs.test, x_obs.train=x_obs.train, x_obs.test=x_obs.test, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx,
#        out.dir = "/tmp/unit-test/simulated-mtx-uvw-10K", model.name=c("uvw1", "uvw2"), nFactors=c(1,2), nIter=10);


# (3) Compare to the reference run
warnings()
out.dir = "/tmp/unit-test/simulated-mtx-uvw-10K"
setting = data.frame(
                name          = c("uvw1", "uvw2"),
                nFactors      = c(     1,      2), # number of interaction factors
                has.u         = c(  TRUE,   TRUE), # whether to use u_i' v_j or v_i' v_j
                has.gamma     = c( FALSE,  FALSE), # just set to F
                nLocalFactors = c(     0,      0), # just set to 0
                is.logistic   = c( FALSE,  FALSE)  # whether to use the logistic model for binary rating
);
ok = TRUE;
for(i in 1:nrow(setting)){
        name = setting[i,"name"];
        smry.1 = read.table(paste(input.dir,"/mtx-",name,".summary.txt",sep=""), header=TRUE, as.is=TRUE);
        smry.2 = read.table(paste(out.dir,"_",name,"/summary",sep=""), header=TRUE, as.is=TRUE);
        for(metric in c("TestLoss", "LossInTrain")){
                diff = abs(smry.1[,metric] - smry.2[,metric]);
                if(any(diff > 1e-10)){
                        ok=FALSE; cat("Observe Difference in ",metric,"\n",sep="");
                        print(data.frame(smry.1["Iter"], smry.1[metric], smry.2[metric], Diff=diff, "."=c("", "***")[(diff>0)+1]));
                }
        }
}
if(ok){
        cat("\nNo Problem Found!!\n\n",sep="")
}else{
        cat("\nSome Problems Found!!\n\n",sep="");
}
