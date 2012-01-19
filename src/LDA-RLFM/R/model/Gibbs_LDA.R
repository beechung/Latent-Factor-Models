### Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### Gibbs_LDA
###
### IN&OUT: corpus_topic
###         INPUT:  initial topic assignment for each term
###         OUTPUT: last sample of corpus_topic
###         THE VALUES in this array WILL BE CHANGED BY THIS FUNCTION (CALL BY REFERENCE!!)
###
### INPUT:  corpus  = data.frame(item, term, weight);
###         param   = list(lambda, eta);
###         try     = list(lambda, eta);
###         
### OUTPUT: mean    = list(z_avg, phi);
###         param   = list(eta, lambda);
###         objval  = list(eta, lambda);
###
Gibbs_LDA <- function(
    corpus_topic=NULL, corpus, param=NULL, try, nTopics, nIter, nSamples, nBurnIn=1,
    debug=0, verbose=0
){
    if(is.null(corpus_topic)){
        corpus_topic = as.integer(floor(runif(nrow(corpus),1,nTopics+1-1e-20)));
    }
    if(is.null(param)){
        param = list(lambda=try$lambda[1], eta=try$eta[1]);
    }
    
    if(length(corpus_topic) != nrow(corpus)) stop("length(corpus_topic) != nrow(corpus)");
    if(max(corpus_topic) > nTopics) stop("max(corpus_topic) > nTopics");
    if(max(corpus_topic) != nTopics) warning("max(corpus_topic) != nTopics");

    if(!is.double(try$lambda)) stop("!is.double(try$lambda)");
    if(!is.double(try$eta)) stop("!is.double(try$eta)");
    if(!is.double(param$lambda)) stop("!is.double(param$lambda)");
    if(!is.double(param$eta)) stop("!is.double(param$eta)");
    
    for(i in 1:nIter){
        if(verbose > 0) cat("Iteration ",i,"\n",sep="");
        ans = Gibbs_LDA_EStep.C(corpus_topic, corpus, param, try, nTopics, nSamples, nBurnIn, 
                                debug, verbose-1);
        param$lambda = try$lambda[which.min(ans$objval$lambda)];
        param$eta    = try$eta[   which.min(ans$objval$eta)];
        if(verbose > 0){
            temp = data.frame(lambda=try$lambda, lambda.obj=ans$objval$lambda);
            print(temp);
            cat("lambda = ",param$lambda,"\n",sep="");
            temp = data.frame(eta=try$eta, eta.obj=ans$objval$eta);
            print(temp);
            cat("eta = ",param$eta,"\n",sep="");
        }
    }
    
    out = list(corpus_topic=corpus_topic, mean=ans$mean, param=param, objval=ans$objval);
    return(out);
}
