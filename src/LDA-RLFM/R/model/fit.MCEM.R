### Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### Fit a LDA-based regression latent factor model using
### Monte-Carlo EM (through Gibbs sampling)
###
### Initialization of parameters should be done before calling this method.
###
### See README.txt for naming convention and other details
###
###  INPUT: Training Dataset
###         factor  = list(alpha, beta, gamma, u, v, s, corpus_topic); # Initial factor values
###         obs     = data.frame(y, user, item);
###         feature = list(x_dyad, x_user, x_item);
###         param   = list(b, g0, d0, c0, G, D, H, var_y, var_alpha, var_beta, var_gamma, var_u, var_v, var_s, lambda, eta);
###         corpus  = data.frame(item, term, weight); # weight can be NULL
###         try     = list(lambda, eta);
###         
###         IMPORTANT NOTE: In the training dataset, every user and item has to have observed
###         rating, feature and terms in corpus.  In the test data, it is fine to have new users
###         and items.
###
###         Test Dataset (OPTIONAL: used to compute test-set metrics over iterations)
###         Test dataset may include new items and new users
###         test.obs     = data.frame(y, user, item);
###         test.feature = list(x_dyad, x_user, x_item);
###         test.corpus  = data.frame(item, term, weight);
###
###   OPTIONS:
###     (drawTopicSample=FALSE, useTopicProb=FALSE): Disable the LDA part and use the input factor$z_avg to draw samples for factor$s.
###     (drawTopicSample=TRUE,  useTopicProb=FALSE):
###         z_avg[j,k] =  avg_w I{term w in item j belongs to topic k} from samples
###         corpus_topic[w,k] = I{term w in corpus$item[w] belongs to topic k} from samples
###     (drawTopicSample=TRUE,  useTopicProb=TRUE):
###         z_avg[j,k] =  avg_w Pr(term w in item j belongs to topic k)
###         corpus_topic[w,k] = I{term w in corpus$item[w] belongs to topic k} from samples
###     (drawTopicSample=FALSE, useTopicProb=TRUE):
###         z_avg[j,k] =  avg_w Pr(term w in item j belongs to topic k)
###         corpus_topic[w,k] = Pr(term w in corpus$item[w] belongs to topic k) (no actual multinomial sampling at all)
###
###     lambda.exponentOnNumRatings = k (default 0)
###         lambda_j = lambda * (#ratingsForItem_j + 1)^k
###
###   NOTE: (1) If you want an INTERCEPT in regression, please add the intercept column in
###             feature${x_dyad,x_user,x_item} by yourself.
###         (2) To DISABLE fitting a FACTOR, set param$var_FACTOR == NULL. Both sampling for the
###             FACTOR and regression for the FACTOR will be disabled.
###         (3) To DISABLE the RLFM part, set factor$u, factor$v, param$G, param$D to NULL.
###         (4) To DISABLE the LDA  part, set factor$s, factor$corpus_topic, param$H, and corpus to NULL.
###
### OUTPUT:
###   output$est.last:    Parameter estimates at the end of the last iteration
###   output$est.minRMSE: Parameter estimates with the lowest RMSE computed using the test dataset
###                       if is.null(test.obs), then this output will NOT be available.
###
fit.MCEM <- function(
    nIter,        # Number of EM iterations
    nSamples,     # Number of samples drawn in each E-step: could be a vector of size nIter.
    nBurnIn,      # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
    factor,       # Initial factor values
    obs,          # Dyadic observations
    feature,      # Feature values
    param,        # Initial parameter values
    corpus=NULL,  # The text corpus
    try=NULL,     # Values of lambda and eta that you want to try
    test.obs=NULL,     # Test data: Dyadic observations for testing
    test.feature=NULL, #            Feature values for testing
    test.corpus=NULL,  #            Text corpus for testing
    drawTopicSample=TRUE, # Set to TRUE, if you want to draw LDA samples; otherwise, factor$z_avg will be used.
    useTopicProb=FALSE,   # Set to TRUE, if you want z_avg[j,k] = avg_w Pr(term w in item j belongs to topic k), instead of using samples.
    lambda.exponentOnNumRatings=0,  # lambda * (#ratingsForItem_j + 1)^lambda.exponentOnNumRatings (default: 0) 
	use.MC.predict=FALSE,   # Whether to user MC.predict to test
	MC.predict.nSamples=200, # Used only when use.MC.predict==TRUE
	MC.predict.nBurnIn=10,   # Used only when use.MC.predict==TRUE
	MC.predict.useTopicProb=TRUE, # Used only when use.MC.predict==TRUE
	out.level=0,  # out.level=1: Save the factor & parameter values to out.dir/est.last and out.dir/est.minRMSE
    out.dir=NULL, # out.level=2: Save the factor & parameter values of each iteration i to out.dir/est.i
    out.append=FALSE,
    debug=0,      # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
    verbose=0,    # Set to 0 to disable console output; Set to 100 to print everything to the console
    verbose.E=verbose-1,
    verbose.M=verbose-1,
	verbose.test=0,
    use.C=TRUE,   # Whether to use the C implementation (R implementation does not have full functionalities)
    lm=F,         # Whether to use lm to fit linear regression
    error.handler=stop, # You can change it to warning so that the program won't stop
    ...           # Additional parameters passing into the regression functions (e.g., bayesglm)
){
    topicOption = 0;
    if(drawTopicSample){
        if(useTopicProb) topicOption = 2
        else             topicOption = 1;
    }else{
        if(useTopicProb) topicOption = 3;
    }
    
    obs$user = as.integer(obs$user);
    obs$item = as.integer(obs$item);
    if(!is.null(corpus)){
        corpus$item = as.integer(corpus$item);
        corpus$term = as.integer(corpus$term);
    }
    
    if(length(nSamples)!=nIter){
        if(length(nSamples) != 1) stop("length(nSamples)!=nIter && length(nSamples) != 1");
        nSamples <- rep(nSamples,nIter);
    }
    if(length(nBurnIn)!=nIter){
        if(length(nBurnIn) != 1) stop("length(nBurnIn)!=nIter && length(nBurnIn) != 1");
        nBurnIn <- rep(nBurnIn,nIter);
    }
    if(is.null(factor$gamma)) factor$gamma = rep(1.0, length(factor$alpha));
    
    warning.any.not.in(c("alpha", "beta", "gamma"), names(factor), "The following components in factor are required: ", stop=TRUE);
    size = syncheck.LDA_RLFM.spec(factor=factor, obs=obs, corpus=corpus, feature=feature, param=param, 
                                  is.corpus_topic.matrix=(topicOption == 3), warning=10, print=TRUE);
    check.ObsAndCorpus(obs, corpus, size$nUsers, size$nItems, size$nTerms, error=error.handler);
    
    nObs     = size$nObs;
    nUsers   = size$nUsers;
    nItems   = size$nItems;
    nFactors = size$nFactors;
    nTopics  = size$nTopics;

    if(!is.null(test.obs)){
        if(is.null(test.feature)) stop("test.obs != NULL, but test.feature is NULL");
        if(nTopics && is.null(test.corpus)) stop("test.obs != NULL, but test.feature is NULL");
		data.train = list(obs=obs, feature=feature, corpus=corpus);
		data.test = list(obs=test.obs, feature=test.feature, corpus=test.corpus);
	}
    
    factor = deepCopy(factor);
    if(nTopics > 0){
        if(topicOption != 0){
            count = getTopicCounts(corpus, factor$corpus_topic, nItems, nTopics, size$nTerms);
            factor$z_avg = count$z_avg;
        }else{
            if(is.null(factor$z_avg)) stop("factor$z_avg should not be NULL");
            if(is.null(factor$phi)) stop("factor$phi should not be NULL");
            if(nrow(factor$z_avg) != nItems) stop("nrow(factor$z_avg) != nItems");
            if(ncol(factor$z_avg) != nTopics) stop("ncol(factor$z_avg) != nTopics");
            try = NULL;
        }
    }
    
    LL = rep(NA, nIter+1); # LL records the complete data logLikelihood of each iteration
    bestLL = logLikelihood(obs, factor, feature, param, corpus);
    LL[1] = bestLL;

    if(!is.null(test.obs)){
        RMSE = rep(NA, nIter+1); # RMSE records the RMSE of each iteration
        if(is.null(factor$phi)) factor$phi = matrix(1/size$nTerms, nrow=nTopics, ncol=size$nTerms);
        minRMSE = pred.gauss(fit=list(factor=factor, param=param), obs=test.obs, feature=test.feature, corpus=test.corpus)$rmse;
        RMSE[1] = minRMSE;
        # Factor & paramter estimates with min RMSE
        est.minRMSE = list(factor=factor, param=param);
    }
    
    if(verbose >= 1){
        cat("START fit.MCEM (initial complete data logLikelihood: ",bestLL," + constant)\n",sep="");
        if(!is.null(test.obs))
        cat("                initial RMSE: ",minRMSE,"\n",sep="");
    }
    
    ###
    ### Output to out.dir
    ###
    if(out.level > 0){
        if(is.null(out.dir)) stop("Please specify out.dir");
        if(file.exists(out.dir) && !out.append){
            cat("Output File '",out.dir,"' EXISTS. Append? [y/n] ",sep="");
            ans = readLines(n=1);
            if(ans != "y"){
                cat("Exit!!\n");
                return(NULL);
            }
        }else if(!file.exists(out.dir)){
            dir.create(out.dir, recursive=TRUE);
        }
        save(file=paste(out.dir,"/est.0",sep=""), list=c("factor", "param"));
        thisRMSE = -1;
        if(!is.null(test.obs)){
            thisRMSE = RMSE[1];
        }
        summary = data.frame(Method="MCEM", Iter=0, nSteps=nSamples[1], CDlogL=LL[1], RMSE=thisRMSE, RMSE.simple=thisRMSE, TimeUsed1=0, TimeUsed2=0);
        file = paste(out.dir,"/summary",sep="");
        if(file.exists(file))
            write.table(summary, file=file, append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
        else
            write.table(summary, file=file, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE);
    }
    
    begin.time = proc.time();
    
    for(iter in 1:nIter){
    
        if(verbose >= 1){
            cat("---------------------------------------------------------\n",
                "Iteration ",iter,"\n",
                "---------------------------------------------------------\n",
                "start E-STEP\n",sep="");
        }
        b.time = proc.time();
        
        ###
        ### E-STEP
        ###
        if(use.C){
            mc_e = MCEM_EStep.C(
                factor, obs, corpus, feature, param, try, nSamples[iter], nBurnIn=nBurnIn[iter], 
                drawTopicSample=topicOption, nRatingExponent=lambda.exponentOnNumRatings, debug=debug, verbose=verbose.E
            );
            # IMPORTANT NOTE: factor is call-by-reference in MCEM_EStep.C
            # Now, factor = mc_e$mean with the same memory address!!!
        }else{
            mc_e = MCEM_EStep.R(
                factor, obs, corpus, feature, param, try, nSamples[iter], nBurnIn=nBurnIn[iter],
                debug=debug, verbose=verbose.E
            );
        }
        factor = mc_e$mean;
        
        time.used = proc.time() - b.time;
        time.used.1 = time.used;
        if(verbose >= 1){
            b.time.LL = proc.time();
            ll = logLikelihood(obs, factor, feature, param, corpus);
            time.used.LL = proc.time() - b.time.LL;
            cat("end   E-STEP (complete data logLikelihood = ",ll," + constant,  used ",time.used[3]," sec,  LL used ",time.used.LL[3]," sec)\n",
                "start M-STEP\n",sep="");
        }
        b.time = proc.time();
        
        ###
        ### M-STEP
        ###
        param.new = MCEM_MStep(
            factor.mean=mc_e$mean, factor.sumvar=mc_e$sumvar, obs=obs, 
            corpus=corpus, feature=feature, param=param, try=try, 
            lda.objval=mc_e$objval, debug=0, verbose=verbose.M, lm=lm, ...
        );
        
        if(verbose >= 1){
            
        }
        param = param.new;
        time.used = proc.time() - b.time;
        time.used.2 = time.used;

        b.time.LL = proc.time();
        LL[iter+1] = logLikelihood(obs, factor, feature, param, corpus);
        time.used.LL = proc.time() - b.time.LL;
        
		rmse.simple = -1;
        if(!is.null(test.obs)){
            b.time.RMSE = proc.time();
			rmse.simple = pred.gauss(
					fit=list(factor=factor, param=param), obs=test.obs, 
					feature=test.feature, corpus=test.corpus
			)$rmse;
			if(use.MC.predict){
				RMSE[iter+1] = pred.MC.gauss(
						fit=list(factor=factor, param=param), data.train=data.train, data.test=data.test,
						nSamples=min(MC.predict.nSamples,nSamples[iter]), nBurnIn=MC.predict.nBurnIn, 
						useTopicProb=MC.predict.useTopicProb, verbose=verbose.test
				)$rmse;
			}else{
				RMSE[iter+1] = rmse.simple;
			}
            time.used.RMSE = proc.time() - b.time.RMSE;
        }

        if(verbose >= 1){
            cat("end   M-STEP (complete data logLikelihood = ",LL[iter+1]," + constant,  used ",time.used[3]," sec,  LL used ",time.used.LL[3]," sec)\n",sep="");
            if(!is.null(test.obs))
            cat("              RMSE = ",RMSE[iter+1],"  used ",time.used.RMSE[3]," sec (RMSE.simple=",rmse.simple,"\n",sep="");
        }

        ###
        ### Update the est.minRMSE model if the RMSE decreases
        ###
        if(!is.null(test.obs)){ if(RMSE[iter+1] < minRMSE){
            minRMSE = RMSE[iter+1];
            est.minRMSE$factor=deepCopy(factor);
            est.minRMSE$param=param;
            
            if(verbose >= 10){
                cat("RMSE decreases!\n");
            }
        }}
        
        ###
        ### Output to out.dir
        ###
        if(out.level > 0){
            file = paste(out.dir,"/est.last",sep="");
            if(out.level >= 2){
                file.prev = paste(out.dir,"/est.",(iter-1),sep="");
                if(file.exists(file.prev)) file.remove(file.prev);
                file.rename(file, file.prev);
            }
            save(file=file, list=c("factor", "param"));
            thisRMSE = -1;
            if(!is.null(test.obs)){
                if(RMSE[iter+1] == minRMSE) file.copy(file, paste(out.dir,"/est.minRMSE",sep=""), overwrite=TRUE);
                thisRMSE = RMSE[iter+1];
            }
            summary = data.frame(Method="MCEM", Iter=iter, nSteps=nSamples[iter], CDlogL=LL[iter+1], RMSE=thisRMSE, RMSE.simple=rmse.simple, TimeUsed1=time.used.1[3], TimeUsed2=time.used.2[3]);
            file = paste(out.dir,"/summary",sep="");
            write.table(summary, file=file, append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE);
        }
    }
    
    if(is.null(test.obs)){
        RMSE = NULL;  minRMSE = NULL;  est.minRMSE = NULL;
    }
    
    output = list(
        trainingCDL = LL,
        testRMSE    = RMSE,
        minRMSE     = minRMSE,
        est.minRMSE = est.minRMSE,
        est.last = list(factor=factor, param=param)
    );

    time.used = proc.time() - begin.time;
    if(verbose >= 1){
        cat("END   fit.MCEM (used ",time.used[3]," sec)\n",sep="");
        if(!is.null(test.obs))
        cat("                best RMSE: ",minRMSE,"\n",sep="");
    }
    
    return(output);
}
