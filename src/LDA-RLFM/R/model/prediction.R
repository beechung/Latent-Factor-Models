### Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### Predict for the Gaussian case using Monte-Carlo mean
###
###     See reindexData and indexTestData for how to create data.train and data.test.
###
### NOTES:
###   1. Users that are only in data.test must be indexed after all the old users
###   2. Items that are only in data.test must be indexed after all the old items
###   3. data.train$termLevels and data.test$termLevels should match
###   4. fit$factor are the factors of the users and items in data.train only.
###
pred.MC.gauss <- function(
    fit,        #  Fitted model: list(param, factor)
    data.train, # Training data: list(obs, feature, corpus, userIDs, itemIDs, termLevels)
    data.test,  #     Test data: list(obs, feature, corpus, userIDs, itemIDs, termLevels)
    nSamples, nBurnIn=1, useTopicProb=TRUE,
    test.term.weight=1, transduction=FALSE,
    use.C=TRUE, lognorm.sd=0.5,
    debug=0, verbose=0
){
    param = fit$param;
    factor.train = fit$factor;
    
    # sanity check
    size = syncheck.LDA_RLFM.spec(factor.train, data.train$obs, data.train$corpus, data.train$feature, param, is.corpus_topic.matrix=is.matrix(factor.train$corpus_topic));
    check.included(data.train$feature$x_user, data.test$feature$x_user, TRUE, "data.train$feature$x_user", "data.test$feature$x_user");
    check.included(data.train$feature$x_item, data.test$feature$x_item, TRUE, "data.train$feature$x_item", "data.test$feature$x_item");
    if(max(data.train$obs$user) > size$nUsers) stop("max(data.train$obs$user) > size$nUsers");
    if(max(data.train$obs$item) > size$nItems) stop("max(data.train$obs$item) > size$nItems");
    if(max(data.test$obs$user) > nrow(data.test$feature$x_user)) stop("max(data.test$obs$user) > nrow(data.test$feature$x_user)");
    if(max(data.test$obs$item) > nrow(data.test$feature$x_item)) stop("max(data.test$obs$item) > nrow(data.test$feature$x_item)");
    check.corpus(data.train$corpus, data.test$corpus, size$nTopics, data.train$termLevels, data.test$termLevels);
    check.included(data.train$userIDs, data.test$userIDs, FALSE, "data.train$userIDs", "data.test$userIDs");
    check.included(data.train$itemIDs, data.test$itemIDs, FALSE, "data.train$itemIDs", "data.test$itemIDs");
    
    # initialize factors
    factor = factor.train;
    nTestUsers = nrow(data.test$feature$x_user);
    nTestItems = nrow(data.test$feature$x_item);
    if(nTestUsers > size$nUsers){
        factor$alpha = c(factor.train$alpha, rep(0.0, nTestUsers-size$nUsers));
        factor$gamma = c(factor.train$gamma, rep(1.0, nTestUsers-size$nUsers));
        if(size$nFactors > 0){
            factor$u = matrix(0.0, nrow=nTestUsers, ncol=size$nFactors);
            factor$u[1:size$nUsers,] = factor.train$u;
        }
        if(size$nTopics > 0){
            factor$s = matrix(0.0, nrow=nTestUsers, ncol=size$nTopics);
            factor$s[1:size$nUsers,] = factor.train$s;
        }
    }
    corpus = data.train$corpus;
    if(nTestItems > size$nItems){
        factor$beta = c(factor.train$beta, rep(0.0, nTestItems-size$nItems));
        if(size$nFactors > 0){
            factor$v = matrix(0.0, nrow=nTestItems, ncol=size$nFactors);
            factor$v[1:size$nItems,] = factor.train$v;
        }
        if(size$nTopics > 0){
            
            corpus.new = data.test$corpus[!(data.test$corpus$item %in% data.train$corpus$item),];
            
            if(is.null(corpus.new$weight)) corpus.new$weight = test.term.weight
            else                           corpus.new$weight = corpus.new$weight * test.term.weight;
            
            corpus = rbind(data.train$corpus, corpus.new);
            if(nrow(corpus) != nrow(data.test$corpus)) stop("error");
            
            rnd_topics = as.integer(floor(runif(nrow(data.test$corpus)-size$corpusSize,1,size$nTopics+1-1e-20)));
            
            if(is.matrix(factor.train$corpus_topic)){
                factor$corpus_topic = matrix(as.integer(0), nrow=nrow(corpus), ncol=size$nTopics);
                factor$corpus_topic[1:size$corpusSize,] = factor.train$corpus_topic;
                index.new = (size$corpusSize+1):nrow(factor$corpus_topic);
                factor$corpus_topic[index.new,] = exp(rnorm(length(index.new)*size$nTopics, 0, sd=lognorm.sd));
                factor$corpus_topic[index.new,] = normalize_sumToOne2D(factor$corpus_topic[index.new,,drop=FALSE], 1);
            }else{
                factor$corpus_topic = c(factor.train$corpus_topic, rnd_topics);
            }
        }
    }
    param$var_y = mean(param$var_y);
    size = syncheck.LDA_RLFM.spec(factor, data.test$obs, corpus, data.test$feature, param, is.corpus_topic.matrix=is.matrix(factor$corpus_topic));
    
    if(use.C) pred.func = MC_predict.C
    else      pred.func = MC_predict.R;
    
    ans = pred.func(
        factor=factor, obs.train=data.train$obs, obs.test=data.test$obs, 
        corpus=corpus, feature=data.test$feature, param=param, 
        x_dyad.train=data.train$feature$x_dyad, x_dyad.test=data.test$feature$x_dyad,
        nSamples=nSamples, nBurnIn=nBurnIn, useTopicProb=useTopicProb, transduction=transduction,
        debug=debug, verbose=verbose
    );
    return(ans);
}
check.included <- function(x, y, not.null, x.name, y.name){
    if(!not.null){
        if(is.null(x) && !is.null(y)) stop(x.name," is null, but ",y.name," is not");
        if(is.null(y) && !is.null(x)) stop(y.name," is null, but ",x.name," is not");
        if(is.null(x) && is.null(y)) return(0);
    }
    if(is.vector(x)){
        if(!is.vector(y)) stop(x.name," is a vector, but ",y.name," is not");
        if(length(x) > length(y)) stop(x.name," and ",y.name," do not match! (length)");
        if(any(x != y[1:length(x)])) stop(x.name," and ",y.name," do not match! (content)");
    }else if(is.matrix(x)){
        if(!is.matrix(y)) stop(x.name," is a matrix, but ",y.name," is not");
        if(ncol(x) != ncol(y)) stop(x.name," and ",y.name," do not match! (ncol)");
        if(nrow(x) >  nrow(y)) stop(x.name," and ",y.name," do not match! (nrow)");
        if(any(x != y[1:nrow(x),])) stop(x.name," and ",y.name," do not match! (content)");
    }
}
check.corpus <- function(x, y, nTopics, levels1, levels2){
    if(nTopics > 0){
        if(is.null(x)) stop("training corpus is null");
        if(is.null(y)) stop("test corpus is null");
        if(is.null(x$weight) && !is.null(y$weight)) stop("training corpus does not match test corpus (weight)");
        if(is.null(y$weight) && !is.null(x$weight)) stop("training corpus does not match test corpus (weight)");
        if(nrow(x) > nrow(y)) stop("training corpus does not match test corpus (nrow)");
        c1 = x[order(x$item, x$term),];
        temp = y[y$item %in% x$item,];
        c2 = temp[order(temp$item, temp$term),];
        if(nrow(c1) != nrow(c2)) stop("training corpus does not match test corpus (nrow)");
        if(any(c1$item != c2$item)) stop("training corpus does not match test corpus (item)");
        if(any(c1$term != c2$term)) stop("training corpus does not match test corpus (term)");
        if(!is.null(c1$weight) && any(c1$weight != c2$weight)) stop("training corpus does not match test corpus (weight)");
        if(is.null(levels1) && !is.null(levels2)) stop("training corpus does not match test corpus (termLevels)");
        if(is.null(levels2) && !is.null(levels1)) stop("training corpus does not match test corpus (termLevels)");
        if(!is.null(levels1) && (length(levels1) != length(levels2) || any(levels1 != levels2))) stop("training corpus does not match test corpus (termLevels)");
    }
}

###
### Predict for the Gaussian case
###
###     User 1 to length(fit$factor$alpha) are old users
###     Item 1 to length(fit$factor$beta)  are old items
###
### To disable a FACTOR, set var_FACTOR = NULL
###
pred.gauss <- function(
    fit,          # Fitted model
    obs,          # Observation
    feature,      # Feature values
    corpus=NULL,  # The text corpus
    factorOnly=F, featureOnly=F,
    metrics=c(),  # an array of the following "kendall", "pearson"
    useGlobalTopicProb=FALSE, # Whether you Pr[topic] when normalizing phi
    output.factor=FALSE
){
    userIDs = unique(obs$user);
    itemIDs = unique(obs$item);
    
    nOldUsers = length(fit$factor$alpha);  nNewUsers = max(userIDs);
    nOldItems = length(fit$factor$beta);   nNewItems = max(itemIDs);
    nFactors = if(is.null(fit$factor$u)) 0 else ncol(fit$factor$u);
    nTopics  = if(is.null(fit$factor$s)) 0 else ncol(fit$factor$s);
    
    oldUserIDs <- userIDs[userIDs <= nOldUsers];
    newUserIDs <- userIDs[userIDs > nOldUsers];
    oldItemIDs <- itemIDs[itemIDs <= nOldItems];
    newItemIDs <- itemIDs[itemIDs > nOldItems];

    if(!factorOnly){
        if(is.null(feature$x_user) || nrow(feature$x_user) < nNewUsers) stop("some new users in test don't have covariates");
        if(is.null(feature$x_item) || nrow(feature$x_item) < nNewItems) stop("some new items in test don't have covariates");
    }
    if(is.null(feature$x_dyad) || nrow(feature$x_dyad) != nrow(obs)) stop("is.null(feature$x_dyad) || nrow(feature$x_dyad) != nrow(obs)");
    if(!is.null(corpus)){
        temp = sort(unique(corpus$item));
        if(length(temp) < nNewItems || any(temp[1:nNewItems] != 1:nNewItems)) stop("Some items do not have terms in the corpus.");
    }

    userIDs = c(oldUserIDs, newUserIDs);
    itemIDs = c(oldItemIDs, newItemIDs);
    
    factor = list(
        alpha = getReg1DFactor(fit$factor$alpha, feature$x_user, fit$param$g0, oldUserIDs, newUserIDs, factorOnly, featureOnly, disable=is.null(fit$param$var_alpha)),
        beta  = getReg1DFactor(fit$factor$beta,  feature$x_item, fit$param$d0, oldItemIDs, newItemIDs, factorOnly, featureOnly, disable=is.null(fit$param$var_beta)),
        gamma = getReg1DFactor(fit$factor$gamma, feature$x_user, fit$param$c0, oldUserIDs, newUserIDs, factorOnly, featureOnly, disable=is.null(fit$param$var_gamma), defaultValue=1),
        u     = getReg2DFactor(fit$factor$u, feature$x_user, fit$param$G, oldUserIDs, newUserIDs, factorOnly, featureOnly, disable=is.null(fit$param$var_u)),
        v     = getReg2DFactor(fit$factor$v, feature$x_item, fit$param$D, oldItemIDs, newItemIDs, factorOnly, featureOnly, disable=is.null(fit$param$var_v)),
        s     = getReg2DFactor(fit$factor$s, feature$x_user, fit$param$H, oldUserIDs, newUserIDs, factorOnly, featureOnly, disable=is.null(fit$param$var_s)),
        z_avg = getTopicFactor(fit$factor$z_avg, fit$factor$corpus_topic, corpus, fit$factor$phi, oldItemIDs, newItemIDs, disable=(nTopics <= 0), useGlobalTopicProb=useGlobalTopicProb)
    );
    
    if(is.null(factor$gamma)) factor$gamma = rep(1.0, nNewUsers);
    obs$user = match(obs$user, userIDs);
    obs$item = match(obs$item, itemIDs);
    
    pred.y = predict.y.from.factors(obs, factor, feature, fit$param)
    rmse   = sqrt(mean( (obs$y - pred.y)^2 ));
    mae    = mean( abs(obs$y - pred.y) );

    out = list(pred.y=pred.y, rmse=rmse, mae=mae);

    if("kendall" %in% metrics){
        out$kendall = cor(obs$y, pred.y, method="kendall");
    }else if("pearson" %in% metrics){
        out$kendall = cor(obs$y, pred.y, method="pearson");
    }
    if(output.factor) out$factor = factor;
  
    return(out);
}

###
### Input (factor, feature, param) is one of:
###       (alpha, x_user, g0), (beta, x_item, d0), (gamma, x_user, c0)
###
### output[1:length(oldIDs)] the factor for oldIDs
### output[length(oldIDs)+(1:length(newIDs))] the factor for newIDs
###
getReg1DFactor <- function(factor, feature, param, oldIDs, newIDs, factorOnly, featureOnly, disable, defaultValue=0){
    if(factorOnly && featureOnly) stop("factorOnly && featureOnly");
    if(disable) return(NULL);
    
    out = rep(NA, length(oldIDs)+length(newIDs));
    if(featureOnly){
        out[] = feature[c(oldIDs, newIDs),,drop=FALSE] %*% param;
    } else {
        if(length(oldIDs) > 0) out[1:length(oldIDs)] = factor[oldIDs];
        if(length(newIDs) > 0){
            out[length(oldIDs)+(1:length(newIDs))] <- if(factorOnly) defaultValue else feature[newIDs,,drop=FALSE] %*% param;
        }
    }
    return(out);
}

###
### Input (factor, feature, param) is one of:
###       (u,     x_user, G ), (v,    x_item  D ), (s,     x_user, H )
###
### output[1:length(oldIDs),] the factor for oldIDs
### output[length(oldIDs)+(1:length(newIDs)),] the factor for newIDs
###
getReg2DFactor <- function(factor, feature, param, oldIDs, newIDs, factorOnly, featureOnly, disable, defaultValue=0){
    if(factorOnly && featureOnly) stop("factorOnly && featureOnly");
    nFactors = 0;
    if(!is.null(factor)) nFactors = ncol(factor);
    if(!is.null(param)){
        if(nFactors != 0 && nFactors != ncol(param)) stop("error");
        nFactors = ncol(param);
    }
    
    if(disable || nFactors == 0) return(NULL);
    
    out = matrix(NA, nrow=length(oldIDs)+length(newIDs), ncol=nFactors);
    if(featureOnly){
        out[,] = feature[c(oldIDs, newIDs),,drop=FALSE] %*% param;
    } else {
        if(length(oldIDs) > 0) out[1:length(oldIDs),] = factor[oldIDs,];
        if(length(newIDs) > 0){
            out[length(oldIDs)+(1:length(newIDs)),] <- if(factorOnly) defaultValue else feature[newIDs,,drop=FALSE] %*% param;
        }
    }
    return(out);
}

###
### Get z_avg
###
### output[1:length(oldIDs),] the factor for oldIDs
### output[length(oldIDs)+(1:length(newIDs)),] the factor for newIDs
###
getTopicFactor <- function(z_avg, corpus_topic, corpus, phi, oldIDs, newIDs, disable, useGlobalTopicProb=TRUE){
    if(disable) return(NULL);
    if(ncol(z_avg) != nrow(phi)) stop("error!!");
    out = matrix(NA, nrow=length(oldIDs)+length(newIDs), ncol=ncol(z_avg));
    if(length(oldIDs) > 0) out[1:length(oldIDs),] = z_avg[oldIDs,];
    if(length(newIDs) > 0){
        
        corpus = corpus[corpus$item %in% newIDs,];
        corpus_pos = match(corpus$item, newIDs);  # corpus_pos[x] = j iff newID[j] = corpus$item[x]
        
        if(useGlobalTopicProb){
            
            if(is.matrix(corpus_topic)){
                # corpus_topic: corpusSize x nTopics
                pr_topic = sum_margin(corpus_topic, side=2);
            }else{
                # temp = corpus$weight;
                # if(is.null(temp)) temp = rep(1, length(corpus_topic));
                temp = rep(1, length(corpus_topic));  # no weights for now
                
                prob_topic = aggregate(temp, list(topic=corpus_topic), sum);
                pr_topic = prob_topic$x / sum(prob_topic$x);
                
                if(length(prob_topic$topic) != nrow(phi) || any(prob_topic$topic != 1:nrow(phi))){
                    warning("Some topics have no probability!");
                    temp = pr_topic;
                    pr_topic = rep(0, nrow(phi));
                    pr_topic[prob_topic$topic] = temp;
                }
            }
            if(length(pr_topic) != nrow(phi)) stop("length(pr_topic) != nrow(phi)");
            
            phi = pr_topic * phi;
        }
        
        phi = normalize_sumToOne2D(phi, 2);
        temp = selectColumn_agg_sum(phi, corpus_pos, corpus$term, corpus$weight);
        temp = t(normalize_sumToOne2D(temp, 2));
        out[length(oldIDs)+(1:length(newIDs)),] = temp;
    }
    return(out);
}

MC_predict.R <- function(
    factor, obs.train, obs.test, corpus, feature, param, x_dyad.train, x_dyad.test,
    nSamples, nBurnIn=1, useTopicProb=FALSE, transduction=FALSE,
    use.C=c(),
    debug=0, verbose=0
){
    if(useTopicProb) stop("useTopicProb not supported");
    if(!transduction) stop("transduction = FALSE is not supported");
    feature$x_dyad = x_dyad.train; # just for syncheck.LDA_RLFM.spec
    size = syncheck.LDA_RLFM.spec(factor=factor, obs=obs.train, corpus=corpus, feature=feature, param=param, 
                                  is.corpus_topic.matrix=useTopicProb);
    
    if(length(param$var_y) == 1) param$var_y = rep(param$var_y, size$nObs);
    if(length(param$var_alpha) == 1) param$var_alpha = rep(param$var_alpha, size$nUsers);
    if(length(param$var_beta) == 1) param$var_beta = rep(param$var_beta, size$nItems);
    if(length(param$var_gamma) == 1) param$var_gamma = rep(param$var_gamma, size$nUsers);
    if(length(param$var_u) == 1) param$var_u = rep_matrix(diag(param$var_u, size$nFactors), size$nUsers);
    if(length(param$var_v) == 1) param$var_v = rep_matrix(diag(param$var_v, size$nFactors), size$nItems);
    if(length(param$var_s) == 1) param$var_s = rep_matrix(diag(param$var_s, size$nTopics), size$nUsers);
    
    nTestCases = nrow(obs.test);
        
    if(max(obs.test$item) > size$nItems) stop("max(obs.test$item) > size$nItems");
    if(max(obs.test$user) > size$nUsers) stop("max(obs.test$user) > size$nUsers");
    if(!is.integer(obs.test$item)) stop("!is.integer(obs.test$item)");
    if(!is.integer(obs.test$user)) stop("!is.integer(obs.test$user)");
    if(nrow(x_dyad.test) != nTestCases)  stop("nrow(x_dyad.test) != nTestCase");
    if(ncol(x_dyad.test) != ncol(x_dyad.train)) stop("ncol(x_dyad.test) != ncol(x_dyad.train)");
        
    xb.train = x_dyad.train %*% param$b;
    xb.test  = x_dyad.test  %*% param$b;
    g0x_user = NULL; c0x_user = NULL; Gx_user = NULL; Hx_user = NULL;
    if(!is.null(param$g0)) g0x_user = feature$x_user %*% param$g0;
    if(!is.null(param$c0)) c0x_user = feature$x_user %*% param$c0;
    if(!is.null(param$G))  Gx_user  = feature$x_user %*% param$G;
    if(!is.null(param$H))  Hx_user  = feature$x_user %*% param$H;
    
    d0x_item = NULL;  Dx_item = NULL;
    if(!is.null(param$d0)) d0x_item = feature$x_item %*% param$d0;
    if(!is.null(param$D))  Dx_item  = feature$x_item %*% param$D;

    # Checks
    if(is.null(g0x_user) && size$nVar_alpha > 0) stop("error!");
    if(is.null(c0x_user) && size$nVar_gamma > 0) stop("error!");
    if(is.null(Gx_user)  && size$nVar_u > 0)     stop("error!");
    if(is.null(Hx_user)  && size$nVar_s > 0)     stop("error!");
    if(is.null(d0x_item) && size$nVar_beta > 0)  stop("error!");
    if(is.null(Dx_item)  && size$nVar_v > 0)     stop("error!");
    
    phi = NULL;
    if(size$nTopics > 0){
        phi = matrix(0.0, nrow=size$nTopics, ncol=size$nTerms);
    }
    
    if(is.null(factor$alpha)) stop("is.null(factor$alpha)");
    if(is.null(factor$beta)) stop("is.null(factor$beta)");
    if(is.null(factor$gamma)) stop("is.null(factor$gamma)");
    
    all.alpha = array(NA, dim=c(size$nUsers, nSamples));
    all.beta  = array(NA, dim=c(size$nItems, nSamples));
    all.gamma = array(NA, dim=c(size$nUsers, nSamples));
    if(size$nFactors > 0){
        all.u     = array(NA, dim=c(size$nUsers, size$nFactors, nSamples));
        all.v     = array(NA, dim=c(size$nItems, size$nFactors, nSamples));
    }
    if(size$nTopics > 0){
    
        all.pred  = array(NA, dim=c(nTestCases, nSamples));
        
        all.s     = array(NA, dim=c(size$nUsers, size$nTopics,  nSamples));
        all.z_avg = array(NA, dim=c(size$nItems, size$nTopics,  nSamples));
        all.corpus_topic = array(NA, dim=c(size$corpusSize, nSamples));
        all.cnt_item_topic = array(NA, dim=c(size$nItems, size$nTopics, nSamples));
        all.cnt_topic_term = array(NA, dim=c(size$nTopics, size$nTerms, nSamples));
        all.cnt_item       = array(NA, dim=c(size$nItems, nSamples));
        all.cnt_topic      = array(NA, dim=c(size$nTopics, nSamples));
        
        z_avg = get_z_avg(corpus, factor$corpus_topic, size);
    }

    for(iter in (1-nBurnIn):nSamples){
        
        if(size$nVar_alpha > 0){
            if("alpha" %in% use.C) factor$alpha = draw_alpha.C(factor, obs.train, param, xb.train, g0x_user, z_avg, size)$sample
            else                   factor$alpha = draw_alpha(factor, obs.train, param, xb.train, g0x_user, z_avg, size)$sample;
        }

        if(size$nVar_beta > 0){
            if("beta" %in% use.C) factor$beta  = draw_beta.C( factor, obs.train, param, xb.train, d0x_item, z_avg, size)$sample
            else                  factor$beta  = draw_beta( factor, obs.train, param, xb.train, d0x_item, z_avg, size)$sample;
        }
        
        if(size$nVar_gamma > 0){
            if("gamma" %in% use.C) factor$gamma = draw_gamma.C(factor, obs.train, param, xb.train, c0x_user, z_avg, size)$sample
            else                   factor$gamma = draw_gamma(factor, obs.train, param, xb.train, c0x_user, z_avg, size)$sample;
        }
        
        if(size$nVar_u > 0 && size$nFactors > 0){
            if("u" %in% use.C) factor$u = draw_u.C(factor, obs.train, param, xb.train, Gx_user, z_avg, size)$sample
            else               factor$u = draw_u(factor, obs.train, param, xb.train, Gx_user, z_avg, size)$sample;
        }
        
        if(size$nVar_v > 0 && size$nFactors > 0){
            if("v" %in% use.C) factor$v = draw_v.C(factor, obs.train, param, xb.train, Dx_item, z_avg, size)$sample
            else               factor$v = draw_v(factor, obs.train, param, xb.train, Dx_item, z_avg, size)$sample;
        }
        
        if(size$nVar_s > 0 && size$nTopics > 0){
            if("s" %in% use.C) factor$s = draw_s.C(factor, obs.train, param, xb.train, Hx_user, z_avg, size)$sample
            else               factor$s = draw_s(factor, obs.train, param, xb.train, Hx_user, z_avg, size)$sample;
        }
        
        if(size$nTopics > 0){
            if("topic" %in% use.C) factor$corpus_topic = draw_topic.C(factor, obs.train, param, xb.train, corpus, size)$corpus_topic
            else                   factor$corpus_topic = draw_topic(factor, obs.train, param, xb.train, corpus, size)$corpus_topic;
            z_avg = get_z_avg(corpus, factor$corpus_topic, size);
        }
        
        if(iter <= 0) next; # Burn-in
        
        all.alpha[,iter] = factor$alpha;
        all.beta[,iter]  = factor$beta;
        all.gamma[,iter] = factor$gamma;
        if(size$nFactors > 0){
            all.u[,,iter]    = factor$u;
            all.v[,,iter]    = factor$v;
        }
        if(size$nTopics > 0){
            all.s[,,iter]     = factor$s;
            all.z_avg[,,iter] = z_avg;
            all.corpus_topic[,iter] = factor$corpus_topic;
            tc = getTopicCounts.R(corpus, factor$corpus_topic, length(factor$beta), ncol(factor$s), max(corpus$term));
            all.cnt_item_topic[,,iter] = tc$cnt_item_topic;
            all.cnt_topic_term[,,iter] = tc$cnt_topic_term;
            all.cnt_topic[,iter] = tc$cnt_topic;
            all.cnt_item[,iter]  = tc$cnt_item;
        }
        
        all.pred[,iter] = xb.test * factor$gamma[obs.test$user]+ factor$alpha[obs.test$user] + factor$beta[obs.test$item];
        if(size$nFactors > 0) all.pred[,iter] = all.pred[,iter] + sum_margin(factor$u[obs.test$user,] * factor$v[obs.test$item,], 1);
        if(size$nTopics  > 0) all.pred[,iter] = all.pred[,iter] + sum_margin(factor$s[obs.test$user,] *    z_avg[obs.test$item,], 1);
    }
    
    out = list();
    out$mean = list(
        alpha = apply(all.alpha, 1, mean),
        beta  = apply(all.beta,  1, mean),
        gamma = apply(all.gamma, 1, mean)
    );
    if(size$nFactors > 0){
        out$mean$u = apply(all.u, c(1,2), mean);
        out$mean$v = apply(all.v, c(1,2), mean);
    }
    if(size$nTopics > 0){
        out$mean$s = apply(all.s, c(1,2), mean);
        out$mean$z_avg = apply(all.z_avg, c(1,2), mean);
        out$mean$corpus_topic = factor$corpus_topic
    }
    out$pred.y = apply(all.pred, 1, mean);
    out$rmse   = sqrt(mean( (obs.test$y - out$pred.y)^2 ));
    out$mae    = mean( abs(  obs.test$y - out$pred.y) );
    
    return(out);
}
