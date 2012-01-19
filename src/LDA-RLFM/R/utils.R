### Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### Re-index users, items and terms
### INPUT: obs    = data.frame(user_id, item_id, y);
###        x_dyad = data.frame(user_id, item_id, feature1, feature2, ...)
###        x_user = data.frame(user_id, feature1, feature2, ...)
###        x_item = data.frame(item_id, feature1, feature2, ...)
###        corpus = data.frame(item_id, term, weight);
### 
reindexData <- function(
    obs, x_dyad=NULL, x_user=NULL, x_item=NULL, corpus=NULL,
    add.intercept=TRUE,          # whether to add an intercept
    UserIDs=sort(unique(obs$user_id)), # selected UserIDs; out$x_user[i,] will correspond to userIDs[i]
    ItemIDs=sort(unique(obs$item_id)), # selected ItemIDs; out$x_item[i,] will correspond to ItemIDs[i]
    TermLevels=NULL              # selected TermLevels if corpus$term is of type factor or character
){
    if(is.null(obs$user_id)) stop("user_id must be a column in obs");
    if(is.null(obs$item_id)) stop("item_id must be a column in obs");
    if(is.null(obs$y)) stop("y must be a column in obs");
    if(!is.null(x_dyad) && is.null(x_dyad$user_id)) stop("user_id must be a column in x_dyad");
    if(!is.null(x_dyad) && is.null(x_dyad$item_id)) stop("item_id must be a column in x_dyad");
    if(!is.null(x_user) && is.null(x_user$user_id)) stop("user_id must be a column in x_user");
    if(!is.null(x_item) && is.null(x_item$item_id)) stop("item_id must be a column in x_item");
    if(!is.null(corpus) && is.null(corpus$item_id)) stop("item_id must be a column in corpus");
    if(!is.null(corpus) && is.null(corpus$term)) stop("term must be a column in corpus");
    if(!is.null(corpus) && is.null(corpus$weight)) stop("weight must be a column in corpus");
    
    out = list();
    out$obs = data.frame(user=match(obs$user_id, UserIDs),
                         item=match(obs$item_id, ItemIDs),
                         y=as.double(obs$y));
    obs.selected = !is.na(out$obs$user) & !is.na(out$obs$item);
    out$obs = out$obs[obs.selected,]
    out$feature = list();
    if(add.intercept){
        regFormula = formula(~.);
    }else{
        regFormula = formula(~.-1);
    }
    if(is.null(x_dyad)){
        out$feature$x_dyad = matrix(1.0, nrow=nrow(out$obs), ncol=1);
    }else{
        if(any(obs$user_id != x_dyad$user_id)) stop("obs$user_id != x_dyad$user_id");
        if(any(obs$item_id != x_dyad$item_id)) stop("obs$item_id != x_dyad$item_id");
        out$feature$x_dyad = model.matrix(regFormula, x_dyad[,-match(c("user_id","item_id"), names(x_dyad)),drop=FALSE]);
        out$feature$x_dyad = out$feature$x_dyad[obs.selected,]
    }
    if(is.null(x_user)){
        out$feature$x_user = matrix(1.0, nrow=length(UserIDs), ncol=1);
    }else{
        select = match(UserIDs, x_user$user_id);
        if(any(is.na(select))) stop("Some user IDs are not found in x_user");
        out$feature$x_user = model.matrix(regFormula, x_user[select,-match(c("user_id"), names(x_user)),drop=FALSE]);
    }
    if(is.null(x_item)){
        out$feature$x_item = matrix(1.0, nrow=length(ItemIDs), ncol=1);
    }else{
        select = match(ItemIDs, x_item$item_id);
        if(any(is.na(select))) stop("Some item IDs are not found in x_item");
        out$feature$x_item = model.matrix(regFormula, x_item[select,-match(c("item_id"), names(x_item)),drop=FALSE]);
    }
    if(!is.null(corpus)){
        if(is.integer(corpus$term)){
            if(!is.null(TermLevels)) stop("corpus$term is an integer array; you should not specify TermLevels!");
            term=corpus$term;
        }else if(is.factor(corpus$term)){
            if(is.null(TermLevels)){
                term = as.integer(corpus$term);
                out$termLevels = levels(corpus$term);
            }else{
                term = match(corpus$term, TermLevels);
                out$termLevels = TermLevels;
            }
        }else if(is.character(corpus$term)){
            if(is.null(TermLevels)){
                out$termLevels = unique(corpus$term);
            }else{
                out$termLevels = TermLevels;
            }
            term = match(corpus$term, out$termLevels);
        }else stop("Unknown type for corpus$term.");
        out$corpus = data.frame(item=match(corpus$item_id, ItemIDs),
                                term=term,
                                weight=as.double(corpus$weight));
        if(any(is.na(out$corpus$term))){
            out$corpus = out$corpus[!is.na(out$corpus$term),];
        }
        if(any(is.na(out$corpus$item))){
            out$corpus = out$corpus[!is.na(out$corpus$item),];
        }
    }
    out$userIDs = UserIDs;
    out$itemIDs = ItemIDs;
    
    return(out);
}

indexTestData <- function(
    data.train,
    obs, x_dyad=NULL, x_user=NULL, x_item=NULL, corpus=NULL,
    add.intercept=TRUE,          # whether to add an intercept
    TermLevels=NULL
){
    UserIDs = c(data.train$userIDs, setdiff(unique(obs$user_id), data.train$userIDs));
    ItemIDs = c(data.train$itemIDs, setdiff(unique(obs$item_id), data.train$itemIDs));
    data.test = reindexData(obs, x_dyad=x_dyad, x_user=x_user, x_item=x_item, corpus=corpus, 
                            UserIDs=UserIDs, ItemIDs=ItemIDs, TermLevels=TermLevels, add.intercept=add.intercept);
    if(!is.null(corpus)){
        data.test$corpus = data.test$corpus[data.test$corpus$term %in% data.train$corpus$term,];
    }
    return(data.test);
}

###
### Check whehter some users or items do not have any rating or some items
### do not have any terms.
###
check.ObsAndCorpus <- function(obs, corpus, nUsers, nItems, nTerms, error=stop){
    temp = sort(unique(obs$user));
    if(length(temp) != nUsers || any(temp != 1:nUsers)) error("Some users do not have observed ratings.");
    temp = sort(unique(obs$item));
    if(length(temp) != nItems || any(temp != 1:nItems)) error("Some items do not have observed ratings.");
    if(!is.null(corpus)){
        temp = sort(unique(corpus$item));
        if(length(temp) != nItems || any(temp != 1:nItems)) error("Some items do not have terms in the corpus.");
        temp = sort(unique(corpus$term));
        if(length(temp) != nTerms || any(temp != 1:nTerms)){
            cat("\n\nWARNING: Some terms do not belong to any items in the corpus.\n\n");
            warning("Some terms do not belong to any items in the corpus.");
        }
    }
}


###
### Checks for user_id, item_id encoding, features and terms
###
### Looks at data.train${userIDs, itemIDs, feature, termLevels}
###
check.compatible <- function(data.train, data.test, strict=FALSE){
    if(strict || (!is.null(data.train$userIDs) && !is.null(data.test$userIDs))){
        len = min(length(data.train$userIDs), length(data.test$userIDs));
        if(any(data.train$userIDs[1:len] != data.test$userIDs[1:len])) stop("The two datasets have incompatible user ids");
    }
    if(strict || (!is.null(data.train$itemIDs) && !is.null(data.test$itemIDs))){
        len = min(length(data.train$itemIDs), length(data.test$itemIDs));
        if(any(data.train$itemIDs[1:len] != data.test$itemIDs[1:len])) stop("The two datasets have incompatible item ids");
    }
    check.feature.comp(data.train$feature$x_dyad, data.test$feature$x_dyad, "x_dyad", strict);
    check.feature.comp(data.train$feature$x_user, data.test$feature$x_user, "x_user", strict);
    check.feature.comp(data.train$feature$x_item, data.test$feature$x_item, "x_item", strict);
    if(strict || (!is.null(data.train$termLevels) && !is.null(data.test$termLevels))){
        len = min(length(data.train$termLevels), length(data.test$termLevels));
        if(any(data.train$termLevels[1:len] != data.test$termLevels[1:len])) stop("The two datasets have incompatible terms");
    }
    if(strict || !is.null(data.train$corpus) || !is.null(data.test$corpus)){
        if(is.null(data.train$corpus) || is.null(data.test$corpus)) stop("one of the corpus is null!");
        t1 = unique(data.train$corpus$term);
        t2 = unique(data.test$corpus$term);
        if(!all(t2 %in% t1)) stop("Some terms in the test set do not exist in the training set!");
    }
}
check.feature.comp <- function(x1, x2, name, strict=FALSE){
    if(strict || !is.null(x1) || !is.null(x2)){
        if(is.null(x1) || is.null(x2)) stop("One dataset has ",name,"; one doesn't");
        name1 = dimnames(x1)[[2]];
        name2 = dimnames(x2)[[2]];
        if(strict || !is.null(name1) || !is.null(name2)){
            if(length(name1) != length(name2) || any(name1 != name2))
                stop("Feature ",name," not compatible!");
        }
    }
}

###
### Perform ridge regression and select the best lambda based on
### test-set MSE
### Y: a vector specifying the response
### X: design matrix (should include the intercept if you want it)
###
lm.ridge.bestMSE <- function(
    Y, X, lambda=c(0,10^(-5:3)),
    fracTest=0.2, fracTrain=NULL, nRepeat=3, nRefine=10, lambda.max=1e8, final=TRUE
){
    if(fracTest >= 1 || fracTest <= 0) stop("fracTest = ",fracTest);
    if(is.null(fracTrain)) fracTrain = 1 - fracTest;
    if(fracTrain + fracTest > 1) stop("fracTrain + fracTest = ", fracTrain + fracTest);
    if(nRepeat < 1) stop("nRepeat = ",nRepeat);
    if(floor(length(Y)*fracTest) < 2) stop("floor(length(Y)*fracTest) = ",floor(length(Y)*fracTest));
    if(floor(length(Y)*fracTrain) < 2) stop("floor(length(Y)*fracTrain) = ",floor(length(Y)*fracTrain));
    if(length(Y) != nrow(X)) stop("length(Y) != nrow(X)");
    lambda = sort(lambda);
    
    mse = rep(0, length(lambda));
    for(i in 1:nRepeat){
        rnd = (1:length(Y))[order(runif(length(Y)))];
        index.train = rnd[1:floor(length(Y)*fracTrain)];
        index.test  = rnd[(length(Y)-floor(length(Y)*fracTest)):length(Y)];
        X.train = X[index.train,];
        Y.train = Y[index.train];
        X.test = X[index.test,];
        Y.test = Y[index.test];
        model = lm.ridge(Y.train ~ X.train - 1, lambda=lambda);
        pred = X.test %*% t(coef(model));
        # browser();
        mse = mse + drop(apply((matrix(rep(Y.test, length(lambda)), ncol=length(lambda)) - pred)^2, 2, sum))/length(Y.test);
    }
    mse = mse / nRepeat;
    
    index.min = which.min(mse);
    if(nRefine > 1){
        if(index.min == length(lambda)){
            prev.min = lambda[length(lambda)];
            lambda.new = 10^( log10(prev.min) + (1:nRefine)*(log10(lambda.max) - log10(prev.min)) / nRefine );
        }else if(index.min == 1){
            lambda.new = 10^( log10(lambda[2]) - (nRefine:1));
        }else{
            b = index.min - 1;
            e = index.min + 1;
            lambda.new = 10^( log10(lambda[b]) + (1:nRefine)*(log10(lambda[e]) - log10(lambda[b])) / (nRefine+1) );
        }
        # browser();
        ans = lm.ridge.bestMSE(Y, X, lambda=lambda.new, fracTest=fracTest, fracTrain=fracTrain, nRepeat=nRepeat, nRefine=0, final=FALSE);
        lambda = c(lambda, lambda.new);
        mse    = c(mse, ans$mse);
        temp = order(lambda);
        lambda = lambda[temp];
        mse    = mse[temp];
    }
    
    out = list(lambda=lambda, mse=mse);
    if(final){
        lambda.best = lambda[which.min(mse)];
        model = lm.ridge(Y ~ X - 1, lambda=lambda.best);
        out$lambda.best = lambda.best;
        out$coef = coef(model);
    }
    return(out);
}

###
### Draw random sample from Dirichlet(alpha)
###
rdirichlet <- function (n, alpha) 
{
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    return(x/as.vector(sm))
}

###
### Create a deep copy of the input object
###
deepCopy <- function(x){
    if(is.list(x)){
        out = list();
        for(i in 1:length(x)){
            out[[i]] = deepCopy(x[[i]]);
        }
        if(length(out) > 0){
            names(out) = names(x)[1:length(out)];
        }
        return(out);
    }
    # if(is.array(x)) out = array(rep(x), dim=dim(x))
    # else if(is.vector(x)) out = rep(x);
    if(is.integer(x)){
        out = x + as.integer(0);
    }else if(is.numeric(x)){
        out = x + 0;
    }else if(is.null(x)){
        out = NULL;
    }else stop("Type not supported");
    
    return(out);
}

###
### Predict the response using the factors
###
predict.y.from.factors <- function(obs, factor, feature, param){
    user = obs$user; item = obs$item;
    ans = feature$x_dyad %*% param$b;
    if(!is.null(factor$gamma)) ans = ans * factor$gamma[user];
    if(!is.null(factor$alpha)) ans = ans + factor$alpha[user];
    if(!is.null(factor$beta))  ans = ans + factor$beta[item];
    if(!is.null(factor$u))     ans = ans + sum_margin(factor$u[user,,drop=FALSE] * factor$v[item,,drop=FALSE], 1);
    if(!is.null(factor$s))     ans = ans + sum_margin(factor$s[user,,drop=FALSE] * factor$z_avg[item,,drop=FALSE], 1);
    return(ans);
}

###
### Complete data log likelihood with the constant term removed
###
logLikelihood <- function(
    obs, factor, feature, param, corpus
){
    size = syncheck.LDA_RLFM.spec(factor=factor, obs=obs, corpus=corpus, feature=feature, param=param, 
                                  is.corpus_topic.matrix=if(!is.null(factor$corpus_topic)) is.matrix(factor$corpus_topic) else FALSE);
    
    nObs     = size$nObs;
    nUsers   = size$nUsers;
    nItems   = size$nItems;
    nFactors = size$nFactors;
    nTopics  = size$nTopics;
    
    ans = 0;
    err = obs$y - predict.y.from.factors(obs, factor, feature, param);
    ans = ans + sum(err^2 / param$var_y) + (if(length(param$var_y) == 1) nObs * log(param$var_y) else sum(log(param$var_y)));
    if(!is.null(param$var_u)){
        err = factor$u - feature$x_user %*% param$G;
        ans = ans + sum(err^2) / param$var_u + nUsers * nFactors * log(param$var_u);
    }
    if(!is.null(param$var_v)){
        err = factor$v - feature$x_item %*% param$D;
        ans = ans + sum(err^2) / param$var_v + nItems * nFactors * log(param$var_v);
    }
    if(!is.null(param$var_s)){
        err = factor$s - feature$x_user %*% param$H;
        ans = ans + sum(err^2) / param$var_s + nUsers * nTopics * log(param$var_s);
    }
    if(!is.null(param$var_alpha)){
        err = factor$alpha - feature$x_user %*% param$g0;
        ans = ans + sum(err^2) / param$var_alpha + nUsers * log(param$var_alpha);
    }
    if(!is.null(param$var_beta)){
        err = factor$beta - feature$x_item %*% param$d0;
        ans = ans + sum(err^2) / param$var_beta  + nItems * log(param$var_beta);
    }
    if(!is.null(param$var_gamma)){
        err = factor$gamma - feature$x_user %*% param$c0;
        ans = ans + sum(err^2) / param$var_gamma  + nUsers * log(param$var_gamma);
    }
    ans = (1/2) * ans;
    
    # LDA part
    if(!is.null(factor$corpus_topic)){
        count = getTopicCounts(corpus, factor$corpus_topic, nItems, nTopics, size$nTerms);
        ans = ans + compute_LDA_negLL(param$eta,    nTopics, size$nTerms, count$cnt_topic, count$cnt_topic_term);
        ans = ans + compute_LDA_negLL(param$lambda, nItems,  nTopics,     count$cnt_item,  count$cnt_item_topic);
    }
    
    return(-ans);
}
compute_LDA_negLL <- function(eta, K, W, Z_k, Z_kl){
    LL = rep(NA, length(eta));
    for(i in 1:length(eta)){
        ans = lgamma(Z_k + W*eta[i]) - apply(lgamma(Z_kl + eta[i]), 1, sum);
        LL[i] = K * (W * lgamma(eta[i]) - lgamma(W*eta[i])) + sum(ans);
    }
    return(LL);
}


###
### E[logLikelihood] over the factors
###
### NOT YET FINISHED!!!
###
ElogLik <- function(
    obs, factor, feature, param, corpus
){
    size = syncheck.LDA_RLFM.spec(factor=factor, obs=obs, corpus=corpus, feature=feature, param=param);
    
    nObs     = size$nObs;
    nUsers   = size$nUsers;
    nItems   = size$nItems;
    nFactors = size$nFactors;
    nTopics  = size$nTopics;
    
    ###
    ### TODO: the following are not correct and to be implemented.
    ###
    
    ans = 0;
    err = obs$y - predict.y.from.factors(obs, factor, feature, param);
    ans = ans + sum(err^2 / param$var_y) + (if(length(param$var_y) == 1) nObs * log(param$var_y) else sum(log(param$var_y)));
    if(!is.null(param$var_u)){
        err = factor$u - feature$x_user %*% param$G;
        ans = ans + sum(err^2) / param$var_u + nUsers * nFactors * log(param$var_u);
    }
    if(!is.null(param$var_v)){
        err = factor$v - feature$x_item %*% param$D;
        ans = ans + sum(err^2) / param$var_v + nItems * nFactors * log(param$var_v);
    }
    if(!is.null(param$var_s)){
        err = factor$s - feature$x_user %*% param$H;
        ans = ans + sum(err^2) / param$var_s + nUsers * nTopics * log(param$var_s);
    }
    if(!is.null(param$var_alpha)){
        err = factor$alpha - feature$x_user %*% param$g0;
        ans = ans + sum(err^2) / param$var_alpha + nUsers * log(param$var_alpha);
    }
    if(!is.null(param$var_beta)){
        err = factor$beta - feature$x_item %*% param$d0;
        ans = ans + sum(err^2) / param$var_beta  + nItems * log(param$var_beta);
    }
    if(!is.null(param$var_gamma)){
        err = factor$gamma - feature$x_user %*% param$c0;
        ans = ans + sum(err^2) / param$var_gamma  + nUsers * log(param$var_gamma);
    }
    ans = (1/2) * ans;
    
    # LDA part
    if(!is.null(factor$corpus_topic)){
        count = getTopicCounts(corpus, factor$corpus_topic, nItems, nTopics, size$nTerms);
        ans = ans + compute_LDA_negLL(param$eta,    nTopics, size$nTerms, count$cnt_topic, count$cnt_topic_term);
        ans = ans + compute_LDA_negLL(param$lambda, nItems,  nTopics,     count$cnt_item,  count$cnt_item_topic);
    }
    
    return(-ans);
}

###
### Add noise to the data$factor and data$param
###
add.noise <- function(data, factor.sd, factor.prob, param.coef.sd, param.var.sd,issequ=F){
    size = syncheck.LDA_RLFM.spec(factor=data$factor, obs=data$obs, corpus=data$corpus, feature=data$feature, param=data$param);
    data$factor$alpha = data$factor$alpha + rnorm(length(data$factor$alpha), 0, factor.sd);
    data$factor$beta  = data$factor$beta  + rnorm(length(data$factor$beta),  0, factor.sd);
    data$factor$gamma = data$factor$gamma + rnorm(length(data$factor$gamma), 0, factor.sd);
    data$factor$u = data$factor$u + rnorm(length(data$factor$u), 0, factor.sd);
    data$factor$v = data$factor$v + rnorm(length(data$factor$v), 0, factor.sd);
    data$factor$s = if(issequ)data$factor$u else data$factor$s + rnorm(length(data$factor$s), 0, factor.sd);
    select = runif(length(data$factor$corpus_topic), 0, 1) < factor.prob;
    data$factor$corpus_topic[select] = as.integer(floor(runif(sum(select), 1, size$nTopics+0.999999999)));

    data$param$b  = data$param$b  + rnorm(length(data$param$b),  0, param.coef.sd);
    data$param$g0 = data$param$g0 + rnorm(length(data$param$g0), 0, param.coef.sd);
    data$param$d0 = data$param$d0 + rnorm(length(data$param$d0), 0, param.coef.sd);
    data$param$c0 = data$param$c0 + rnorm(length(data$param$c0), 0, param.coef.sd);
    data$param$G  = data$param$G  + rnorm(length(data$param$G),  0, param.coef.sd);
    data$param$D  = data$param$D  + rnorm(length(data$param$D),  0, param.coef.sd);
    data$param$H  = if(issequ) data$param$G else data$param$H  + rnorm(length(data$param$H),  0, param.coef.sd);
    
    data$param$var_y = data$param$var_y + rnorm(length(data$param$var_y), 0, param.var.sd);
    data$param$var_alpha = data$param$var_alpha + rnorm(length(data$param$var_alpha), 0, param.var.sd);
    data$param$var_beta = data$param$var_beta + rnorm(length(data$param$var_beta), 0, param.var.sd);
    data$param$var_gamma = data$param$var_gamma + rnorm(length(data$param$var_gamma), 0, param.var.sd);
    data$param$var_u = data$param$var_u + rnorm(length(data$param$var_u), 0, param.var.sd);
    data$param$var_v = data$param$var_v + rnorm(length(data$param$var_v), 0, param.var.sd);
    data$param$var_s = if(issequ)data$param$var_u else data$param$var_s + rnorm(length(data$param$var_s), 0, param.var.sd);
    
    return(data);
}


###
### Utility function to get topic counters
###
getTopicCounts.R <- function(corpus, corpus_topic, nItems, nTopics, nTerms){
    cnt_item_topic = matrix(0, nrow=nItems,  ncol=nTopics); 
    cnt_topic_term = matrix(0, nrow=nTopics, ncol=nTerms);
    cnt_topic = rep(0, nTopics);
    cnt_item  = rep(0, nItems);
    for(n in 1:nrow(corpus)){
        thisItem   = corpus[n, "item"];
        thisTerm   = corpus[n, "term"];
        thisWeight = corpus[n, "weight"];
        if(is.matrix(corpus_topic)){
            thisTopic  = corpus_topic[n,];
            cnt_item_topic[thisItem,] = cnt_item_topic[thisItem,] + thisWeight*thisTopic;
            cnt_topic_term[,thisTerm] = cnt_topic_term[,thisTerm] + thisWeight*thisTopic;
            cnt_topic = cnt_topic + thisWeight*thisTopic;
        }else{
            thisTopic  = corpus_topic[n];
            cnt_item_topic[thisItem, thisTopic] = cnt_item_topic[thisItem, thisTopic] + thisWeight;
            cnt_topic_term[thisTopic, thisTerm] = cnt_topic_term[thisTopic, thisTerm] + thisWeight;
            cnt_topic[thisTopic] = cnt_topic[thisTopic] + thisWeight;
        }
        cnt_item[thisItem] = cnt_item[thisItem] + thisWeight;
    }
    output = list(
        cnt_item_topic = cnt_item_topic,
        cnt_topic_term = cnt_topic_term,
        cnt_topic = cnt_topic,
        cnt_item = cnt_item);
    return(output);
}

get_z_avg <- function(corpus, corpus_topic, size){
    z_avg = matrix(0, nrow=size$nItems, ncol=size$nTopics);
    for(n in 1:nrow(corpus)){
        thisItem   = corpus[n, "item"];
        thisWeight = corpus[n, "weight"];
        if(is.matrix(corpus_topic)){
            thisTopic  = corpus_topic[n,];
            z_avg[thisItem,] = z_avg[thisItem,] + thisWeight * thisTopic;
        }else{
            thisTopic  = corpus_topic[n];
            z_avg[thisItem, thisTopic] = z_avg[thisItem, thisTopic] + thisWeight;
        }
    }
    for(i in 1:nrow(z_avg)){
        z_avg[i,] = z_avg[i,] / sum(z_avg[i,]);
    }
    return(z_avg);
}


rep_matrix <- function(matrix, num){
    ans = array(NA, dim=c(num, nrow(matrix), ncol(matrix)));
    for(i in 1:num){
        ans[i,,] = matrix;
    }
    return(ans);
}

###
### Syntactic check of model specification
###
###     factor  = list(alpha, beta, gamma, u, v, s, corpus_topic);
###     obs     = data.frame(y, user, item);
###     corpus  = data.frame(item, term, weight);
###     feature = list(x_dyad, x_user, x_item);
###     param   = list(b, g0, d0, c0, G, D, H, var_y, var_alpha, var_beta, var_gamma, var_u, var_v, var_s, lambda, eta);
###
syncheck.LDA_RLFM.spec <- function(
    factor, obs, corpus, feature, param, warning=0, print=FALSE,
     factor.name = c("alpha", "beta", "gamma", "u", "v", "s", "corpus_topic"),
        obs.name = c("y", "user", "item"),
     corpus.name = c("item", "term", "weight"),
    feature.name = c("x_dyad", "x_user", "x_item"),
      param.name = c("b", "g0", "d0", "c0", "G", "D", "H", "var_y", "var_alpha", "var_beta", "var_gamma", "var_u", "var_v", "var_s", "lambda", "eta"),
      is.corpus_topic.matrix = FALSE
){
    if(warning > 1){
        warning.any.not.in(names(factor),  c(factor.name, "z_avg"), "You specified the following unnecessary components in factor: ", print=print);
        warning.any.not.in(names(obs),     obs.name,                "You specified the following unnecessary components in obs: ", print=print);
        warning.any.not.in(names(corpus),  corpus.name,             "You specified the following unnecessary components in corpus: ", print=print);
        warning.any.not.in(names(feature), feature.name,            "You specified the following unnecessary components in feature: ", print=print);
        warning.any.not.in(names(param),   c(param.name, "theta", "phi"), "You specified the following unnecessary components in param: ", print=print);
    }
    if(warning > 2){
        # warning.any.not.in(c("alpha", "beta", "gamma"), names(factor), "The following components in factor are required: ", stop=TRUE, print=print);
        warning.any.not.in(factor.name,  names(factor),  "You did not specify the following components in factor: ", print=print);
        warning.any.not.in(obs.name,     names(obs),     "You did not specify the following components in obs: ", print=print);
        warning.any.not.in(corpus.name,  names(corpus),  "You did not specify the following components in corpus: ", print=print);
        warning.any.not.in(feature.name, names(feature), "You did not specify the following components in feature: ", print=print);
        warning.any.not.in(param.name,   names(param),   "You did not specify the following components in param: ", print=print);
    }
    out = list(
        nObs=get.size(nrow(obs)),
        nUsers=get.size(length(factor[["alpha"]]), length(factor[["gamma"]]), nrow(factor[["u"]]), nrow(factor[["s"]]), nrow(feature[["x_user"]])), 
        nItems=get.size(length(factor[["beta"]]), nrow(factor[["v"]]), nrow(feature[["x_item"]])),
        corpusSize = if(is.corpus_topic.matrix) get.size(nrow(factor[["corpus_topic"]]), nrow(corpus))
                     else                       get.size(length(factor[["corpus_topic"]]), nrow(corpus)), 
        nTerms=if(is.null(corpus)) as.integer(0) else get.size(max(corpus[,"term"])), 
        nCorpusWeights=get.size(length(corpus[,"weight"])),
        nFactors=get.size(ncol(factor[["u"]]),ncol(factor[["v"]]),ncol(param[["G"]]),ncol(param[["D"]])),
        nTopics=get.size(ncol(factor[["s"]]),ncol(factor[["H"]])),
        nDyadicFeatures=get.size(ncol(feature[["x_dyad"]]),length(param[["b"]])),
        nUserFeatures=get.size(ncol(feature[["x_user"]]), length(param[["g0"]]), length(param[["c0"]]), nrow(param[["G"]]), nrow(param[["H"]])),
        nItemFeatures=get.size(ncol(feature[["x_item"]]), length(param[["d0"]]), nrow(param[["D"]])),
        nVar_y=as.integer(length(param[["var_y"]])), 
        nVar_alpha=as.integer(length(param[["var_alpha"]])), 
        nVar_beta=as.integer(length(param[["var_beta"]])), 
        nVar_gamma=as.integer(length(param[["var_gamma"]])), 
        nVar_u=as.integer(length(param[["var_u"]])), 
        nVar_v=as.integer(length(param[["var_v"]])), 
        nVar_s=as.integer(length(param[["var_s"]]))
    );
    if(out$nVar_u > 1) out$nVar_u = out$nUsers;
    if(out$nVar_v > 1) out$nVar_v = out$nItems;
    if(out$nVar_s > 1) out$nVar_s = out$nUsers;
    
    if(!is.null(factor[["alpha"]])) check.individual("factor$alpha",is.double(factor[["alpha"]]),length(factor[["alpha"]]),out$nUsers, 
                                    stopIfNull=c("obs$y"=is.null(obs[["y"]]),"feature$x_user"=is.null(feature[["x_user"]]),"param$g0"=is.null(param[["g0"]])));
    if(!is.null(factor[["beta"]]))  check.individual("factor$beta", is.double(factor[["beta"]]), length(factor[["beta"]]), out$nItems, 
                                    stopIfNull=c("obs$y"=is.null(obs[["y"]]),"feature$x_item"=is.null(feature[["x_item"]]),"param$d0"=is.null(param[["d0"]])));
    if(!is.null(factor[["gamma"]])){
        if(is.null(param[["var_gamma"]])){
            check.individual("factor$gamma",is.double(factor[["gamma"]]),length(factor[["gamma"]]),out$nUsers, 
                             stopIfNull=c("obs$y"=is.null(obs[["y"]])));
        }else{
            check.individual("factor$gamma",is.double(factor[["gamma"]]),length(factor[["gamma"]]),out$nUsers, 
                             stopIfNull=c("obs$y"=is.null(obs[["y"]]),"feature$x_user"=is.null(feature[["x_user"]]),"param$c0"=is.null(param[["c0"]])));
        }
    }
    if(!is.null(factor[["u"]])) check.individual("factor$u",is.double(factor[["u"]]),dim(factor[["u"]]),c(out$nUsers,out$nFactors), 
                                stopIfNull=c("obs$y"=is.null(obs[["y"]]),"feature$x_user"=is.null(feature[["x_user"]]),"factor$v"=is.null(factor[["v"]]),"param$G"=is.null(param[["G"]])));
    if(!is.null(factor[["v"]])) check.individual("factor$v",is.double(factor[["v"]]),dim(factor[["v"]]),c(out$nItems,out$nFactors), 
                                stopIfNull=c("obs$y"=is.null(obs[["y"]]),"feature$x_item"=is.null(feature[["x_item"]]),"factor$u"=is.null(factor[["u"]]),"param$D"=is.null(param[["D"]])));
    if(!is.null(factor[["s"]])) check.individual("factor$s",is.double(factor[["s"]]),dim(factor[["s"]]),c(out$nUsers,out$nTopics),  
                                stopIfNull=c("obs$y"=is.null(obs[["y"]]),"feature$x_user"=is.null(feature[["x_user"]]),"factor$corpus_topic"=is.null(factor[["corpus_topic"]]),"param$H"=is.null(param[["H"]])));
    if(!is.null(factor[["corpus_topic"]])){
        if(is.corpus_topic.matrix){
            check.individual("factor$corpus_topic",is.double(factor[["corpus_topic"]]),dim(factor[["corpus_topic"]]), c(out$corpusSize, out$nTopics), 
                             stopIfNull=c("factor$s"=is.null(factor[["s"]]),"corpus"=(length(corpus)==0),"param$lambda"=is.null(param[["lambda"]]),"param$eta"=is.null(param[["eta"]])));
        }else{
            check.individual("factor$corpus_topic",is.integer(factor[["corpus_topic"]]),length(factor[["corpus_topic"]]),out$corpusSize, 
                             stopIfNull=c("factor$s"=is.null(factor[["s"]]),"corpus"=(length(corpus)==0),"param$lambda"=is.null(param[["lambda"]]),"param$eta"=is.null(param[["eta"]])));
        }
    }
    if(length(obs) > 0){
        warning.any.not.in(obs.name, names(obs), "You did not specify the following components in obs: ", stop=TRUE);
        check.individual("obs$y",is.double(obs[["y"]]),length(obs[["y"]]),out$nObs);
        check.individual("obs$user",is.integer(obs[["user"]]),length(obs[["user"]]),out$nObs);
        check.individual("obs$item",is.integer(obs[["item"]]),length(obs[["item"]]),out$nObs);
    }
    if(length(corpus) > 0){
        warning.any.not.in(c("item", "term"), names(corpus), "You did not specify the following components in corpus: ", stop=TRUE);
        check.individual("corpus$item",is.integer(corpus[["item"]]),length(corpus[["item"]]),out$corpusSize);
        check.individual("corpus$term",is.integer(corpus[["term"]]),length(corpus[["term"]]),out$corpusSize);
        if(!is.null(corpus[["weight"]])) check.individual("corpus$weight",is.double(corpus[["weight"]]),length(corpus[["weight"]]),out$corpusSize);
    }
    if(!is.null(feature[["x_dyad"]])) check.individual("feature$x_dyad",is.double(feature[["x_dyad"]]),dim(feature[["x_dyad"]]),c(out$nObs,out$nDyadicFeatures));
    if(!is.null(feature[["x_user"]])) check.individual("feature$x_user",is.double(feature[["x_user"]]),dim(feature[["x_user"]]),c(out$nUsers,out$nUserFeatures));
    if(!is.null(feature[["x_item"]])) check.individual("feature$x_item",is.double(feature[["x_item"]]),dim(feature[["x_item"]]),c(out$nItems,out$nItemFeatures));

    if(!is.null(param[["b"]]))  check.individual("param$b", is.double(param[["b"]]), length(param[["b"]]), out$nDyadicFeatures);
    if(!is.null(param[["g0"]])) check.individual("param$g0",is.double(param[["g0"]]),length(param[["g0"]]),out$nUserFeatures);
    if(!is.null(param[["d0"]])) check.individual("param$d0",is.double(param[["d0"]]),length(param[["d0"]]),out$nItemFeatures);
    if(!is.null(param[["c0"]])) check.individual("param$c0",is.double(param[["c0"]]),length(param[["c0"]]),out$nUserFeatures);
    if(!is.null(param[["G"]]))  check.individual("param$G", is.double(param[["G"]]), dim(param[["G"]]), c(out$nUserFeatures,out$nFactors));
    if(!is.null(param[["D"]]))  check.individual("param$D", is.double(param[["D"]]), dim(param[["D"]]), c(out$nItemFeatures,out$nFactors));
    if(!is.null(param[["H"]]))  check.individual("param$H", is.double(param[["H"]]), dim(param[["H"]]), c(out$nUserFeatures,out$nTopics));

    if(!is.null(param[["var_y"]]))     check.individual("param$var_y",    is.double(param[["var_y"]]),    length(param[["var_y"]]),    out$nObs,  okDim=1);
    if(!is.null(param[["var_alpha"]])) check.individual("param$var_alpha",is.double(param[["var_alpha"]]),length(param[["var_alpha"]]),out$nUsers,okDim=1);
    if(!is.null(param[["var_beta"]]))  check.individual("param$var_beta", is.double(param[["var_beta"]]), length(param[["var_beta"]]), out$nItems,okDim=1);
    if(!is.null(param[["var_gamma"]])) check.individual("param$var_gamma",is.double(param[["var_gamma"]]),length(param[["var_gamma"]]),out$nUsers,okDim=1);
    if(!is.null(param[["var_u"]]))     check.individual("param$var_u",    is.double(param[["var_u"]]),    first.not.null(dim(param[["var_u"]]),length(param[["var_u"]])), c(out$nUsers,out$nFactors,out$nFactors), okDim=1);
    if(!is.null(param[["var_v"]]))     check.individual("param$var_v",    is.double(param[["var_v"]]),    first.not.null(dim(param[["var_v"]]),length(param[["var_v"]])), c(out$nItems,out$nFactors,out$nFactors), okDim=1);
    if(!is.null(param[["var_s"]]))     check.individual("param$var_s",    is.double(param[["var_s"]]),    first.not.null(dim(param[["var_s"]]),length(param[["var_s"]])), c(out$nUsers,out$nTopics,out$nTopics), okDim=1);
    
    if(!is.null(param[["lambda"]])) check.individual("param$lambda",is.double(param[["lambda"]]),length(param[["lambda"]]),1);
    if(!is.null(param[["eta"]]))    check.individual("param$eta",   is.double(param[["eta"]]),   length(param[["eta"]]),   1);
    
    return(out);
}
warning.any.not.in <- function(a, b, msg, stop=FALSE, print=FALSE){
    if(any(!(a %in% b))){
        temp = a[!(a %in% b)];
        if(stop) stop(msg,paste(temp, collapse=", "))
        else     warning(msg,paste(temp, collapse=", "),call.=FALSE);
        if(print) cat("\nWARNING: ",msg,paste(temp, collapse=", "),"\n\n",sep="");
    }
}
get.size <- function(...){
    temp = c(...);
    if(is.null(temp)) return(as.integer(0));
    s = max(temp);
    return(as.integer(s));
}
check.individual <- function(name, typeOK, dim1, dim2, okDim=NULL, stopIfNull=NULL){
    if(!typeOK) stop(name,"'s data type is not correct!");
    if(!is.null(stopIfNull)){
        temp = names(stopIfNull)[stopIfNull];
        if(length(temp) > 0) stop("When ",name," is specified, the following cannot be NULL: ",paste(temp,collapse=", "));
    }
    if(!is.null(okDim) && length(dim1) == length(okDim) && all(dim1 == okDim)) return();
    if(length(dim1) != length(dim2)) stop(name,"'s dim (or length) is not correct: (",paste(dim1,collapse=","),") vs (",paste(dim2,collapse=","),")");
    if(any(dim1 != dim2)) stop(name,"'s dim (or length) is not correct: (",paste(dim1,collapse=","),") vs (",paste(dim2,collapse=","),")");
}
first.not.null <- function(...){
    temp = list(...);
    for(i in 1:length(temp)){
        if(!is.null(temp[[i]])) return(temp[[i]]);
    }
    return(NULL);
}

###
### Check whether two objects x1 and x2 are different, where x1 and x2 can be lists of lists.
###
is.diff <- function(x1, x2, precision=1e-10, prefix=""){
    if(length(x1) != length(x2)){
        cat(prefix,": Different length! (",length(x1)," vs ",length(x2),")\n",sep="");
        return(TRUE);
    }
    if(length(x1) == 0) return(FALSE);
    if(is.list(x1)){
        if(!is.list(x2)){
            cat(prefix,": Different types! (list vs non-list)\n",sep="");
            return(TRUE);
        }
        name1 = sort(names(x1));
        name2 = sort(names(x2));
        if(is.null(name1) || is.null(name2)){
            if(!is.null(name2) || !is.null(name1)){
                cat(prefix,": One has no names; the other has names!\n",sep="");
                return(TRUE);
            }
            ans = FALSE;
            for(i in 1:length(x1)){
                ret = is.diff(x1[[i]], x2[[i]], precision, prefix=paste(prefix,"[[",i,"]]",sep=""));
                ans = ans || ret;
            }
            return(ans);
        }else{
            if(any(name1 != name2)){
                cat(prefix,": Different names!\n",sep="");
                return(TRUE);
            }
            ans = FALSE;
            for(i in 1:length(name1)){
                name = name1[i];
                ret = is.diff(x1[[name]], x2[[name]], precision, prefix=paste(prefix,"$",name,sep=""));
                ans = ans || ret;
            }
            return(ans);
        }
    }else{
        indexDiff = (1:length(x1))[x1 != x2];
        # print(indexDiff);
        if(length(indexDiff) == 0) return(FALSE);
        
        value = cbind(x1[indexDiff], x2[indexDiff]);
        denom = apply(abs(value), 1, max);
        diff  = abs(value[,1]-value[,2]) / denom;
        
        # print(diff);
        
        temp = (1:length(indexDiff))[diff > precision];
        indexDiff = indexDiff[temp];
        diff      = diff[temp];
        if(length(indexDiff) == 0) return(FALSE);
        
        for(i in 1:length(indexDiff)){
            index = indexDiff[i];
            cat(prefix,"[",index,"]: ",x1[index]," vs ",x2[index]," (diff=",diff[i],")\n",sep="");
        }
        return(TRUE);
    }
}

###
### Draw a random sample from a multivariate normal distribtuion
###     This function is implemented to make sure the C implementation generates exactly
###     the same output as the R implementation
###
my_rmvnorm <- function(n, mu, Sigma=NULL, Sigma.inv=NULL, debug=10, tol=1e-8){
    nDim <- length(mu)
    if (!all(dim(Sigma) == c(nDim, nDim))) 
        stop("incompatible arguments")
    if(is.null(Sigma.inv)){
        Sigma.inv = solve(Sigma);
        dim.names = dimnames(Sigma);
    }else if(is.null(Sigma)){
        dim.names = dimnames(Sigma.inv);
    }else{
        error("Please specify one of Sigma or Sigma.inv (not both)");
    }
    eigen <- sym_eigen(Sigma.inv);
    eigen.value <- 1/eigen$values;
        
    if(debug >= 3){
        if(is.null(Sigma)) Sigma = solve(Sigma.inv);
        if(max(abs(eigen$vectors %*% diag(eigen.value, nDim) %*% t(eigen$vectors) - Sigma)) > tol * abs(eigen.value[1]))
            stop("sym_eigen(Sigma) seems to have some problems!!");
    }
    
    if (!all(eigen.value >= -tol * abs(eigen.value[1]))) stop("'Sigma' is not positive definite");
    Z <- matrix(rnorm(nDim * n), n);
    
    # cat("eigenvalue: ");   print(drop(eigen.value));
    # cat("eigenvector:\n"); print(eigen$vectors, nDim);
    # cat("temp:\n");
    # print(eigen$vectors %*% diag(sqrt(pmax(eigen.value, 0)), nDim));
    # cat("rnd: ");
    # print(drop(Z));
    
    out <- drop(mu) + eigen$vectors %*% diag(sqrt(pmax(eigen.value, 0)), nDim) %*% t(Z);
    mu.names <- names(mu);
    if (is.null(mu.names) && !is.null(temp <- dim.names)) mu.names <- temp[[1]];
    dimnames(out) <- list(mu.names, NULL);
    if(n == 1) return(drop(out))
    else       return(t(out));
}

# Get multivariate normal sample
#   mean[k,] is the mean vector for the kth sample point
#   var[k,,] is the variance-covariance matrix for the kth sample point
# output[k,] is the kth sample point
getMVNSample <- function(mean, var=NULL, var.inv=NULL, FUN=my_rmvnorm){
    
    nPoints = nrow(mean);
    nDim    = ncol(mean);
    if(!is.null(var)) temp = dim(var)
    else temp = dim(var.inv) ;
    if(temp[1] != nPoints || temp[2] != nDim || temp[3] != nDim) stop("size mismatch");
    
    if(nDim == 1) return(matrix(rnorm(nPoints, mean, sqrt(var)), nrow=nPoints, ncol=1));
    
    output = matrix(NA, nrow=nPoints, ncol=nDim);
    for(k in 1:nPoints){
        if(!is.null(var.inv)){
            output[k,] = FUN(1, mu=mean[k,], Sigma.inv=var.inv[k,,]);
        }else{
            output[k,] = FUN(1, mu=mean[k,], Sigma=var[k,,]);
        }
    }
    return(output);
}

subsample.ROC <- function(perf, nPoints=1000){
    n = length(perf@x.values[[1]]);
    size = floor(n / nPoints);
    select = ((1:n) %% size) == 1;
    select[n] = TRUE;
    perf@x.values[[1]] = perf@x.values[[1]][select];
    perf@y.values[[1]] = perf@y.values[[1]][select];
    perf@alpha.values[[1]] = perf@alpha.values[[1]][select];
    return(perf);
}
