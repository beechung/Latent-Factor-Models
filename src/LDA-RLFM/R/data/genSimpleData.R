### Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

# Generate some normal data for testing
#   Rating/Observation for user i, item j
#     y[i,j] ~ N( (x_dyad[i,j]'b)*gamma[i] + alpha[i] + beta[j] + u[i]'v[j] + s[i]'z_avg[j], var_y)
#   alpha[i] ~ N(g0'x_user[i], var_alpha)
#    beta[j] ~ N(d0'x_item[j], var_beta)
#   gamma[i] ~ N(c0'x_user[i], var_gamma)
#       s[i] ~ N(H'x_user[i],  var_s)
#       u[i] ~ N(G'x_user[i],  var_u)
#       v[j] ~ N(D'x_item[j],  var_v)
#   z_avg[j] = topic distribution averaged over terms
#
#   All feature values (x_dyad, x_user, x_item) follow N(0,1)
#   
#   output$obs is a data.frame with columns (user, item, rating)
#                                              i     j      y
#   output$feature$x_dyad[k,] is the feature vector for the kth observation
#   output$feature$x_user[i,] is the feature vector for user i
#   output$feature$x_item[j,] is the feature vector for item j
#   
#       x_user[output$obs[k,"User"],] gives the user feature vector for the kth observation
#       x_item[output$obs[k,"Item"],] gives the item feature vector for the kth observation
#       x_dyad[k,] gives the dyadic feature vector fot the kth observation
#
#   output$factor$alpha[i]: The main effect for user i
#   output$factor$beta[j]:  The main effect for item j
#   output$factor$gamma[i]: The main effect for user i
#   output$factor$u[i,]:    The factor vector for user i
#   output$factor$s[i,]:    The factor vector for user i
#   output$factor$v[j,]:    The factor vector for item j
#   output$factor$corpus_topic[n]: The topic of the nth term in the corpus
#
#   output$corpus = data.frame(item, term, weight)
#
#   output$param = a copy of the input parameters
#
genNormalData <- function(
    nUsers, nItems, nObs, nTerms, corpusSize,
    b, var_y,
    g0, d0, c0, var_alpha, var_beta, var_gamma,
    G,  D,  H,  var_u,     var_v,    var_s,
    eta, lambda
){
    if(!is.vector(b))  stop("b should be a vector");
    if(!is.vector(g0)) stop("g0 should be a vector");
    if(!is.vector(d0)) stop("d0 should be a vector");
    if(!is.null(c0) && !is.vector(c0)) stop("c0 should be a vector");
    if(!is.null(G) && !is.matrix(G))  stop("G should be a matrix");
    if(!is.null(D) && !is.matrix(D))  stop("D should be a matrix");
    if(!is.matrix(H))  stop("H should be a matrix");
    
    nDyadicFeatures = length(b);
    nUserFeatures   = length(g0);
    nItemFeatures   = length(d0);
    nFactors        = if(is.null(G)) 0 else ncol(G);
    nTopics         = ncol(H);
    
    if(nFactors > 0 && ncol(D) != nFactors) stop("ncol(D) != nFactors");
    if(!is.null(G) && nrow(G) != nUserFeatures) stop("nrow(G) != nUserFeatures");
    if(nrow(H) != nUserFeatures) stop("nrow(H) != nUserFeatures");
    if(!is.null(D) && nrow(D) != nItemFeatures) stop("nrow(D) != nItemFeatures");
    if(!is.null(c0) && length(c0) != nUserFeatures) stop("length(c0) != nUserFeatures");
    if(nObs < nUsers || nObs < nItems) stop("nObs < nUsers || nObs < nItems");
    if(corpusSize < nTerms || corpusSize < nItems) stop("corpusSize < nTerms || corpusSize < nItems");

    ###
    ### Sample topics 
    ###
    cps_item = c(floor(runif(corpusSize-nItems, min=1, max=nItems+1)), 1:nItems);
    cps_weight = runif(corpusSize, min=0.2, max=2);
    
    # theta: nItems x nTopics
    # phi: nTopics x nTerms
    theta = rdirichlet(nItems,  rep(lambda, nTopics));
    phi   = rdirichlet(nTopics, rep(eta,    nTerms));

    z_avg = matrix(NA, nrow=nItems, ncol=nTopics);
    cps_topic = rep(NA, corpusSize);
    
    index = tapply(1:corpusSize, cps_item, c, simplify=F);
    if(any(names(index) != 1:nItems)) stop("error: any(names(index) != 1:nItems)");
    
	for(j in 1:length(index)){
		select = index[[j]];
        temp = rmultinom(1,length(select),theta[j,]);
        z_avg[j,] = temp / length(select);
        cps_topic[select] = indexWithQuantities(temp);
	}
    
    cps_term = rep(NA, corpusSize);

    index = tapply(1:corpusSize, cps_topic, c, simplify=F);
	for(k in 1:length(index)){
		select = index[[k]];
        temp = indexWithQuantities(rmultinom(1,length(select),phi[k,]));
        cps_term[select] = temp[order(runif(length(temp)))];
	}
    
    x_dyad = matrix(rnorm(nObs*nDyadicFeatures), nrow=nObs,   ncol=nDyadicFeatures,dimnames=list(NULL, sprintf("x_dyad_%03d", 1:nDyadicFeatures)));
    x_user = matrix(rnorm(nUsers*nUserFeatures), nrow=nUsers, ncol=nUserFeatures,  dimnames=list(NULL, sprintf("x_user_%03d", 1:nUserFeatures)));
    x_item = matrix(rnorm(nItems*nItemFeatures), nrow=nItems, ncol=nItemFeatures,  dimnames=list(NULL, sprintf("x_item_%03d", 1:nItemFeatures)));
    
    alpha = rnorm(nUsers, mean=x_user %*% g0, sd=sqrt(var_alpha));
    beta  = rnorm(nItems, mean=x_item %*% d0, sd=sqrt(var_beta));
    if(is.null(c0)){
        gamma = rep(1.0, nUsers);
    }else{
        gamma = rnorm(nUsers, mean=x_user %*% c0, sd=sqrt(var_gamma));
    }
    
    if(!is.null(G)) u = x_user %*% G + rnorm(nUsers*nFactors, mean=0, sd=sqrt(var_u))
    else            u = matrix(0.0, nrow=nUsers, ncol=1);
    
    if(!is.null(D)) v = x_item %*% D + rnorm(nItems*nFactors, mean=0, sd=sqrt(var_v))
    else            v = matrix(0.0, nrow=nItems, ncol=1);
    
    s = x_user %*% H + rnorm(nUsers*nTopics,  mean=0, sd=sqrt(var_s));
    
    user = c(1:nUsers, floor(runif(nObs-nUsers, min=1, max=nUsers+1)));
    item = c(floor(runif(nObs-nItems, min=1, max=nItems+1)), 1:nItems);
        
    y = (x_dyad %*% b) * gamma[user] + alpha[user] + beta[item] + 
        apply(u[user,,drop=FALSE] * v[item,,drop=FALSE], 1, sum) + 
        apply(s[user,,drop=FALSE] * z_avg[item,,drop=FALSE], 1, sum) +
        rnorm(nObs, mean=0, sd=sqrt(var_y));
    
    if(is.null(G)) u = NULL;
    if(is.null(D)) v = NULL;
    
    output=list();
    output$obs     = data.frame(user=as.integer(user), item=as.integer(item), y=y);
    output$corpus  = data.frame(item=as.integer(cps_item), term=as.integer(cps_term), weight=cps_weight);
    output$feature = list(x_user=x_user, x_item=x_item, x_dyad=x_dyad);
    output$factor  = list(alpha=alpha, beta=beta, gamma=gamma, u=u, v=v, s=s, corpus_topic=cps_topic, z_avg=z_avg);
    output$param   = list(b=b, g0=g0, d0=d0, c0=c0, G=G, H=H, D=D, var_y=var_y,
                          var_alpha=var_alpha, var_beta=var_beta, var_gamma=var_gamma,
                          var_u=var_u, var_v=var_v, var_s=var_s, eta=eta, lambda=lambda, theta=theta, phi=phi);
    
    return(output);
}


###
### Split a dataset into train and test
###
split.forTesting <- function(data, nUsers.testOnly, nItems.testOnly, nObs.test){
    size = syncheck.LDA_RLFM.spec(data$factor, data$obs, data$corpus, data$feature, data$param);
    user.test = ((1:size$nUsers)[order(runif(size$nUsers))])[1:nUsers.testOnly];
    item.test = ((1:size$nItems)[order(runif(size$nItems))])[1:nItems.testOnly];
    obs.forTest = (data$obs$user %in% user.test) | (data$obs$item %in% item.test);
    obsID.train = (1:size$nObs)[!obs.forTest];
    obsID.test  = (1:size$nObs)[ obs.forTest];
    need = nObs.test - length(obsID.test);
    if(need < 0){
        warning("nObs.test = ",nObs.test,", but we have already allocated ", length(obsID.test), " for testing because of test-only users and items!");
    }else if(need > 0){
        obs.forTest = ((1:length(obsID.train))[order(runif(length(obsID.train)))])[1:need];
        obsID.test  = c(obsID.test, obsID.train[obs.forTest]);
        obsID.train = obsID.train[-obs.forTest];
    }
    
    userID.train = sort(unique(data$obs$user[obsID.train]));
    userID.test  = sort(unique(data$obs$user[obsID.test]));
    userID.test  = c(userID.train, userID.test[!(userID.test %in% userID.train)]);
    
    itemID.train = sort(unique(data$obs$item[obsID.train]));
    itemID.test  = sort(unique(data$obs$item[obsID.test]));
    itemID.test  = c(itemID.train, itemID.test[!(itemID.test %in% itemID.train)]);
    
    out = list();
    out$train = selectDataSubset(data, obsID.train, userID.train, itemID.train);
    out$test  = selectDataSubset(data, obsID.test,  userID.test,  itemID.test);
    
    size.train = syncheck.LDA_RLFM.spec(out$train$factor, out$train$obs, out$train$corpus, out$train$feature, out$train$param);
    size.test  = syncheck.LDA_RLFM.spec(out$test$factor, out$test$obs, out$test$corpus, out$test$feature, out$test$param);
    
    if(any(sort(unique(out$train$obs$user)) != 1:size.train$nUsers)) stop("error");
    if(any(sort(unique(out$train$obs$item)) != 1:size.train$nItems)) stop("error");
    if(max(out$test$obs$user) != size.test$nUsers) stop("error");
    if(max(out$test$obs$item) != size.test$nItems) stop("error");
    if(size.test$nUsers != size$nUsers) stop("error");
    if(size.test$nItems != size$nItems) stop("error");
    if(size.train$nUsers + nUsers.testOnly > size$nUsers) stop("error");
    if(size.train$nItems + nItems.testOnly > size$nItems) stop("error");
    
    return(out);
}

###
### Select a subset of data by obsID (row number), userID and itemID
###
###                 [1]  [2]  [3]  [4]  [5]  <- New ID in the output data
### E.g. userID = c( 1,   2,   4,   8,   9)  <- Original ID in input data
###                                             (ID in data$user, row number in data$factor$alpha, ...)
###
### obsID: the selected observation IDs (row number in the input data$obs)
###
selectDataSubset <- function(data, obsID, userID, itemID){
    # user.newID[oldID] = newID
    # item.newID[oldID] = newID
    user.newID = c(); item.newID = c();
    user.newID[userID] = as.integer(1:length(userID));
    item.newID[itemID] = as.integer(1:length(itemID));
    newObs = data.frame(
        user = user.newID[data$obs$user[obsID]],
        item = item.newID[data$obs$item[obsID]],
        y    = data$obs$y[obsID]
    );
    corpus.select = data$corpus$item %in% itemID;
    newCorpus = data$corpus[corpus.select,];
    newCorpus$item = item.newID[newCorpus$item];
    newFeature = list(
        x_user=data$feature$x_user[userID,,drop=FALSE], x_item=data$feature$x_item[itemID,,drop=FALSE], x_dyad=data$feature$x_dyad[obsID,,drop=FALSE]
    );
    newFactor = list(
        alpha=data$factor$alpha[userID], beta=data$factor$beta[itemID], gamma=data$factor$gamma[userID], 
        u=data$factor$u[userID,], v=data$factor$v[itemID,], s=data$factor$s[userID,], 
        corpus_topic=data$factor$corpus_topic[corpus.select], z_avg=data$factor$z_avg[itemID,]
    );
    newData = list(
        obs = newObs, corpus = newCorpus, feature = newFeature,
        factor = newFactor, param = data$param
    );
    return(newData);
}
