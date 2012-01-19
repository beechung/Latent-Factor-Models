### Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### MCEM_EStep.R
###
### INPUT:  factor  = list(alpha, beta, gamma, u, v, s, corpus_topic);
###                   Initial factor values
###         obs     = data.frame(y, user, item);
###         corpus  = data.frame(item, term, weight);
###         feature = list(x_dyad, x_user, x_item);
###         param   = list(b, g0, d0, c0, G, D, H, var_y, var_alpha, var_beta, var_gamma, var_u, var_v, var_s, lambda, eta);
###         try     = list(lambda, eta);
###         
### OUTPUT: mean    = list(alpha, beta, gamma, u, v, s, z_avg, gamma2, o_gamma, corpus_topic);  # gamma2 = gamma^2
###         sumvar  = list(alpha, beta, gamma, u, v, s, o_adj);
###         objval  = list(eta, lambda);
###         sampvar = list(alpha, beta, gamma, u, v, s, z_avg);
###
###   isOldUser[i]: [TRUE/FALSE] Whehter the ith user is an old user; default: all FALSE
###                 For old users, we set g0x_user[i] = alpha[i], c0x_user[i] = gamma[i], Gx_user[i,] = u[i,] and Hx_user[i,] = s[i,]
###   isOldItem[j]: [TRUE/FALSE] Whehter the jth item is an old item; default: all FALSE
###                 For old items, we set d0x_item[j] = beta[j], Dx_item[j] = v[j,]
###
MCEM_EStep.R <- function(
    factor, obs, corpus, feature, param, try, nSamples, nBurnIn=1,
    userFactorVar=0, itemFactorVar=0, userTopicVar=0, itemTopicVar=0, drawTopicSample=1,
    isOldUser=NULL, isOldItem=NULL,
    debug=0, verbose=0, use.C=c()
){
    size = syncheck.LDA_RLFM.spec(factor=factor, obs=obs, corpus=corpus, feature=feature, param=param);
    
    if(length(param$var_y) == 1) param$var_y = rep(param$var_y, size$nObs);
    if(length(param$var_alpha) == 1) param$var_alpha = rep(param$var_alpha, size$nUsers);
    if(length(param$var_beta) == 1) param$var_beta = rep(param$var_beta, size$nItems);
    if(length(param$var_gamma) == 1) param$var_gamma = rep(param$var_gamma, size$nUsers);
    if(length(param$var_u) == 1) param$var_u = rep_matrix(diag(param$var_u, size$nFactors), size$nUsers);
    if(length(param$var_v) == 1) param$var_v = rep_matrix(diag(param$var_v, size$nFactors), size$nItems);
    if(length(param$var_s) == 1) param$var_s = rep_matrix(diag(param$var_s, size$nTopics), size$nUsers);
    
    xb  = feature$x_dyad %*% param$b;
    if(is.null(isOldUser)){
        g0x_user = feature$x_user %*% param$g0;
        c0x_user = feature$x_user %*% param$c0;
        Gx_user  = feature$x_user %*% param$G;
        Hx_user  = feature$x_user %*% param$H;
    }else{
        if(length(isOldUser) != size$nUsers) stop("length(isOldUser) != nUsers");
        x_user.new = feature$x_user[!isOldUser,,drop=FALSE];
        g0x_user = factor$alpha;
        c0x_user = factor$gamma;
        Gx_user  = factor$u;
        Hx_user  = factor$s;
        g0x_user[!isOldUser] = x_user.new %*% param$g0;
        c0x_user[!isOldUser] = x_user.new %*% param$c0;
        Gx_user[!isOldUser,] = x_user.new %*% param$G;
        Hx_user[!isOldUser,] = x_user.new %*% param$H;
    }
    
    if(is.null(isOldItem)){
        d0x_item = feature$x_item %*% param$d0;
        Dx_item  = feature$x_item %*% param$D;
    }else{
        if(length(isOldItem) != size$nItems) stop("length(isOldItem) != nItems");
        x_item.new = feature$x_item[!isOldItem,,drop=FALSE];
        d0x_item = factor$beta;
        Dx_item  = factor$v;
        d0x_item[!isOldItem] = x_item.new %*% d0;
        Dx_item[!isOldItem,] = x_item.new %*% D;
    }
    
    all.alpha = array(NA, dim=c(size$nUsers, nSamples));
    all.beta  = array(NA, dim=c(size$nItems, nSamples));
    all.gamma = array(NA, dim=c(size$nUsers, nSamples));
    all.u     = array(NA, dim=c(size$nUsers, size$nFactors, nSamples));
    all.v     = array(NA, dim=c(size$nItems, size$nFactors, nSamples));
    all.s     = array(NA, dim=c(size$nUsers, size$nTopics,  nSamples));
    all.z_avg = array(NA, dim=c(size$nItems, size$nTopics,  nSamples));
    all.corpus_topic = array(NA, dim=c(size$corpusSize, nSamples));
    all.o_gamma = array(NA, dim=c(size$nObs, nSamples));
    all.o = array(NA, dim=c(size$nObs, nSamples));
    
    all.cnt_item_topic = array(NA, dim=c(size$nItems, size$nTopics, nSamples));
    all.cnt_topic_term = array(NA, dim=c(size$nTopics, size$nTerms, nSamples));
    all.cnt_item       = array(NA, dim=c(size$nItems, nSamples));
    all.cnt_topic      = array(NA, dim=c(size$nTopics, nSamples));

    z_avg = get_z_avg(corpus, factor$corpus_topic, size);

    for(iter in (1-nBurnIn):nSamples){
        
        if("alpha" %in% use.C) factor$alpha = draw_alpha.C(factor, obs, param, xb, g0x_user, z_avg, size)$sample
        else                   factor$alpha = draw_alpha(factor, obs, param, xb, g0x_user, z_avg, size)$sample;

        if("beta" %in% use.C) factor$beta  = draw_beta.C( factor, obs, param, xb, d0x_item, z_avg, size)$sample
        else                  factor$beta  = draw_beta( factor, obs, param, xb, d0x_item, z_avg, size)$sample;
        
        if("gamma" %in% use.C) factor$gamma = draw_gamma.C(factor, obs, param, xb, c0x_user, z_avg, size)$sample
        else                   factor$gamma = draw_gamma(factor, obs, param, xb, c0x_user, z_avg, size)$sample;
        
        if("u" %in% use.C) factor$u = draw_u.C(factor, obs, param, xb, Gx_user, z_avg, size)$sample
        else               factor$u = draw_u(factor, obs, param, xb, Gx_user, z_avg, size)$sample;
        
        if("v" %in% use.C) factor$v = draw_v.C(factor, obs, param, xb, Dx_item, z_avg, size)$sample
        else               factor$v = draw_v(factor, obs, param, xb, Dx_item, z_avg, size)$sample;
        
        if("s" %in% use.C) factor$s = draw_s.C(factor, obs, param, xb, Hx_user, z_avg, size)$sample
        else               factor$s = draw_s(factor, obs, param, xb, Hx_user, z_avg, size)$sample;
        
        if("topic" %in% use.C) factor$corpus_topic = draw_topic.C(factor, obs, param, xb, corpus, size)$corpus_topic
        else                   factor$corpus_topic = draw_topic(factor, obs, param, xb, corpus, size)$corpus_topic;
        
        z_avg = get_z_avg(corpus, factor$corpus_topic, size);
        
        # cat("s = \n"); print(factor$s);
        
        if(iter <= 0) next; # Burn-in
        
        all.alpha[,iter] = factor$alpha;
        all.beta[,iter]  = factor$beta;
        all.gamma[,iter] = factor$gamma;
        all.u[,,iter]    = factor$u;
        all.v[,,iter]    = factor$v;
        all.s[,,iter]    = factor$s;
        all.z_avg[,,iter] = z_avg;
        all.corpus_topic[,iter] = factor$corpus_topic;
        o = obs$y - factor$alpha[obs$user] - factor$beta[obs$item] -
            sum_margin(factor$u[obs$user,] * factor$v[obs$item,], 1) -
            sum_margin(factor$s[obs$user,] *    z_avg[obs$item,], 1);
        all.o_gamma[,iter] = o * factor$gamma[obs$user];
        all.o[,iter] = o;
        
        tc = getTopicCounts.R(corpus, factor$corpus_topic, length(factor$beta), ncol(factor$s), max(corpus$term));
        all.cnt_item_topic[,,iter] = tc$cnt_item_topic;
        all.cnt_topic_term[,,iter] = tc$cnt_topic_term;
        all.cnt_topic[,iter] = tc$cnt_topic;
        all.cnt_item[,iter]  = tc$cnt_item;
    }
    
    out = list();
    out$mean = list(
        alpha = apply(all.alpha, 1, mean),
        beta  = apply(all.beta,  1, mean),
        gamma = apply(all.gamma, 1, mean),
        u = apply(all.u, c(1,2), mean),
        v = apply(all.v, c(1,2), mean),
        s = apply(all.s, c(1,2), mean),
        z_avg = apply(all.z_avg, c(1,2), mean),
        gamma2 = apply(all.gamma^2, 1, mean),
        o_gamma = apply(all.o_gamma, 1, mean),
        corpus_topic = factor$corpus_topic
    );
    out$sumvar = list(
        alpha = sum(apply(all.alpha, 1, var)),
        beta  = sum(apply(all.beta,  1, var)),
        gamma = sum(apply(all.gamma, 1, var)),
        u = sum(apply(all.u, c(1,2), var)),
        v = sum(apply(all.v, c(1,2), var)),
        s = sum(apply(all.s, c(1,2), var)),
        o_adj = sum(apply(all.o^2, 1, mean) - out$mean$o_gamma^2 / out$mean$gamma2[obs$user])
    );
    out$objval = list(
        eta = compute_LDAobj(try$eta, size$nTopics, size$nTerms, all.cnt_topic, all.cnt_topic_term),
        lambda = compute_LDAobj(try$lambda, size$nItems, size$nTopics, all.cnt_item, all.cnt_item_topic)
    );
    
    return(out);
}

compute_LDAobj <- function(etas, K, W, Z_k, Z_kl){
    obj = rep(NA, length(etas));
    for(i in 1:length(etas)){
        eta = etas[i];
        ans = apply(lgamma(Z_kl + eta), c(1,2), mean);
        ans = apply(lgamma(Z_k + W*eta), 1, mean) - apply(ans, 1, sum);
        obj[i] = K * (W * lgamma(eta) - lgamma(W*eta)) + sum(ans);
    }
    return(obj);
}

###
### Functions to draw samples
###

draw_alpha <- function(factor, obs, param, xb, g0x_user, z_avg, size){
    user = obs$user;
    item = obs$item;
    o = obs$y - xb * factor$gamma[user] - factor$beta[item];
    if(size$nFactors > 0) o = o - sum_margin(factor$u[user,] * factor$v[item,], 1);
    if(size$nTopics  > 0) o = o - sum_margin(factor$s[user,] *    z_avg[item,], 1);

    # print(z_avg);
    # print(o);
    
    mean = g0x_user;
    var  = param$var_alpha;
    
    # compute var
    temp = aggregate(1/param$var_y, list(user), sum);
    userIDs = temp[,1];
    var[userIDs] = 1 / (1/param$var_alpha[userIDs] + temp[,2]);
    # compute mean
    temp = aggregate(o/param$var_y, list(user), sum);
    userIDs = temp[,1];
    mean[userIDs] = var[userIDs] * (g0x_user[userIDs]/param$var_alpha[userIDs] + temp[,2]);

    # draw sample
    out = list(
        sample = rnorm(size$nUsers, mean=mean, sd=sqrt(var)),
        mean = mean,
        var  = var
    );
    return(out);
}

draw_beta <- function(factor, obs, param, xb, d0x_item, z_avg, size){
    user = obs$user;
    item = obs$item;
    o = obs$y - xb * factor$gamma[user] - factor$alpha[user];
    if(size$nFactors > 0) o = o - sum_margin(factor$u[user,] * factor$v[item,], 1);
    if(size$nTopics  > 0) o = o - sum_margin(factor$s[user,] *    z_avg[item,], 1);

    mean = d0x_item;
    var  = param$var_beta;
    
    # compute var
    temp = aggregate(1/param$var_y, list(item), sum);
    itemIDs = temp[,1];
    var[itemIDs] = 1 / (1/param$var_beta[itemIDs] + temp[,2]);
    # compute mean
    temp = aggregate(o/param$var_y, list(item), sum);
    itemIDs = temp[,1];
    mean[itemIDs] = var[itemIDs] * (d0x_item[itemIDs]/param$var_beta[itemIDs] + temp[,2]);
    # draw sample
    out = list(
        sample = rnorm(size$nItems, mean=mean, sd=sqrt(var)),
        mean = mean,
        var  = var
    );
    return(out);
}

draw_gamma <- function(factor, obs, param, xb, c0x_user, z_avg, size){
    user = obs$user;
    item = obs$item;
    o = obs$y - factor$alpha[user] - factor$beta[item];
    if(size$nFactors > 0) o = o - sum_margin(factor$u[user,] * factor$v[item,], 1);
    if(size$nTopics  > 0) o = o - sum_margin(factor$s[user,] *    z_avg[item,], 1);

    mean = c0x_user;
    var  = param$var_gamma;
    
    # compute var
    temp = aggregate((xb^2)/param$var_y, list(user), sum);
    userIDs = temp[,1];
    var[userIDs] = 1 / (1/param$var_gamma[userIDs] + temp[,2]);
    # compute mean
    temp = aggregate((o*xb)/param$var_y, list(user), sum);
    userIDs = temp[,1];
    mean[userIDs] = var[userIDs] * (c0x_user[userIDs]/param$var_gamma[userIDs] + temp[,2]);
    # draw sample
    out = list(
        sample = rnorm(size$nUsers, mean=mean, sd=sqrt(var)),
        mean = mean,
        var  = var
    );
    return(out);
}

draw_u <- function(factor, obs, param, xb, Gx_user, z_avg, size){
    user = obs$user;
    item = obs$item;
    o = obs$y - xb * factor$gamma[user] - factor$alpha[user] - factor$beta[item];
    if(size$nTopics  > 0) o = o - sum_margin(factor$s[user,] *    z_avg[item,], 1);

    out = draw_vectors(o, user, item, Gx_user, factor$v, param$var_u, param$var_y);
    return(out);
}

draw_v <- function(factor, obs, param, xb, Dx_item, z_avg, size){
    user = obs$user;
    item = obs$item;
    o = obs$y - xb * factor$gamma[user] - factor$alpha[user] - factor$beta[item];
    if(size$nTopics  > 0) o = o - sum_margin(factor$s[user,] *    z_avg[item,], 1);
    
    out = draw_vectors(o, item, user, Dx_item, factor$u, param$var_v, param$var_y);
    return(out);
}

draw_s <- function(factor, obs, param, xb, Hx_user, z_avg, size){
    user = obs$user;
    item = obs$item;
    o = obs$y - xb * factor$gamma[user] - factor$alpha[user] - factor$beta[item];
    if(size$nFactors > 0) o = o - sum_margin(factor$u[user,] * factor$v[item,], 1);

    out = draw_vectors(o, user, item, Hx_user, z_avg, param$var_s, param$var_y);
    return(out);
}

draw_topic <- function(factor, obs, param, xb, corpus, size){
    user = obs$user;
    item = obs$item;

    topicCounts = getTopicCounts.R(corpus, factor$corpus_topic, length(factor$beta), ncol(factor$s), max(corpus$term));

    rest = obs$y - xb * factor$gamma[user] - factor$alpha[user] - factor$beta[item];
    if(size$nFactors > 0) rest = rest - sum_margin(factor$u[user,] * factor$v[item,], 1);

    if(is.matrix(factor$corpus_topic)) func = condProbSample_topic2.R
    else                               func = condProbSample_topic.R;
    
    ans = func(
        option=3,
        corpus=corpus, corpus_topic=factor$corpus_topic, 
        cnt_item_topic=topicCounts$cnt_item_topic, cnt_topic_term=topicCounts$cnt_topic_term, cnt_topic=topicCounts$cnt_topic,
        userIndex=user, itemIndex=item, rest=rest,
        s=factor$s, var_y=param$var_y, eta=param$eta, lambda=param$lambda, debug=1000
    );
    
    return(ans);
}

draw_vectors <- function(o, thisIndex, otherIndex, thisPredEff, otherEff, thisVar, var_y){
    mean    = array(NA, dim=dim(thisPredEff));
    var     = array(NA, dim=dim(thisVar));
    var.inv = array(NA, dim=dim(thisVar));
    for(i in 1:nrow(mean)){
        selectedObs = thisIndex == i;
        if(any(selectedObs)){
            selectedOther = otherIndex[selectedObs];
            o.selected = o[selectedObs];
            v.selected = otherEff[selectedOther,,drop=FALSE];

            var_y.selected = var_y[selectedObs];
            sum.vv = t(v.selected) %*% diag(1/var_y.selected, nrow=length(var_y.selected)) %*% v.selected;
            sum.ov = t(v.selected) %*% diag(1/var_y.selected, nrow=length(var_y.selected)) %*% o.selected;

            inv_var_ui = sym_inv.cholesky(thisVar[i,,]);
            
            var.inv[i,,] = sum.vv + inv_var_ui;
            eS <- sym_eigen(var.inv[i,,])
            
            var[i,,] = eS$vectors %*% diag(1/eS$values, nrow(eS$vectors)) %*% t(eS$vectors);
            mean[i,] = var[i,,] %*% (sum.ov + inv_var_ui %*% thisPredEff[i,]);
        }else{
            var[i,,] = thisVar[i,,];
            var.inv[i,,] = solve(var[i,,]);
            mean[i,] = thisPredEff[i,];
        }
    }
    
    out = list(
        sample = getMVNSample(mean, var.inv=var.inv),
        mean = mean,
        var  = var
    );
    return(out);
}


###
### Implemented to validate condProbSample_topic.C in c_funcs.R.
###
condProbSample_topic.R <- function(
    option, # 1:Sample, 2:Probabilities, 3:Sample&Probabilities
    corpus, corpus_topic, cnt_item_topic, cnt_topic_term, cnt_topic,
    userIndex, itemIndex, rest, # rest = the o in the paper
    s, var_y, eta, lambda, debug=0
){
    nItems = nrow(cnt_item_topic);
    nTopics= ncol(cnt_item_topic);
    nTerms = ncol(cnt_topic_term);
    cpsSize = nrow(corpus);
    nObs = length(userIndex);
    
    obsIndices = tapply(1:nObs, itemIndex, c, simplify=F);
    cpsIndices = tapply(1:cpsSize, corpus$item, c, simplify=F);
    if(any(names(cpsIndices) != 1:nItems)) stop("names(cpsIndices) != 1:nItems");
    
    probDist = matrix(NA, nrow=cpsSize, ncol=nTopics);
    
    for(j in 1:nItems){
        thisItem = j;
        cpsIndices_thisItem  = cpsIndices[[j]];
        
        this_obsIndex = obsIndices[[as.character(j)]];
        rest_thisItem = NULL;
        
        if(!is.null(this_obsIndex)){
            rest_thisItem = rest[this_obsIndex];
            var_y_thisItem= var_y[this_obsIndex];
            s_thisItem    = s[userIndex[this_obsIndex],];
        }
        
        for(n in 1:length(cpsIndices_thisItem)){
            cpsIndex   = cpsIndices_thisItem[n];
            thisTerm   = corpus[cpsIndex, "term"];
            thisWeight = corpus[cpsIndex, "weight"];
            thisTopic  = corpus_topic[cpsIndex];

            Z_kl = cnt_topic_term[,thisTerm];   Z_kl[thisTopic] = Z_kl[thisTopic] - thisWeight;
            Z_k  = cnt_topic;                   Z_k[ thisTopic] = Z_k[ thisTopic] - thisWeight;
            Z_jk = cnt_item_topic[thisItem,];   Z_jk[thisTopic] = Z_jk[thisTopic] - thisWeight;
            
            LDApart = ((Z_kl + eta) / (Z_k + nTerms * eta)) * (Z_jk + lambda);
            
            LL = rep(0, nTopics);
            if(length(rest_thisItem) > 0){
                for(k in 1:nTopics){
                    z_j = Z_jk;  z_j[k] = z_j[k] + thisWeight;
                    z_j = z_j / sum(z_j);
                    diff = rest_thisItem - s_thisItem %*% z_j;
                    LL[k] = sum(dnorm(diff, mean=0, sd=sqrt(var_y_thisItem), log=TRUE));
                }
            }
            # cat("(",j,",",n,") LDApart = [",LDApart,"]  LL = [",LL,"]\n");
            
            prob = log(LDApart) + LL;
            prob = exp(prob - max(prob));
            prob = prob / sum(prob);
            draw = rmultinom(1, 1, prob);
            
            pickedTopic = (1:nTopics)[draw == 1];
            corpus_topic[cpsIndex] = pickedTopic;
            cnt_topic_term[thisTopic,thisTerm]   = cnt_topic_term[thisTopic,thisTerm]   - thisWeight;
            cnt_topic_term[pickedTopic,thisTerm] = cnt_topic_term[pickedTopic,thisTerm] + thisWeight;
            cnt_item_topic[thisItem,thisTopic]   = cnt_item_topic[thisItem,thisTopic]   - thisWeight;
            cnt_item_topic[thisItem,pickedTopic] = cnt_item_topic[thisItem,pickedTopic] + thisWeight;
            cnt_topic[thisTopic]   = cnt_topic[thisTopic]   - thisWeight;
            cnt_topic[pickedTopic] = cnt_topic[pickedTopic] + thisWeight;
            
            probDist[cpsIndex,] = prob;
        }
    }
    out = list(corpus_topic=corpus_topic, cnt_item_topic=cnt_item_topic, cnt_topic_term=cnt_topic_term, cnt_topic=cnt_topic, probDist=probDist);
    return(out);
}


###
### Implemented to validate condProbSample_topic.C in c_funcs.R.
###
condProbSample_topic2.R <- function(
    option, # 1:Sample, 2:Probabilities, 3:Sample&Probabilities
    corpus, corpus_topic, cnt_item_topic, cnt_topic_term, cnt_topic,
    userIndex, itemIndex, rest, # rest = the o in the paper
    s, var_y, eta, lambda, debug=0
){
    nItems = nrow(cnt_item_topic);
    nTopics= ncol(cnt_item_topic);
    nTerms = ncol(cnt_topic_term);
    cpsSize = nrow(corpus);
    nObs = length(userIndex);
    
    if(!is.matrix(corpus_topic)) stop("!is.matrix(corpus_topic)");
    if(any(dim(corpus_topic) != c(cpsSize, nTopics))) stop("dim(corpus_topic) != c(cpsSize, nTopics)");

    obsIndices = tapply(1:nObs, itemIndex, c, simplify=F);
    cpsIndices = tapply(1:cpsSize, corpus$item, c, simplify=F);
    if(any(names(cpsIndices) != 1:nItems)) stop("names(cpsIndices) != 1:nItems");
    
    probDist = matrix(NA, nrow=cpsSize, ncol=nTopics);
    
    for(j in 1:nItems){
        thisItem = j;
        cpsIndices_thisItem  = cpsIndices[[j]];
        
        this_obsIndex = obsIndices[[as.character(j)]];
        rest_thisItem = NULL;
        
        if(!is.null(this_obsIndex)){
            rest_thisItem = rest[this_obsIndex];
            var_y_thisItem= var_y[this_obsIndex];
            s_thisItem    = s[userIndex[this_obsIndex],];
        }
        
        for(n in 1:length(cpsIndices_thisItem)){
            cpsIndex   = cpsIndices_thisItem[n];
            thisTerm   = corpus[cpsIndex, "term"];
            thisWeight = corpus[cpsIndex, "weight"];
            thisTopic  = corpus_topic[cpsIndex,];

            Z_kl = cnt_topic_term[,thisTerm];   Z_kl = Z_kl - thisWeight*thisTopic;
            Z_k  = cnt_topic;                   Z_k  = Z_k  - thisWeight*thisTopic;
            Z_jk = cnt_item_topic[thisItem,];   Z_jk = Z_jk - thisWeight*thisTopic;
            
            LDApart = ((Z_kl + eta) / (Z_k + nTerms * eta)) * (Z_jk + lambda);
            
            LL = rep(0, nTopics);
            if(length(rest_thisItem) > 0){
                for(k in 1:nTopics){
                    z_j = Z_jk;  z_j[k] = z_j[k] + thisWeight;
                    z_j = z_j / sum(z_j);
                    diff = rest_thisItem - s_thisItem %*% z_j;
                    LL[k] = sum(dnorm(diff, mean=0, sd=sqrt(var_y_thisItem), log=TRUE));
                }
            }
            # cat("(",j,",",n,") LDApart = [",LDApart,"]  LL = [",LL,"]\n");
            
            prob = log(LDApart) + LL;
            prob = exp(prob - max(prob));
            prob = prob / sum(prob);

            corpus_topic[cpsIndex,] = prob;
            cnt_topic_term[,thisTerm] = cnt_topic_term[,thisTerm] + thisWeight * (prob - thisTopic);
            cnt_item_topic[thisItem,] = cnt_item_topic[thisItem,] + thisWeight * (prob - thisTopic);
            cnt_topic                 = cnt_topic                 + thisWeight * (prob - thisTopic);
            
            probDist[cpsIndex,] = prob;
        }
    }
    out = list(corpus_topic=corpus_topic, cnt_item_topic=cnt_item_topic, cnt_topic_term=cnt_topic_term, cnt_topic=cnt_topic, probDist=probDist);
    return(out);
}


draw_alpha.C <- function(factor, obs, param, xb, g0x_user, z_avg, size, debug=0){
    user = obs$user;
    item = obs$item;
    o = obs$y - xb * factor$gamma[user] - factor$beta[item];
    if(size$nFactors > 0) o = o - sum_margin(factor$u[user,] * factor$v[item,], 1);
    if(size$nTopics  > 0) o = o - sum_margin(factor$s[user,] *    z_avg[item,], 1);

    ans = condMeanVarSample_singleDim.C(
        option=3, thisEffIndex=user, rest=o,
        fittedEff=g0x_user, multiplier=rep(1.0,length(o)), var_y=param$var_y, var_eff=param$var_alpha, debug=debug
    );
    return(ans);
}

draw_beta.C <- function(factor, obs, param, xb, d0x_item, z_avg, size, debug=0){
    user = obs$user;
    item = obs$item;
    o = obs$y - xb * factor$gamma[user] - factor$alpha[user];
    if(size$nFactors > 0) o = o - sum_margin(factor$u[user,] * factor$v[item,], 1);
    if(size$nTopics  > 0) o = o - sum_margin(factor$s[user,] *    z_avg[item,], 1);

    ans = condMeanVarSample_singleDim.C(
        option=3, thisEffIndex=item, rest=o,
        fittedEff=d0x_item, multiplier=rep(1.0,length(o)), var_y=param$var_y, var_eff=param$var_beta, debug=debug
    );
    return(ans);
}

draw_gamma.C <- function(factor, obs, param, xb, c0x_user, z_avg, size, debug=0){
    user = obs$user;
    item = obs$item;
    o = obs$y - factor$alpha[user] - factor$beta[item];
    if(size$nFactors > 0) o = o - sum_margin(factor$u[user,] * factor$v[item,], 1);
    if(size$nTopics  > 0) o = o - sum_margin(factor$s[user,] *    z_avg[item,], 1);
        
    ans = condMeanVarSample_singleDim.C(
        option=3, thisEffIndex=user, rest=o,
        fittedEff=c0x_user, multiplier=xb, var_y=param$var_y, var_eff=param$var_gamma, debug=debug
    );
    return(ans);
}

draw_u.C <- function(factor, obs, param, xb, Gx_user, z_avg, size, debug=0){
    user = obs$user;
    item = obs$item;
    o = obs$y - xb * factor$gamma[user] - factor$alpha[user] - factor$beta[item];
    if(size$nTopics  > 0) o = o - sum_margin(factor$s[user,] *    z_avg[item,], 1);
        
    ans = condMeanVarSample_multiDim.C(
        option=3, thisEffIndex=user, otherEffIndex=item, rest=o,
        fittedEff=Gx_user, otherEff=factor$v, var_y=param$var_y, var_eff=param$var_u, oi=NULL, debug=debug
    );
    return(ans);
}

draw_v.C <- function(factor, obs, param, xb, Dx_item, z_avg, size, debug=0){
    user = obs$user;
    item = obs$item;
    o = obs$y - xb * factor$gamma[user] - factor$alpha[user] - factor$beta[item];
    if(size$nTopics  > 0) o = o - sum_margin(factor$s[user,] *    z_avg[item,], 1);

    ans = condMeanVarSample_multiDim.C(
        option=3, thisEffIndex=item, otherEffIndex=user, rest=o,
        fittedEff=Dx_item, otherEff=factor$u, var_y=param$var_y, var_eff=param$var_v, oi=NULL, debug=debug
    );
    return(ans);
}

draw_s.C <- function(factor, obs, param, xb, Hx_user, z_avg, size, debug=0){
    user = obs$user;
    item = obs$item;
    o = obs$y - xb * factor$gamma[user] - factor$alpha[user] - factor$beta[item];
    if(size$nFactors > 0) o = o - sum_margin(factor$u[user,] * factor$v[item,], 1);

    ans = condMeanVarSample_multiDim.C(
        option=3, thisEffIndex=user, otherEffIndex=item, rest=o,
        fittedEff=Hx_user, otherEff=z_avg, var_y=param$var_y, var_eff=param$var_s, oi=NULL, debug=debug
    );
    return(ans);
}

draw_topic.C <- function(factor, obs, param, xb, corpus, size, debug=0){
    user = obs$user;
    item = obs$item;

    topicCounts = getTopicCounts(corpus, factor$corpus_topic, length(factor$beta), ncol(factor$s), max(corpus$term));

    rest = obs$y - xb * factor$gamma[user] - factor$alpha[user] - factor$beta[item];
    if(size$nFactors > 0) rest = rest - sum_margin(factor$u[user,] * factor$v[item,], 1);

    ans = condProbSample_topic.C(
        option=3,
        corpus=corpus, corpus_topic=factor$corpus_topic, 
        cnt_item_topic=topicCounts$cnt_item_topic, cnt_topic_term=topicCounts$cnt_topic_term, cnt_topic=topicCounts$cnt_topic,
        userIndex=user, itemIndex=item, rest=rest,
        s=factor$s, var_y=param$var_y, eta=param$eta, lambda=param$lambda, debug=debug
    );
    return(ans);
}


