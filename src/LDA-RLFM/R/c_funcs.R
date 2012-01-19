### Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### Inversion of a symmetric metrix using Cholesky
###     This function is implemented to make sure the C implementation generates exactly
###     the same output as the R implementation
###
sym_inv.cholesky <- function(x){
    if(!is.double(x)) stop("x is not a double matrix");
    if(!is.matrix(x)) stop("x is not a double matrix");
    if(nrow(x) != ncol(x)) stop("x is not a square matrix");
    n = as.integer(nrow(x));
    ans = .C("sym_inv_byCholesky", x, n, as.integer(1), DUP=TRUE);
    return(ans[[1]]);
}

###
### Eigen decomposition of a symmetric metrix
###     This function is implemented to make sure the C implementation generates exactly
###     the same output as the R implementation
###
sym_eigen <- function(x){
    if(!is.matrix(x)) stop("x is not matrix");
    if(nrow(x) != ncol(x)) stop("nrow(x) != ncol(x)");
    if(!is.double(x)) stop("x is not a double-precision matrix");
    output = list(values  = rep(0.0, nrow(x)),
                  vectors = matrix(0.0, nrow=nrow(x), ncol=ncol(x)));
    .C("sym_eigen", x, nrow(x), output$values, output$vectors, DUP=FALSE);
    return(output);
}

###
### Sum up each row (or column) of a matrix
###     side=1: Sum up each row and return a vector with length nrow
###     side=2: Sum up each column and return a vector with length ncol
###
sum_margin <- function(A, side){
    if(!is.matrix(A)) stop("'A' is not a matrix");
    if(!is.double(A)) stop("'A' is not double");
    if(side == 1) out = double(nrow(A))
    else if(side == 2) out = double(ncol(A))
    else stop("Unknown side=",side," (side=1 for rows, side=2 for columns)");
    ans = .C("sum_margin", out, A, as.integer(nrow(A)), as.integer(ncol(A)), as.integer(side), DUP=FALSE);
    return(out);
}

###
### out[k,j] = sum_{x s.t. groupBy[x] = j} mat[k,select[x]] * weight[x]
###
selectColumn_agg_sum <- function(mat, groupBy, select, weight=NULL){
    out = matrix(0.0, nrow=nrow(mat), ncol=max(groupBy));
    if(!is.double(mat)) stop("!is.double(mat)");
    if(!is.integer(groupBy)) stop("!is.integer(groupBy)");
    if(!is.integer(select)) stop("!is.integer(select)");
    if(length(select) != length(groupBy)) stop("length(select) != length(groupBy)");
    nWeightw = 0;
    if(!is.null(weight)){
        nWeights = length(weight);
        if(nWeights != length(select)) stop("length(weight) != length(select)");
        if(!is.double(weight)) stop("!is.double(weight)");
    }
    ans = .C("selectColumn_agg_sum",
        out, as.integer(nrow(out)), as.integer(ncol(out)),
        mat, as.integer(nrow(mat)), as.integer(ncol(mat)),
        select, groupBy, as.integer(length(select)),
        weight, as.integer(nWeightw),
        DUP=FALSE
    );
    return(out);
}

###
### margin = 1: each row sum up to one
### margin = 2: each column sum up to one
###
normalize_sumToOne2D <- function(mat, margin){
    out = matrix(0.0, nrow=nrow(mat), ncol=ncol(mat));
    if(!is.double(mat)) stop("!is.double(mat)");
    ans = .C("normalize_sumToOne2D",
        out, mat, as.integer(nrow(mat)),as.integer(ncol(mat)), as.integer(margin),
        DUP=FALSE
    );
    return(out);
}


###
### See condMeanVarSample_singleDim in MCEM_EStep.c
###
condMeanVarSample_singleDim.C <- function(
    option, # 1:Sample, 2:Mean&Var, 3:Sample&Mean&Var
    thisEffIndex, rest, # rest = the o in the paper
    fittedEff, multiplier, var_y, var_eff, debug=0
){
    nObs = as.integer(length(thisEffIndex));
    nThisEff  = as.integer(length(fittedEff));
    nVar_y    = as.integer(length(var_y));
    nVar_eff  = as.integer(length(var_eff));
    
    if(length(rest) != nObs) stop("length(rest) != nObs");
    if(!is.null(multiplier) && length(multiplier) != nObs) stop("length(multiplier) != nObs");
    if(!(nVar_y == 1 || nVar_y == nObs)) stop("length(var_y) has problem");
    if(!(nVar_eff == 1 || nVar_eff == nThisEff)) stop("length(var_eff) has problem");
    
    out = list(sample=as.double(NULL), mean=as.double(NULL), var=as.double(NULL));
    if(option == 1 || option == 3)  out$sample = rep(0.0, nThisEff);
    if(option == 2 || option == 3){ out$mean   = rep(0.0, nThisEff);   out$var = rep(0.0, nThisEff);}

    if(!is.double(rest)) stop("!is.double(rest)");
    if(!is.null(multiplier) && !is.double(multiplier)) stop("!is.double(multiplier)");
    if(!is.double(fittedEff)) stop("!is.double(fittedEff)");
    if(!is.double(var_y)) stop("!is.double(var_y)");
    if(!is.double(var_eff)) stop("!is.double(var_eff)");
    
    ans = .C("condMeanVarSample_singleDim",
        # OUTPUT
        out$sample, out$mean, out$var,
        # INPUT
        as.integer(option), as.integer(thisEffIndex), rest,
        fittedEff, multiplier, var_y, var_eff,
        nObs, nThisEff, nVar_y, nVar_eff,
        as.integer(debug), DUP=FALSE
    );

    return(out);
}

###
### See generateObsIndex in util.c
###
generateObsIndex <- function(
    effIndex, nEff, debug
){
    nEff = as.integer(nEff);
    nObs   = as.integer(length(effIndex));
    
    out = list(
        obsIndex = integer(nObs),
        start    = integer(nEff),
        num      = integer(nEff)
    );
    
    if(!is.integer(effIndex)) stop("!is.integer(effIndex)");
    
    ans = .C("generateObsIndex",
        out$obsIndex, out$start, out$num,
        # INPUT
        effIndex, nObs, nEff,
        # OTHER
        as.integer(debug),
        DUP=FALSE
    );
    return(out);
}

###
### See condMeanVarSample_multiDim in MCEM_EStep.c
###
condMeanVarSample_multiDim.C <- function(
    option, # 1:Sample, 2:Mean&Var, 3:Sample&Mean&Var
    thisEffIndex, otherEffIndex, rest, # rest = the o in the paper
    fittedEff, otherEff, var_y, var_eff, oi=NULL, debug=0
){
    nObs = as.integer(length(thisEffIndex));
    nThisEff  = as.integer(nrow(fittedEff));
    nOtherEff = as.integer(nrow(otherEff));
    nFactors  = as.integer(ncol(fittedEff));
    
    nVar_y = as.integer(length(var_y));
    nVar_eff = as.integer(length(var_eff));
    if(nVar_eff > 1){
        nVar_eff = as.integer(dim(var_eff)[1]);
    }
    
    if(ncol(otherEff) != nFactors) stop("ncol(otherEff) != nFactors");
    if(length(otherEffIndex) != nObs) stop("length(otherEffIndex) != nObs");
    if(length(rest) != nObs) stop("length(rest) != nObs");
    
    if(!(nVar_y == 1 || nVar_y == nObs)) stop("length(var_y) has problem");
    if(!(nVar_eff == 1 || nVar_eff == nThisEff)) stop("length(var_eff) has problem");
    
    out = list(sample=as.double(NULL), mean=as.double(NULL), var=as.double(NULL));
    if(option == 1 || option == 3)  out$sample = matrix(0.0, nrow=nThisEff, ncol=nFactors);
    if(option == 2 || option == 3){ out$mean   = matrix(0.0, nrow=nThisEff, ncol=nFactors);  out$var = array(0.0, dim=c(nThisEff, nFactors, nFactors));}

    if(!is.double(rest)) stop("!is.double(rest)");
    if(!is.double(fittedEff)) stop("!is.double(fittedEff)");
    if(!is.double(otherEff)) stop("!is.double(otherEff)");
    if(!is.double(var_y)) stop("!is.double(var_y)");
    if(!is.double(var_eff)) stop("!is.double(var_eff)");
    
    if(is.null(oi)){
        cat("obsIndex is null; built it now!\n");
        oi = generateObsIndex(thisEffIndex, nThisEff, debug);
    }
    
    if(length(oi$obsIndex) != nObs) stop("length(oi$obsIndex) != nObs");
    if(length(oi$start) != nThisEff) stop("length(oi$start) != nThisEff");
    if(length(oi$num) != nThisEff) stop("length(oi$num) != nThisEff");
    if(!is.integer(oi$obsIndex)) stop("!is.integer(oi$obsIndex)");
    if(!is.integer(oi$start)) stop("!is.integer(oi$start)");
    if(!is.integer(oi$num)) stop("!is.integer(oi$num)");

    ans = .C("condMeanVarSample_multiDim",
        # OUTPUT
        out$sample, out$mean, out$var,
        # INPUT
        as.integer(option), as.integer(thisEffIndex), as.integer(otherEffIndex), rest,
        fittedEff, otherEff, var_y, var_eff,
        nObs, nThisEff, nOtherEff, nFactors, nVar_y, nVar_eff,
        oi$obsIndex, oi$start, oi$num,
        as.integer(debug),
        DUP=FALSE
    );
    
    return(out);
}


###
### See condProbSample_topic in MCEM_EStep.c
###
###     corpus = data.frame(item, term, weight), where weight can be null
###
condProbSample_topic.C <- function(
    option, # 1:Sample, 2:Probabilities, 3:Sample&Probabilities
    corpus, corpus_topic, cnt_item_topic, cnt_topic_term, cnt_topic,
    userIndex, itemIndex, rest, # rest = the o in the paper
    s, var_y, eta, lambda, oi=NULL, ci=NULL, nRatingExponent=0,
    debug=0, verbose=0
){
    nObs    = as.integer(length(userIndex));
    nUsers  = as.integer(nrow(s));
    nItems  = as.integer(nrow(cnt_item_topic));
    nTopics = as.integer(ncol(cnt_item_topic));
    nTerms  = as.integer(ncol(cnt_topic_term));
    corpusSize = as.integer(nrow(corpus));
    
    nVar_y  = as.integer(length(var_y));
    nLambda = as.integer(length(lambda));

    if(is.matrix(corpus_topic) && any(dim(corpus_topic) != c(corpusSize, nTopics)))
        stop("!is.matrix(corpus_topic) || dim(corpus_topic) != c(corpusSize, nTopics)");
    if(!is.matrix(corpus_topic) && length(corpus_topic) != corpusSize) stop("length(corpus_topic) != corpusSize");
    if(any(dim(cnt_item_topic) != c(nItems, nTopics))) stop("dim(cnt_item_topic) != c(nItems, nTopics)");
    if(any(dim(cnt_topic_term) != c(nTopics, nTerms))) stop("dim(cnt_topic_term) != c(nTopics, nTerms)");
    if(length(cnt_topic) != nTopics) stop("length(cnt_topic) != nTopics");
    if(length(itemIndex) != nObs) stop("length(itemIndex) != nObs");
    if(length(rest) != nObs) stop("length(rest) != nObs");
    if(any(dim(s) != c(nUsers, nTopics))) stop("dim(s) != c(nUsers, nTopics)");
    if(!(nVar_y == 1 || nVar_y == nObs)) stop("length(var_y) has problem");
    if(!(nLambda == 1 || nLambda == nTopics)) stop("length(lambda) has problem");
    if(length(eta) != 1) stop("length(eta) != 1");
    if(length(nRatingExponent) != 1) stop("length(nRatingExponent) != 1");
    
    out = list(corpus_topic=corpus_topic, cnt_item_topic=cnt_item_topic, cnt_topic_term=cnt_topic_term, cnt_topic=cnt_topic, probDist=as.double(NULL));
    if(option == 2 || option == 3){ out$probDist = matrix(0.0, nrow=corpusSize, ncol=nTopics); }

    if(is.matrix(corpus_topic)  && !is.double(out$corpus_topic)) stop("!is.double(out$corpus_topic)");
    if(!is.matrix(corpus_topic) && !is.integer(out$corpus_topic)) stop("!is.integer(out$corpus_topic)");
    if(!is.double(out$cnt_item_topic)) stop("!is.double(out$cnt_item_topic)");
    if(!is.double(out$cnt_topic_term)) stop("!is.double(out$cnt_topic_term)");
    if(!is.double(out$cnt_topic)) stop("!is.double(out$cnt_topic)");
    if(!is.double(rest)) stop("!is.double(rest)");
    if(!is.double(s)) stop("!is.double(s)");
    if(!is.double(var_y)) stop("!is.double(var_y)");
    if(!is.double(eta)) stop("!is.double(eta)");
    if(!is.double(lambda)) stop("!is.double(lambda)");
    if(!is.integer(userIndex)) stop("!is.integer(userIndex)");
    if(!is.integer(itemIndex)) stop("!is.integer(itemIndex)");
    if(!is.integer(corpus$item)) stop("!is.integer(corpus$item)");
    if(!is.integer(corpus$term)) stop("!is.integer(corpus$term)");
    if(!is.null(corpus$weight) && !is.double(corpus$weight)) stop("!is.double(corpus$weight)");
    
    if(is.null(oi)){
        cat("obsIndex is null; built it now!\n");
        oi = generateObsIndex(itemIndex, nItems, debug);
    }
    if(is.null(ci)){
        cat("cpsIndex is null; built it now!\n");
        ci = generateObsIndex(corpus$item, nItems, debug);
    }

    if(!is.integer(oi$obsIndex)) stop("!is.integer(oi$obsIndex)");
    if(!is.integer(oi$start)) stop("!is.integer(oi$start)");
    if(!is.integer(oi$num)) stop("!is.integer(oi$num)");
    if(!is.integer(ci$obsIndex)) stop("!is.integer(ci$obsIndex)");
    if(!is.integer(ci$start)) stop("!is.integer(ci$start)");
    if(!is.integer(ci$num)) stop("!is.integer(ci$num)");
    
    if(length(oi$obsIndex) != nObs) stop("length(oi$obsIndex) != nObs");
    if(length(oi$start) != nItems) stop("length(oi$start) != nItems");
    if(length(oi$num) != nItems) stop("length(oi$num) != nItems");
    if(length(ci$obsIndex) != corpusSize) stop("length(ci$obsIndex) != corpusSize");
    if(length(ci$start) != nItems) stop("length(ci$start) != nItems");
    if(length(ci$num) != nItems) stop("length(ci$num) != nItems");

    if(is.matrix(corpus_topic)) fn.name = "condProbSample_topic2"
    else                        fn.name = "condProbSample_topic";
    
    ans = .C(fn.name,
        # INPUT & OUTPUT
        out$corpus_topic, out$cnt_item_topic, out$cnt_topic_term, out$cnt_topic,
        # OUTPUT
        out$probDist,
        # INPUT
        as.integer(option),
        rest, s, var_y, eta, lambda,
        userIndex, itemIndex, corpus$item, corpus$term, corpus$weight,
        nUsers, nItems, nObs, corpusSize, nTopics, nTerms, nVar_y, nLambda,
        oi$obsIndex, oi$start, oi$num, ci$obsIndex, ci$start, ci$num,
        as.double(nRatingExponent),
        # OTHER
        as.integer(debug), as.integer(verbose), as.integer(1),
        DUP=FALSE
    );

    return(out);
}

###
### MCEM_EStep.C. See MCEM_EStep(...) in MCEM_EStep.c
###
### IN&OUT: factor  = list(alpha, beta, gamma, u, v, s, corpus_topic);
###                    input: initial factor values
###                   output: sample means of alpha, beta, gamma, u, v, s
###                           last sample of corpus_topic
###                   THE VALUES in this list WILL BE CHANGED BY THIS FUNCTION (CALL BY REFERENCE!!)
###
### INPUT:  obs     = data.frame(y, user, item);
###         corpus  = data.frame(item, term, weight);
###         feature = list(x_dyad, x_user, x_item);
###         param   = list(b, g0, d0, c0, G, D, H, var_y, var_alpha, var_beta, var_gamma, var_u, var_v, var_s, lambda, eta);
###         try     = list(lambda, eta);
###         
###         NOTE: When var_FACTOR == NULL, no sample will be drawn for that FACTOR
###         
### OUTPUT: mean    = list(alpha, beta, gamma, u, v, s, z_avg, gamma2, o_gamma, corpus_topic, phi);  # gamma2 = gamma^2
###         sumvar  = list(alpha, beta, gamma, u, v, s, o_adj);
###         objval  = list(eta, lambda);
###         sampvar = list(alpha, beta, gamma, u, v, s, z_avg);
###
### OPTION: userFactorVar=0, itemFactorVar=0, userTopicVar=0, itemTopicVar=0, drawTopicSample=1 by default
###         Whether to output variance components
###
###   isOldUser[i]: [TRUE/FALSE] Whehter the ith user is an old user; default: all FALSE
###                 For old users, we set g0x_user[i] = alpha[i], c0x_user[i] = gamma[i], Gx_user[i,] = u[i,] and Hx_user[i,] = s[i,]
###   isOldItem[j]: [TRUE/FALSE] Whehter the jth item is an old item; default: all FALSE
###                 For old items, we set d0x_item[j] = beta[j], Dx_item[j] = v[j,]
###
###   If drawTopicSample == 0, no topic sample will be drawn, but use the input factor$z_avg.
###   If drawTopicSample == 1, draw topic samples z_avg[j,k] is avg_w 1{term w in item j belongs to topic k}
###   If drawTopicSample == 2, z_avg[j,k] = avg_w Pr(term w in item j belongs to topic k)
###                                         instead of drawing samples from the probabilities
###                            However, we still draw multinomial samples for corpus_topic
###   If drawTopicSample == 3, z_avg[j,k] =  avg_w Pr(term w in item j belongs to topic k)
###                            corpus_topic[w,k] = Pr(term w in corpus$item[w] belongs to topic k)
###                            No multinomial sampling at all!!
### NOTE: factor$alpha, factor$beta and factor$gamma cannot be NULL!!
###
MCEM_EStep.C <- function(
    factor, obs, corpus, feature, param, try, nSamples, nBurnIn=1,
    userFactorVar=0, itemFactorVar=0, userTopicVar=0, itemTopicVar=0, drawTopicSample=1,
    isOldUser=NULL, isOldItem=NULL,
    nRatingExponent=0,
    debug=0, verbose=0
){
    size = syncheck.LDA_RLFM.spec(factor=factor, obs=obs, corpus=corpus, feature=feature, param=param, 
                                  is.corpus_topic.matrix=(drawTopicSample == 3));
    
    if(size$nTopics > 0 && drawTopicSample == 0 && is.null(factor$z_avg)) stop("factor$z_avg should not be NULL");
    
    z_avg = NULL;
    if(!is.null(factor$z_avg)) z_avg = factor$z_avg;
    if(is.null(z_avg)) z_avg = matrix(as.double(0), nrow=size$nItems, ncol=size$nTopics);
    
    if(nrow(z_avg) != size$nItems || ncol(z_avg) != size$nTopics) stop("z_avg has a wrong dim");
    if(!is.double(z_avg)) stop("z_avg must be double");
    
    if(is.null(try)) try = list();
    if(is.null(try[["lambda"]])) try[["lambda"]] = double(0)
    else if(!is.double(try[["lambda"]])) stop("!is.double(try$lambda)");
    if(is.null(try[["eta"]])) try[["eta"]] = double(0)
    else if(!is.double(try[["eta"]])) stop("!is.double(try$eta)");
    
    gamma2  = double(size$nUsers);
    o_gamma = double(size$nObs);
    
    sumvar  = list(alpha=double(1), beta=double(1), gamma=double(1), u=double(1), v=double(1), s=double(1), o_adj=double(1));
    objval  = list(eta=double(length(try$eta)), lambda=double(length(try$lambda)));
    sampvar = list();
    if(userFactorVar != 0){
        sampvar$alpha = double(size$nUsers);
        sampvar$gamma = double(size$nUsers);
        sampvar$u     = array(double(1),dim=c(size$nUsers, size$nFactors, size$nFactors));
    }
    if(itemFactorVar != 0){
        sampvar$beta = double(size$nItems);
        sampvar$v    = array(double(1),dim=c(size$nItems, size$nFactors, size$nFactors));
    }
    if(userTopicVar != 0) sampvar$s = array(double(1),dim=c(size$nUsers, size$nTopics, size$nTopics));
    if(itemTopicVar != 0) sampvar$z_avg = array(double(1),dim=c(size$nItems, size$nTopics, size$nTopics));
    
    xb  = feature$x_dyad %*% param$b;
    g0x_user = NULL; c0x_user = NULL; Gx_user = NULL; Hx_user = NULL;
    if(is.null(isOldUser)){
        if(!is.null(param$g0)) g0x_user = feature$x_user %*% param$g0;
        if(!is.null(param$c0)) c0x_user = feature$x_user %*% param$c0;
        if(!is.null(param$G))  Gx_user  = feature$x_user %*% param$G;
        if(!is.null(param$H))  Hx_user  = feature$x_user %*% param$H;
    }else{
        if(length(isOldUser) != size$nUsers) stop("length(isOldUser) != nUsers");
        x_user.new = feature$x_user[!isOldUser,,drop=FALSE];
        if(!is.null(param$g0)){ g0x_user = factor$alpha;    g0x_user[!isOldUser] = x_user.new %*% param$g0;}
        if(!is.null(param$c0)){ c0x_user = factor$gamma;    c0x_user[!isOldUser] = x_user.new %*% param$c0;}
        if(!is.null(param$G)){  Gx_user  = factor$u;        Gx_user[!isOldUser,] = x_user.new %*% param$G;}
        if(!is.null(param$H)){  Hx_user  = factor$s;        Hx_user[!isOldUser,] = x_user.new %*% param$H;}
    }
    
    d0x_item = NULL;  Dx_item = NULL;
    if(is.null(isOldItem)){
        if(!is.null(param$d0)) d0x_item = feature$x_item %*% param$d0;
        if(!is.null(param$D))  Dx_item  = feature$x_item %*% param$D;
    }else{
        if(length(isOldItem) != size$nItems) stop("length(isOldItem) != nItems");
        x_item.new = feature$x_item[!isOldItem,,drop=FALSE];
        if(!is.null(param$d0)){ d0x_item = factor$beta;     d0x_item[!isOldItem] = x_item.new %*% d0;}
        if(!is.null(param$D)){  Dx_item  = factor$v;        Dx_item[!isOldItem,] = x_item.new %*% D;}
    }

    # Checks
    if(is.null(g0x_user) && size$nVar_alpha > 0) stop("error!");
    if(is.null(c0x_user) && size$nVar_gamma > 0) stop("error!");
    if(is.null(Gx_user)  && size$nVar_u > 0)     stop("error!");
    if(is.null(Hx_user)  && size$nVar_s > 0)     stop("error!");
    if(is.null(d0x_item) && size$nVar_beta > 0)  stop("error!");
    if(is.null(Dx_item)  && size$nVar_v > 0)     stop("error!");
    
    if(length(nRatingExponent) != 1) stop("length(nRatingExponent) != 1");
    
    problem.dim = c(size$nObs, size$corpusSize, size$nUsers, size$nItems, size$nTerms, 
                    size$nFactors, size$nTopics, size$nCorpusWeights,
                    size$nVar_y, size$nVar_alpha, size$nVar_beta, size$nVar_gamma, size$nVar_u, size$nVar_v, size$nVar_s,
                    length(try$eta), length(try$lambda));
    if(!is.integer(problem.dim)) stop("!is.integer(problem.dim)");
    
    phi = NULL;
    if(drawTopicSample == 0) phi = factor$phi;
    if(size$nTopics > 0 && drawTopicSample != 0){
        phi = matrix(0.0, nrow=size$nTopics, ncol=size$nTerms);
    }
    
    if(drawTopicSample %in% c(0,1,2)){
        ans = .C("MCEM_EStep",
            # INPUT (initial factor values) & OUTPUT (Monte Carlo mean of factor values)
            factor$alpha, factor$beta, factor$gamma, factor$u, factor$v, factor$s, 
            factor$corpus_topic, z_avg,
            # OUTPUT
            sumvar$alpha, sampvar$alpha, sumvar$beta, sampvar$beta, sumvar$gamma, sampvar$gamma, gamma2,
            sumvar$u,     sampvar$u,     sumvar$v,    sampvar$v,    sumvar$s,     sampvar$s,
            sampvar$z_avg, objval$eta, objval$lambda, o_gamma, sumvar$o_adj,
            phi,
            # INPUT
            as.integer(nSamples), as.integer(nBurnIn),
            obs$user, obs$item, corpus$item, corpus$term, corpus$weight, 
            obs$y, xb, g0x_user, d0x_item, c0x_user, Gx_user, Dx_item, Hx_user,
            param$var_y, param$var_alpha, param$var_beta, param$var_gamma, param$var_u, param$var_v, param$var_s,
            param$eta, param$lambda, try$eta, try$lambda, as.double(nRatingExponent),
            problem.dim, as.integer(length(problem.dim)),
            # OPTIONS
            as.integer(userFactorVar), as.integer(itemFactorVar), as.integer(userTopicVar), as.integer(itemTopicVar),
            as.integer(drawTopicSample), as.integer(debug), as.integer(verbose),
            DUP=FALSE
        );
    }else if(drawTopicSample == 3){
        ans = .C("MCEM_EStep2",
            # INPUT (initial factor values) & OUTPUT (Monte Carlo mean of factor values)
            factor$alpha, factor$beta, factor$gamma, factor$u, factor$v, factor$s, 
            factor$corpus_topic, z_avg,
            # OUTPUT
            sumvar$alpha, sampvar$alpha, sumvar$beta, sampvar$beta, sumvar$gamma, sampvar$gamma, gamma2,
            sumvar$u,     sampvar$u,     sumvar$v,    sampvar$v,    sumvar$s,     sampvar$s,
            sampvar$z_avg, objval$eta, objval$lambda, o_gamma, sumvar$o_adj,
            phi,
            # INPUT
            as.integer(nSamples), as.integer(nBurnIn),
            obs$user, obs$item, corpus$item, corpus$term, corpus$weight, 
            obs$y, xb, g0x_user, d0x_item, c0x_user, Gx_user, Dx_item, Hx_user,
            param$var_y, param$var_alpha, param$var_beta, param$var_gamma, param$var_u, param$var_v, param$var_s,
            param$eta, param$lambda, try$eta, try$lambda, as.double(nRatingExponent),
            problem.dim, as.integer(length(problem.dim)),
            # OPTIONS
            as.integer(userFactorVar), as.integer(itemFactorVar), as.integer(userTopicVar), as.integer(itemTopicVar),
            as.integer(drawTopicSample), as.integer(debug), as.integer(verbose),
            DUP=FALSE
        );
    }else stop("unknown drawTopicSample: ",drawTopicSample);
    
    output = list(
        mean=list(alpha=factor$alpha, beta=factor$beta, gamma=factor$gamma, 
                  u=factor$u, v=factor$v, s=factor$s, z_avg=z_avg, 
                  gamma2=gamma2, o_gamma=o_gamma, corpus_topic=factor$corpus_topic, phi=phi),
        sumvar=sumvar, objval=objval, sampvar=sampvar
    );
    return(output);
}


###
### MC_predict: Prediction using Monte-Carlo mean
###
### IN&OUT: factor  = list(alpha, beta, gamma, u, v, s, corpus_topic);
###                    input: initial factor values
###                   output: sample means of alpha, beta, gamma, u, v, s
###                           last sample of corpus_topic
###                   THE VALUES in this list WILL BE CHANGED BY THIS FUNCTION (CALL BY REFERENCE!!)
###
### INPUT:  obs.train = data.frame(y, user, item);
###         obs.test  = data.frame(y, user, item);
###         corpus    = data.frame(item, term, weight);
###         feature   = list(x_user, x_item);
###         x_dyad.train = dyadic features for the training data
###         x_dyad.test  = dyadic features for the test data
###         param     = list(b, g0, d0, c0, G, D, H, var_y, var_alpha, var_beta, var_gamma, var_u, var_v, var_s, lambda, eta);
###         
###         NOTE: When var_FACTOR == NULL, no sample will be drawn for that FACTOR
###         
### OUTPUT: mean    = list(alpha, beta, gamma, u, v, s, z_avg, phi);
###         pred.y  = predicted mean for obs.test
###         rmse    = predictive root mean squared error
###         mae     = predictive mean absolute error
###
### OPTION: useTopicProb [TRUE/FALSE]
###         useTopicProb == FALSE: draw topic samples z_avg[j,k] is avg_w 1{term w in item j belongs to topic k}
###         useTopicProb == TRUE:  z_avg[j,k] =  avg_w Pr(term w in item j belongs to topic k)
###                                corpus_topic[w,k] = Pr(term w in corpus$item[w] belongs to topic k)
###                                No multinomial sampling at all!!
###
### NOTE: factor$alpha, factor$beta and factor$gamma cannot be NULL!!
###
MC_predict.C <- function(
    factor, obs.train, obs.test, corpus, feature, param, x_dyad.train, x_dyad.test,
    nSamples, nBurnIn=1, useTopicProb=TRUE, nRatingExponent=0, transduction=FALSE,
    debug=0, verbose=0
){
    if(useTopicProb && !is.null(factor$corpus_topic) && !is.matrix(factor$corpus_topic)){
        corpus_topic = matrix(as.double(0),nrow=nrow(corpus),ncol=ncol(factor$s));
        corpus_topic[cbind(1:nrow(corpus_topic), factor$corpus_topic)] = 1;
        factor$corpus_topic = corpus_topic;
    }
    feature$x_dyad = x_dyad.train; # just for syncheck.LDA_RLFM.spec
    size = syncheck.LDA_RLFM.spec(factor=factor, obs=obs.train, corpus=corpus, feature=feature, param=param, 
                                  is.corpus_topic.matrix=useTopicProb);
    nTestCases = nrow(obs.test);
    
    corpus_topic_vector = NULL;
    corpus_topic_matrix = NULL;
    if(size$nTopics > 0){
        if(useTopicProb) corpus_topic_matrix = factor$corpus_topic
        else             corpus_topic_vector = factor$corpus_topic;
    }
    
    if(max(obs.test$item) > size$nItems) stop("max(obs.test$item) > size$nItems");
    if(max(obs.test$user) > size$nUsers) stop("max(obs.test$user) > size$nUsers");
    if(!is.integer(obs.test$item)) stop("!is.integer(obs.test$item)");
    if(!is.integer(obs.test$user)) stop("!is.integer(obs.test$user)");
    if(nrow(x_dyad.test) != nTestCases)  stop("nrow(x_dyad.test) != nTestCase");
    if(ncol(x_dyad.test) != ncol(x_dyad.train)) stop("ncol(x_dyad.test) != ncol(x_dyad.train)");
    if(length(nRatingExponent) != 1) stop("length(nRatingExponent) != 1");
    
    z_avg = matrix(as.double(0), nrow=size$nItems, ncol=size$nTopics);
    pred.y  = double(nTestCases);
    
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
    
    if(transduction){
        problem.dim = c(size$nObs, size$corpusSize, size$nUsers, size$nItems, size$nTerms, 
                    size$nFactors, size$nTopics, size$nCorpusWeights,
                    size$nVar_y, size$nVar_alpha, size$nVar_beta, size$nVar_gamma, size$nVar_u, size$nVar_v, size$nVar_s,
                    nTestCases);
        if(!is.integer(problem.dim)) stop("!is.integer(problem.dim)");

        ans = .C("MC_predict",
            # INPUT (initial factor values) & OUTPUT (Monte Carlo mean of factor values)
            #  Required:
            factor$alpha, factor$beta, factor$gamma,
            #  Optional:
            factor$u, factor$v, factor$s, corpus_topic_matrix, corpus_topic_vector,
            # OUTPUT
            z_avg, pred.y, phi,
            # INPUT
            as.integer(nSamples), as.integer(nBurnIn),
            obs.train$user, obs.train$item, corpus$item, corpus$term, corpus$weight, 
            obs.train$y, xb.train, g0x_user, d0x_item, c0x_user, Gx_user, Dx_item, Hx_user,
            param$var_y, param$var_alpha, param$var_beta, param$var_gamma, param$var_u, param$var_v, param$var_s,
            param$eta, param$lambda, 
            obs.test$user, obs.test$item, xb.test,
            as.double(nRatingExponent),
            problem.dim, as.integer(length(problem.dim)),
            as.integer(if(useTopicProb) 1 else 0),
            # OTHER
            as.integer(debug), as.integer(verbose),
            DUP=FALSE
        );
    }else{

        corpus.train = NULL; corpus.newItem = NULL;
        corpus_topic_vector.train = NULL;
        corpus_topic_matrix.train = NULL;
        corpus_topic_vector.newItem = NULL;
        corpus_topic_matrix.newItem = NULL;
        corpusSize.train = as.integer(0);    nCorpusWeights.train = as.integer(0);
        corpusSize.newItem = as.integer(0);  nCorpusWeights.newItem = as.integer(0);
        
        if(size$nTopics > 0){
            itemIDs.train = unique(obs.train$item);
            select.corpus.train = corpus$item %in% itemIDs.train;
            select.corpus.newItem = !select.corpus.train;
            
            corpus.train    = corpus[select.corpus.train,];   corpusSize.train  = nrow(corpus.train);
            corpus.newItem = corpus[select.corpus.newItem,];  corpusSize.newItem = nrow(corpus.newItem);
            if(!is.null(corpus$weight)){
                nCorpusWeights.train   = nrow(corpus.train);
                nCorpusWeights.newItem = nrow(corpus.newItem);
            }
            if(useTopicProb){
                corpus_topic_matrix.train   = corpus_topic_matrix[select.corpus.train,];
                corpus_topic_matrix.newItem = corpus_topic_matrix[select.corpus.newItem,];
            }else{
                corpus_topic_vector.train   = corpus_topic_vector[select.corpus.train];
                corpus_topic_vector.newItem = corpus_topic_vector[select.corpus.newItem];
            }
        }
        
        problem.dim = c(size$nObs, corpusSize.train, size$nUsers, size$nItems, size$nTerms, 
                    size$nFactors, size$nTopics, nCorpusWeights.train,
                    size$nVar_y, size$nVar_alpha, size$nVar_beta, size$nVar_gamma, size$nVar_u, size$nVar_v, size$nVar_s,
                    nTestCases, corpusSize.newItem, nCorpusWeights.newItem);
        if(!is.integer(problem.dim)) stop("!is.integer(problem.dim)");
        

        ans = .C("MC_predict2",
            # INPUT (initial factor values) & OUTPUT (Monte Carlo mean of factor values)
            #  Required:
            factor$alpha, factor$beta, factor$gamma,
            #  Optional:
            factor$u, factor$v, factor$s, 
            corpus_topic_matrix.train, corpus_topic_vector.train,
            corpus_topic_matrix.newItem, corpus_topic_vector.newItem,
            # OUTPUT
            z_avg, pred.y, phi,
            # INPUT
            as.integer(nSamples), as.integer(nBurnIn),
            obs.train$user, obs.train$item, corpus.train$item, corpus.train$term, corpus.train$weight, 
            obs.train$y, xb.train, g0x_user, d0x_item, c0x_user, Gx_user, Dx_item, Hx_user,
            param$var_y, param$var_alpha, param$var_beta, param$var_gamma, param$var_u, param$var_v, param$var_s,
            param$eta, param$lambda, 
            obs.test$user, obs.test$item, xb.test,
            corpus.newItem$item, corpus.newItem$term, corpus.newItem$weight, 
            as.double(nRatingExponent),
            problem.dim, as.integer(length(problem.dim)),
            as.integer(if(useTopicProb) 1 else 0),
            # OTHER
            as.integer(debug), as.integer(verbose),
            DUP=FALSE
        );
    }

    if(useTopicProb) corpus_topic = corpus_topic_matrix
    else             corpus_topic = corpus_topic_vector;
    
    rmse   = sqrt(mean( (obs.test$y - pred.y)^2 ));
    mae    = mean( abs(  obs.test$y - pred.y) );

    output = list(
        mean=list(alpha=factor$alpha, beta=factor$beta, gamma=factor$gamma, 
                  u=factor$u, v=factor$v, s=factor$s, z_avg=z_avg, 
                  corpus_topic=corpus_topic, phi=phi),
        pred.y=pred.y,
        rmse=rmse, mae=mae
    );
    return(output);
}


###
### Gibbs_LDA_EStep.C: See Gibbs_LDA(...) in MCEM_EStep.c
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
###         objval  = list(eta, lambda);
###
Gibbs_LDA_EStep.C <- function(
    corpus_topic, corpus, param, try, nTopics, nSamples, nBurnIn=1,
    debug=0, verbose=0
){
    corpusSize = length(corpus_topic);
    nItems     = as.integer(max(corpus$item));
    nTerms     = as.integer(max(corpus$term)); 
    nCorpusWeights    = if(is.null(corpus$weight)) as.integer(0) else corpusSize;
    nEtaCandidates    = length(try$eta);
    nLambdaCandidates = length(try$lambda);
    
    if(corpusSize != nrow(corpus)) stop("length(corpus_topic) != nrow(corpus)");
    if(max(corpus_topic) > nTopics) stop("max(corpus_topic) > nTopics");
    if(max(corpus_topic) != nTopics) warning("max(corpus_topic) != nTopics");
    
    z_avg = matrix(as.double(0), nrow=nItems, ncol=nTopics);
    phi = matrix(0.0, nrow=nTopics, ncol=nTerms);

    if(!is.double(z_avg)) stop("z_avg must be double");
    
    if(is.null(try)) try = list();
    if(is.null(try[["lambda"]])) try[["lambda"]] = double(0)
    else if(!is.double(try[["lambda"]])) stop("!is.double(try$lambda)");
    if(is.null(try[["eta"]])) try[["eta"]] = double(0)
    else if(!is.double(try[["eta"]])) stop("!is.double(try$eta)");
    
    if(!is.double(param$lambda)) stop("!is.double(param$lambda)");
    if(!is.double(param$eta)) stop("!is.double(param$eta)");
    
    objval  = list(eta=double(length(try$eta)), lambda=double(length(try$lambda)));
        
    problem.dim = as.integer(c(corpusSize, nItems, nTerms, nTopics, nCorpusWeights, nEtaCandidates, nLambdaCandidates));
    if(!is.integer(problem.dim)) stop("!is.integer(problem.dim)");
    if(!is.integer(corpus_topic)) stop("!is.integer(corpus_topic)");
    
    if(!is.integer(corpus$item)) stop("!is.integer(corpus$item)");
    if(!is.integer(corpus$term)) stop("!is.integer(corpus$term)");
    if(!is.null(corpus$weight) && !is.double(corpus$weight)) stop("!is.double(corpus$weight)");
    
    ans = .C("Gibbs_LDA",
        # INPUT (initial topic values) & OUTPUT
        corpus_topic,
        # OUTPUT
        z_avg, phi, objval$eta, objval$lambda,
        # INPUT
        as.integer(nSamples), as.integer(nBurnIn),
        corpus$item,  corpus$term,   corpus$weight, 
        param$eta,    param$lambda,  try$eta,        try$lambda,
        problem.dim, as.integer(length(problem.dim)),
        # OTHER
        as.integer(debug), as.integer(verbose),
        DUP=FALSE
    )
        
    output = list(
        corpus_topic=corpus_topic,
        mean=list(z_avg=z_avg, phi=phi),
        objval=objval
    );
    return(output);
}


getTopicCounts <- function(corpus, corpus_topic, nItems, nTopics, nTerms, debug=0){
    cnt_item_topic = matrix(0.0, nrow=nItems,  ncol=nTopics); 
    cnt_topic_term = matrix(0.0, nrow=nTopics, ncol=nTerms);
    cnt_topic = rep(0.0, nTopics);
    z_avg = matrix(0.0, nrow=nItems, ncol=nTopics); 
    nCorpusWeights = 0;
    if(!is.null(corpus$weight)){
        nCorpusWeights = nrow(corpus);
        if(!is.double(corpus$weight)) stop("error");
    }
    
    if(is.matrix(corpus_topic)){
        if(!is.double(corpus_topic)) stop("error");
        if(nrow(corpus_topic) != nrow(corpus)) stop("error");
        if(ncol(corpus_topic) != nTopics) stop("error");
    }else{
        if(!is.integer(corpus_topic)) stop("error");
    }
    if(!is.integer(corpus$item)) stop("error");
    if(!is.integer(corpus$term)) stop("error");
    
    output = list(
        cnt_item_topic = cnt_item_topic,
        cnt_topic_term = cnt_topic_term,
        cnt_topic = cnt_topic,
        z_avg = z_avg
    );
    if(!is.matrix(corpus_topic)){
        ans = .C("fillInTopicCounts",
            # OUTPUT
            output$cnt_item_topic, output$cnt_topic_term, output$cnt_topic, output$z_avg,
            # INPUT
            corpus_topic, corpus$item, corpus$term, corpus$weight,
            as.integer(nItems), as.integer(nrow(corpus)), as.integer(nTopics), as.integer(nTerms), as.integer(nCorpusWeights),
            as.integer(debug),
            DUP=FALSE
        );
    }else{
        ans = .C("fillInTopicCounts2",
            # OUTPUT
            output$cnt_item_topic, output$cnt_topic_term, output$cnt_topic, output$z_avg,
            # INPUT
            corpus_topic, corpus$item, corpus$term, corpus$weight,
            as.integer(nItems), as.integer(nrow(corpus)), as.integer(nTopics), as.integer(nTerms), as.integer(nCorpusWeights),
            as.integer(debug),
            DUP=FALSE
        );
    }
    output$cnt_item = sum_margin(output$cnt_item_topic, 1);
    
    return(output);
}


indexWithQuantities <- function(x){
    len = sum(x);
    out = integer(len);
    ans = .C("indexWithQuantities",
        out, as.integer(x), as.integer(length(x)),
        DUP=FALSE
    );
    return(out);
}

