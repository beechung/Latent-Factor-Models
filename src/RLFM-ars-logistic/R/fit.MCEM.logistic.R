### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
###
### Author: Liang Zhang
###
### Fit a feature-based factorization model using
### Monte-Carlo EM (and Gibbs sampling)
###
### Initialization of parameters should be done before calling this method.
###
### Input:
###   data.train$obs     = data.frame(src.id, dst.id, y);
###                          src.id: source node (user) ID
###                          dst.id: destination node (item) ID
###   data.train$feature = list(x_obs, x_src, x_dst);
###                          x_obs (nObs x nObsFeatures): observation features
###                          x_src (nSrcNodes x nSrcFeatures): features of source nodes (users)
###                          x_dst (nDstNodes x nDstFeatures): features of destination nodes (items)
###   init.model$factor  = list(alpha, beta, u, v);
###   init.model$param   = list(b, g0, d0, G, D, var_alpha, var_beta, var_u, var_v);
###
### Output:
###   output$est.highestCDL: Parameter estimates with the highest complete data likelihood (CDL)
###   output$est.last:       Parameter estimates at the end of the last iteration
###
fit.ARS.logistic <- function(
    nIter,      # Number of EM iterations
    nSamples, nBurnin,   # Number of samples and burnin drawn in each E-step: could be a vector of size nIter.
    data.train, # Training data = list(obs, feature)
    nFactors,   # Number of factors (i.e., number of dimensions of u)
    init.model=NULL, # Initial model = list(factor, param). Set to NULL to use the default.
    doMstep=TRUE,
    # initialization parameters, should be very very small
    var_alpha=0.1, var_beta=0.1, var_v=0.05, var_u=0.05,
    # others
    out.level=0,  # out.level=1: Save the parameter values out.dir/est.highestCDL and out.dir/est.last
    out.dir=NULL, # out.level=2: Save the parameter values of each iteration i to out.dir/est.i
    out.append=FALSE,
    debug=0,      # Set to 0 to disable internal sanity checking; Set to 10 for most detailed sanity checking
    verbose=0,    # Set to 0 to disable console output; Set to 10 to print everything to the console
    use.C=TRUE,   # Whether to use the C implementation
    use.C.EStep=TRUE, # Whether to use the C implementation of the E-step (for backward compatibility)
    print.path=FALSE, # When using R Estep, shall we print sample path of a few random effects?
    delta1=.001,
    use.glmnet=TRUE,
    # ARS parameters
    ars_ninit=3, ars_qcent=c(5.0,50.0,95.0), # number of initial points and the quantiles of the initial points
    ars_xl=-5, ars_xu=5, # lower bound and upper bound of ARS samples
    ars_alpha=0.5,
    fit.ars.alpha=FALSE, # whether we want to fit ars_alpha in the M-step
    fit.regression=TRUE, # do we want to update the regression parameters?
    beta.int=FALSE, # do we want to put the intercept in the beta prior?
    main.effects=FALSE, # only fit the main effects. Leave u and v set to 0.
    identifiable=TRUE, # Whether we want the model to be identifiable or not
    ...             # Additional parameters passing into the regression functions (e.g., bayesglm)
){
    obs=NULL; feature=NULL;
    if(!is.null(data.train)){
        if(!all(c("obs", "feature") %in% names(data.train))) stop("Please check input parameter 'data.train' when calling function fit.multicontext or run.multicontext: data.train$obs and data.train$feature cannot be NULL");
        obs=data.train$obs;
        feature=data.train$feature;
        if(!all(c("src.id", "dst.id", "y") %in% names(obs))) stop("Input parameter data.train$obs must have the following columns: src.id, dst.id, y");
        if(!all(c("x_obs", "x_src", "x_dst") %in% names(feature))) stop("Input parameter data.train$feature must be a list of x_obs, x_src, x_dst");
    }else{
        stop("Please specify data.train!");
    }

    if(is.null(init.model)){
        init.model = init.simple.random(
                obs=obs, feature=feature,
                nFactors=nFactors, has.u=TRUE, has.gamma=FALSE, is.logistic=TRUE,
                var_alpha=var_alpha, var_beta=var_beta, var_v=var_v, var_u=var_u
        );
        init.model$param$xi = NULL;
    }
    factor = init.model$factor;
    param  = init.model$param;
    if(!all(c("alpha", "beta", "u", "v") %in% names(factor))) stop("init.model$factor must be a list of alpha, beta, u, v");
    if(!all(c("b", "g0", "d0", "G", "D", "var_alpha", "var_beta", "var_u", "var_v") %in% names(param))) stop("init.model$param must be a list of b, g0, d0, G, D, var_alpha, var_beta, var_u, var_v");

    output = fit.MCEM.logistic(
            nIter=nIter,
            doMstep=doMstep,
            nSamples=nSamples, nBurnin=nBurnin,
            user=obs$src.id, item=obs$dst.id, y=obs$y,
            x=feature$x_obs, w=feature$x_src, z=feature$x_dst,
            alpha=factor$alpha, beta=factor$beta, u=factor$u, v=factor$v,
            b=param$b, g0=param$g0, G=param$G, d0=param$d0, D=param$D,
            var_alpha=param$var_alpha, var_beta=param$var_beta, var_u=param$var_u, var_v=param$var_v,
            out.level=out.level, out.dir=out.dir, out.append=out.append,
            debug=debug, verbose=verbose, use.C=use.C, use.C.EStep=use.C.EStep,
            print.path=print.path, delta1=delta1, use.glmnet=use.glmnet,
            ars_ninit=ars_ninit, ars_qcent=ars_qcent,
            ars_xl=ars_xl, ars_xu=ars_xu, ars_alpha=ars_alpha,
            fit.ars.alpha=fit.ars.alpha, fit.regression=fit.regression,
            beta.int=beta.int, center=center, main.effects=main.effects,
            ...
    );

    return(output);
}

fit.MCEM.logistic <- function(
    nIter,      # Number of EM iterations
    doMstep,
    nSamples, nBurnin,   # Number of samples and burnin drawn in each E-step: could be a vector of size nIter.
    user, item, y, x, w, z,  # Observed rating and feature values
    alpha, beta, u, v,       # Main effects and factors
    b, g0, G, d0, D,         # Feature weights (regression parameters)
    var_alpha, var_beta, var_u, var_v=1, # Variances
    out.level=0,  # out.level=1: Save the parameter values out.dir/est.highestCDL and out.dir/est.last
    out.dir=NULL, # out.level=2: Save the parameter values of each iteration i to out.dir/est.i
    out.append=FALSE,
    debug=0,      # Set to 0 to disable internal sanity checking; Set to 10 for most detailed sanity checking
    verbose=0,    # Set to 0 to disable console output; Set to 10 to print everything to the console
    use.C=TRUE,   # Whether to use the C implementation
    use.C.EStep=TRUE,     # Whether to use the C implementation of the E-step (for backward compatibility)
    print.path=F, # when using R Estep, shall we print sample path of a few random effects?
    delta1=.001,
    use.glmnet=F,
    ars_ninit=3, ars_qcent=c(5.0,50.0,95.0),
    ars_xl=-5, ars_xu=5, ars_alpha=0.5,
    fit.ars.alpha=F, # whether we want to fit ars_alpha in the M-step
    fit.regression=T, # do we want to update the regression parameters?
    beta.int=F, # do we want to put the intercept in the beta prior?
    main.effects=F, # only fit the main effects. Leave u and v set to 0.
    identifiable=TRUE, # Whether we want the model to be identifiable or not
    ...         # Additional parameters passing into the regression functions (e.g., bayesglm)
){
    user = as.integer(user);
    item = as.integer(item);

    #if(beta.int && center) stop("Cannot learn intercept in random effects when centered")

    if (main.effects)
    {
        u = matrix(0, dim(u)[1], dim(u)[2])
        v = matrix(0, dim(v)[1], dim(v)[2])
    }

    if (identifiable) {
      # Make sure initialization of v are all positive
      v = abs(v);
      beta.int = FALSE;
      center = TRUE;
    } else {
      center = FALSE;
    }

    if(beta.int && length(d0) != dim(z)[2] + 1) d0 = c(0, d0)

    if(use.C) use.C.EStep = TRUE;

    if(length(nSamples)!=nIter) nSamples <- rep(nSamples[1],nIter);
    if(length(nBurnin)!=nIter) nBurnin <- rep(nBurnin[1],nIter);
    if(debug >= 1) check.input.logistic(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_alpha, var_beta, var_u, var_v);
    nObs     = length(y);
    nUsers   = length(alpha);
    nItems   = length(beta);
    nFactors = ncol(u);

    LL = rep(NA, nIter+1); # LL records the logLikelihood of each iteration
    bestLL = logLikelihood.logistic(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_alpha, var_beta, var_u, var_v, ars_alpha, beta.int, debug, use.C.EStep);
    LL[1] = bestLL;

    if (length(var_u)==1) var_u = rep(var_u,nFactors);
    if (length(var_v)==1) var_v = rep(var_v,nFactors);

    # Paramter estimates with the highest complete data likelihood (CDL)
    est.highestCDL = list(
        alpha=alpha, beta=beta, u=u, v=v,
        b=b, g0=g0, G=G, d0=d0, D=D,
        var_alpha=var_alpha, var_beta=var_beta, var_u=var_u, var_v=var_v
    );
    if(verbose >= 1){
        cat("START fit.MCEM (initial logLikelihood: ",bestLL," + constant)\n",sep="");
        cat("Alpha for logistic spline = ",ars_alpha,"\n");
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
        save(file=paste(out.dir,"/est.0",sep=""),
             list=c("alpha", "beta", "u", "v", "b", "g0", "G", "d0", "D",
                    "var_alpha", "var_beta", "var_u", "var_v"));
        summary = data.frame(Method="MCEM", Iter=0, nSteps=nSamples[1], CDlogL=LL[1], TimeUsed1=0, TimeUsed2=0);
        file = paste(out.dir,"/summary",sep="");
        if(file.exists(file))
            write.table(summary, file=file, append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
        else
            write.table(summary, file=file, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE);
    }
    # Initialize ars_XI_alpha, etc.
    xi = rep(0,ars_ninit);
    for(i in 1:ars_ninit){
      xi[i] <- ars_xl + i*(ars_xu - ars_xl)/(ars_ninit + 1.0);
      if(xi[i] >= ars_xu) xi[i]=xi[i] - .1;
      if(xi[i] <= ars_xl) xi[i] = xi[i] + .1;
    }
    ars_XI_alpha = rep(xi,nUsers);
    ars_XI_beta = rep(xi,nItems);
    ars_XI_u = rep(xi,nUsers*nFactors);
    if (identifiable) {
      # Initialize ars_XI_v as positive
      xi = rep(0,ars_ninit);
      for(i in 1:ars_ninit){
        xi[i] <- i*(ars_xu)/(ars_ninit + 1.0);
        if(xi[i] >= ars_xu) xi[i]=xi[i] - .1;
        if(xi[i] <= 0) xi[i] = xi[i] + .1;
      }
    }
    ars_XI_v = rep(xi,nItems*nFactors);

    begin.time = proc.time();

    for(iter in 1:nIter){

        if(verbose >= 2){
            cat("---------------------------------------------------------\n",
                "Iteration ",iter,"\n",
                "---------------------------------------------------------\n",
                "start E-STEP\n",sep="");
        }
        b.time = proc.time();

        ###
        ### E-STEP
        ###
        if(use.C.EStep){
            mc_e = MC_EStep_logistic_arscid.C(
              nSamples[iter], nBurnin[iter],
              user, item, y, x, w, z,
              alpha, beta, u, v,
              ars_XI_alpha, ars_XI_beta,
              ars_XI_u, ars_XI_v,
              b, g0, G, d0, D,
              var_alpha, var_beta, var_u, var_v,
              ars_ninit, ars_qcent, ars_xl, ars_xu, ars_alpha,
              debug, beta.int, center, main.effects
            );
        }else{
            mc_e = MC_EStep_logistic.R(
              nSamples[iter], nBurnin[iter],
              user, item, y, x, w, z,
              alpha, beta, u, v,
              ars_XI_alpha, ars_XI_beta,
              ars_XI_u, ars_XI_v,
              b, g0, G, d0, D,
              var_alpha, var_beta, var_u, var_v,
              ars_ninit, ars_qcent, ars_xl, ars_xu, ars_alpha,
              debug,
              print.path=print.path
              );
        }
        o     = mc_e$o.mean;
        alpha = mc_e$alpha.mean;
        beta  = mc_e$beta.mean;
        u     = mc_e$u.mean;
        v     = mc_e$v.mean;

        acceptrate.maineff = mc_e$acceptrate.maineff;
        acceptrate.fact = mc_e$acceptrate.fact;
        if(print.path){u11 <- mc_e$u11;u21 <- mc_e$u21;v11 <- mc_e$v11;v21 <- mc_e$v21;}
        time.used = proc.time() - b.time;
        time.used.1 = time.used;
        if(verbose >= 2){
            ll = logLikelihood.logistic(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_alpha, var_beta, var_u, var_v, ars_alpha=ars_alpha, beta.int, debug=debug, use.C=use.C.EStep);
            cat("end   E-STEP (logLikelihood = ",ll," + constant,  used ",time.used[3]," sec)\n",
                "start M-STEP\n",sep="");
        }
        b.time = proc.time();

        ###
        ### M-STEP
        ###

        use.lm=F;

        bold <- b;
        g0old <- g0;      var_alphaold <- var_alpha;
        Gold <- G;        var_uold <- var_u;
        d0old <- d0;      var_betaold <- var_beta;
        Dold <- D;        var_vold <-  var_v;
        if(doMstep==1){
            mc_m = MC_MStep_logistic_arscid(
                user, item, y, x, b, w, z, o,
                alpha=alpha, alpha.sumvar=mc_e$alpha.sumvar, beta=beta, beta.sumvar=mc_e$beta.sumvar,
                u=u, u.sumvar=mc_e$u.sumvar, v=v, v.sumvar=mc_e$v.sumvar,
                debug=debug, lm=use.lm,
                use.glmnet=use.glmnet, fit.ars.alpha=fit.ars.alpha,
                fit.regression=fit.regression, beta.int=beta.int,
                main.effects=main.effects, ...
            );

            if (fit.ars.alpha) {
                b  = mc_m$b;
                ars_alpha = mc_m$ars_alpha;
            }
            if (fit.regression) {
                g0 = mc_m$g0;
                G  = mc_m$G;
                d0 = mc_m$d0;
                D  = mc_m$D;
            }
            if (beta.int) {
                d0[1] = mc_m$d0[1];
            }
            var_alpha = mc_m$var_alpha;
            var_u     = mc_m$var_u;
            var_v     = mc_m$var_v;
            var_beta  = mc_m$var_beta;
            b = mc_m$b;

            if (main.effects) {
                G = matrix(0, nrow=ncol(w), ncol=nFactors)
                D = matrix(0, nrow=ncol(z), ncol=nFactors)
                var_u = rep(1,nFactors)
                var_v = rep(1,nFactors)
            }
        }
        if (identifiable) {
    	    # Order G, D, u and v
    	    ind = order(var_v);
    	    var_u = var_u[ind];
    	    var_v = var_v[ind];
    	    G = G[,ind];
    	    D = D[,ind];
    	    u = u[,ind];
    	    v = v[,ind];
        }
        if(verbose > 0){
          r <- dim(D)[2]
          cat("var_alpha=",var_alpha,"\n");
          cat("var_beta=",var_beta,"\n");
          cat("var_u=",var_u,"\n");
          cat("var_v=",var_v,"\n");
          cat("b=",b,"\n");
          if(beta.int) cat("mean(beta)=",mean(beta),"\n")
          cat("bdiff =",max(abs(b-bold)/(abs(bold) + delta1)),"\n")
          cat("g0diff=",max(abs(g0 - g0old)/(abs(g0old) + delta1)),"\n")
          cat("d0diff=",max(abs(d0 - d0old)/(abs(d0old) + delta1)),"\n")
          cat("varudiff=",max(abs(var_u*var_v - var_uold*var_vold )/(abs(var_uold*var_vold) + delta1)),"\n")
          cat("varmaineffectsdiff=",max(abs(var_alpha*var_beta - var_alphaold*var_betaold )/(abs(var_alphaold*var_betaold) + delta1)),"\n")
          if(!main.effects) {
              trGold <- sum(svd(Gold)$d); trG <- sum(svd(G)$d);trDold <- sum(svd(Dold)$d); trD <- sum(svd(D)$d)
              cat("G0diff=",max(abs(trG - trGold)/(abs(trGold) + delta1)),"\n")
              cat("D0diff=",max(abs(trD - trDold)/(abs(trDold) + delta1)),"\n")
          }
          o2 = x %*% b + o; ep = exp(o2)/(1+exp(o2))
          cat("E(P(click)) = ", mean(ep), "\n")
        }
        LL[iter+1] = logLikelihood.logistic(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_alpha, var_beta, var_u, var_v, ars_alpha, beta.int, debug, use.C.EStep);

        time.used = proc.time() - b.time;
        time.used.2 = time.used;
        if(verbose >= 2){
            cat("end   M-STEP (logLikelihood = ",LL[iter+1]," + constant,  used ",time.used[3]," sec)\n",sep="");
        }
        ###
        ### Update the model (output) if the logLikelihood is increased
        ###
        if(LL[iter+1] > bestLL){
            bestLL = LL[iter+1];
            est.highestCDL$alpha[]=alpha; est.highestCDL$beta[]=beta; est.highestCDL$u[,]=u; est.highestCDL$v[,]=v;
            est.highestCDL$b[]=b; est.highestCDL$g0[]=g0; est.highestCDL$G[,]=G; est.highestCDL$d0[]=d0; est.highestCDL$D[,]=D;
            est.highestCDL$var_alpha=var_alpha; est.highestCDL$var_beta=var_beta; est.highestCDL$var_u=var_u; est.highestCDL$var_v=var_v;

            if(verbose >= 2){
                cat("Found a better fit!\n");
            }
        }

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
            save(file=file,
                 list=c("o","alpha", "beta", "u", "v", "b", "g0", "G", "d0", "D",
                        "var_alpha", "var_beta", "var_u", "var_v"));
            if(LL[iter+1] == bestLL) file.copy(file, paste(out.dir,"/est.highestCDL",sep=""), overwrite=TRUE);
            summary = data.frame(Method="MCEM", Iter=iter, nSteps=nSamples[iter], CDlogL=LL[iter+1], TimeUsed1=time.used.1[3], TimeUsed2=time.used.2[3]);
            file = paste(out.dir,"/summary",sep="");
            write.table(summary, file=file, append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE);
        }
    }
    output = list(
        trainingCDL = LL,
        highestCDL  = bestLL,
        est.highestCDL = est.highestCDL,
        est.last = list(
            alpha=alpha, beta=beta, u=u, v=v,
            b=b, g0=g0, G=G, d0=d0, D=D,
            var_alpha=var_alpha, var_beta=var_beta, var_u=var_u, var_v=var_v
        )
    );
    if(print.path)output$path.out <- list(u11=u11,u21=u21,v11=v11,v21=v21);

    time.used = proc.time() - begin.time;
    if(verbose >= 1){
        cat("END   fit.MCEM (final   logLikelihood: ",bestLL," + constant,  used ",time.used[3]," sec)\n",sep="");
    }

    return(output);
}

fit.MCEM.factorOnly <- function(
    nIter,      # Number of EM iterations
    nSamples,   # Number of samples drawn in each E-step: could be a vector of size nIter.
    user, item, y,           # Observed rating
    alpha, beta, u, v,       # Main effects and factors
    var_y, var_alpha, var_beta, var_u, var_v=1, # Variances
    out.level=0,  # out.level=1: Save the parameter values out.dir/est.highestCDL and out.dir/est.last
    out.dir=NULL, # out.level=2: Save the parameter values of each iteration i to out.dir/est.i
    out.append=FALSE,
    debug=0,      # Set to 0 to disable internal sanity checking; Set to 10 for most detailed sanity checking
    verbose=0,    # Set to 0 to disable console output; Set to 10 to print everything to the console
    use.C=TRUE,   # Whether to use the C implementation
    use.C.EStep=TRUE,     # Whether to use the C implementation of the E-step (for backward compatibility)
    func.rmvnorm=mvrnorm, # A function that generate a random sample of multivariate normal (only matters when use.C.EStep=FALSE)
    print.path=F, # when using R Estep, shall we print sample path of a few random effects?
    delta1=.001,
    ...         # Additional parameters passing into the regression functions (e.g., bayesglm)
){
    nObs     = length(y);
    nUsers   = length(alpha);
    nItems   = length(beta);
    nFactors = ncol(u);

    return(fit.MCEM(
        nIter=nIter, nSamples=nSamples,
        user=user, item=item, y=y, x=matrix(0,nrow=nObs,ncol=1), w=matrix(0,nrow=nUsers,ncol=1), z=matrix(0,nrow=nItems,ncol=1),
        alpha=alpha, beta=beta, u=u, v=v,
        b=0, g0=0, G=matrix(0,nrow=1,ncol=nFactors), d0=0, D=matrix(0,nrow=1,ncol=nFactors),
        var_y=var_y, var_alpha=var_alpha, var_beta=var_beta, var_u=var_u, var_v=var_v,
        out.level=out.level, out.dir=out.dir, out.append=out.append,
        debug=debug, verbose=verbose, use.C=use.C, use.C.EStep=use.C.EStep,
        func.rmvnorm=func.rmvnorm, print.path=print.path, delta1=delta1,
        ...
    ));
}


fit.MCEM.factorhybrid <- function(
                                   nIter,      # Number of EM iterations
    nSamples,   # Number of samples drawn in each E-step: could be a vector of size nIter.
    user, item, y, x, w, z,  # Observed rating and feature values
    alpha, beta, u, v,       # Main effects and factors
    b, g0, G, d0, D,         # Feature weights (regression parameters)
    var_y, var_alpha, var_beta, var_u, var_v=1, # Variances
    out.level=0,  # out.level=1: Save the parameter values out.dir/est.highestCDL and out.dir/est.last
    out.dir=NULL, # out.level=2: Save the parameter values of each iteration i to out.dir/est.i
    out.append=FALSE,
    debug=0,      # Set to 0 to disable internal sanity checking; Set to 10 for most detailed sanity checking
    verbose=0,    # Set to 0 to disable console output; Set to 10 to print everything to the console
    use.C=TRUE,   # Whether to use the C implementation
    use.C.EStep=TRUE,     # Whether to use the C implementation of the E-step (for backward compatibility)
    func.rmvnorm=mvrnorm, # A function that generate a random sample of multivariate normal (only matters when use.C.EStep=FALSE)
    print.path=F, # when using R Estep, shall we print sample path of a few random effects?
    delta1=.001,
    ...         # Additional parameters passing into the regression functions (e.g., bayesglm)
){
  fit1 <- fit.MCEM.factorOnly(nIter=nIter,nSamples=nSamples,user=user,item=item,y=y,
                              alpha=alpha,beta=beta,u=u,v=v,
                              var_y=var_y, var_alpha=var_alpha, var_beta=var_beta,
                              var_u=var_u, var_v=var_v, debug=debug, verbose=verbose,
                              use.C.EStep=use.C.EStep,use.C=use.C,
                              func.rmvnorm=func.rmvnorm,
                              print.path=F,
                              delta1=delta1,...)
  #now run regression on the fitted factors to populate g0,d0,G,D.


}
