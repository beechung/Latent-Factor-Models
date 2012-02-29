###
### Copyright (c) 2012, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### Author: Deepak Agarwal
###

get.sp.des <- function(X){
  #getting a sparse representation of design matrix X: obsid fval nfobs
  ncov <- dim(X)[2];nobs <- dim(X)[1]
  obsid <- fval <- nfobs <- NULL
  for(i in 1:ncov){
    f <- X[,i]
    wk <- c(1:nobs)[f != 0]
    obsid <- c(obsid,wk); fval <- c(fval,f[wk]);nfobs <- c(nfobs,length(wk))
  }
  ans <- list(obsid,fval,nfobs)
  names(ans) <- c("obsid","fval","nfobs")
  ans
}


get.dense.des <- function(feat){
  nobs <- max(feat[,2]) #assumes there are no obs with all features absent.
  ncov <- max(feat[,1]) #assumes there are no features that do not occur at all.
  X <- matrix(0,nrow=nobs,ncol=ncov)
  for(i in 1:ncov)
    X[feat[,2][feat[,1]==i],i] <- feat[,3][feat[,1]==i]
  X
}

sim.logistic <- function(ncov,nobs,mbeta=0,sigbeta=1,zfr=.5){
#simulate the betas (for a fixed ncov): N(m,sig)
  beta <- rnorm(ncov,mean=mbeta,sd=sqrt(sigbeta))
#simulate design matrix nobs x ncov (entry ~ N(0,1/sqrt(ncov)))
  X <- matrix(rnorm(nobs*ncov),ncol=ncov)/sqrt(ncov)
  nzero <- min(ceiling(nobs*zfr),nobs)
  for(i in 1:ncov)X[sample(1:nobs,nzero),i] <- 0
  p <- (1+exp(-X%*%beta))^{-1}
  y <- rbinom(nobs,1,p)
  #also create design matrix in sparse format
  wk <- get.sp.des(X)
  obsid <- wk$obsid
  fval <- wk$fval
  nfobs <- wk$nfobs
  ans <- list(y,X,obsid,fval,nfobs,beta,mbeta,sigbeta)
  names(ans) <- c("y","X","obsid","fval","nfobs","beta","mbeta","sigbeta")
  ans
}

getlgt <- function(x,id,beta,y,X,mbeta,sigb){
  beta[id] <- x
  eta <- X%*%beta
  idx <- c(1:dim(X)[1])[X[,id]!=0]
  llik <- sum(y[idx]*eta[idx] - log1p(exp(eta[idx])))
  lpr <- -.5*sum((x - mbeta[id])^2/sigb[id])
  llik+lpr
}


arslgt <- function(obsid,featval,nfobs,y,intercept=F,betaprior=NULL,vprior=1.0,vprior.intercept=NULL,
                   nsamp=500,burnin=100,erange=10.0,verbose=T,verboseC=0,isC=T){

  cat("vprior=",vprior,"\n")

  y <- as.integer(y);ncov <- as.integer(length(nfobs));nobs <- as.integer(length(y))
  pglb <- length(y[y==1])/nobs
  if(pglb==0 | pglb==1)stop("only one label type in the data")

  off=if(intercept)rep(log(pglb)-log(1.0-pglb),length(y)) else rep(0,length(y))
  
  beta0 <- if(!is.null(betaprior))betaprior else rep(0,ncov)
  if(length(beta0)!=ncov)stop("prior of beta not specified correctly")
  varbeta <- rep(vprior,ncov)
  if(!is.null(vprior.intercept))varbeta[ncov] <- vprior.intercept
  betaout <- beta0

  #parameter settings to run adaptive rejection sampler, no need to change.
  ninit <- as.integer(3); xl <- as.double(-abs(erange));xu <- as.double(abs(erange));xi <- as.double(rep(0.0,ninit))
  eta <- as.double(rep(0,10*nobs));
  
  ncent <- as.integer(3);qcent <- as.double(c(5.0,50.0,95.0));xcent <- as.double(rep(0.0,ncent)); neval <- as.integer(rep(0,ncov))
  for(i in 1:ninit){
    xi[i] <- xl + (i + 1.0)*(xu - xl)/(ninit + 1.0)
    if(xi[i] >= xu)xi[i]=xi[i] - .1;
    if(xi[i] <= xl)xi[i] = xi[i] + .1;
  }
  XI <- as.double(rep(xi,ncov))
  #end of parameter settings.
  Nsamp <- nsamp+burnin
  betamean <- rep(0,ncov); betavar <- rep(0,ncov)
  #sampv <- NEVAL <- matrix(NA,nrow=nsamp,ncol=ncov)

  if(isC){
    if(burnin > 0){
      ans <- .C("ARSSAMP",off,beta0,varbeta,as.integer(obsid),as.double(featval),as.integer(nfobs),y,ncov,nobs,betaout,as.double(betamean),as.double(betavar),qcent,xcent,ncent,ninit,xl,xu,XI,
                eta,neval,as.integer(burnin),as.integer(verboseC))
      betaout <- as.double(ans[[10]])
    }
    betamean <- rep(0,ncov); betavar <- rep(0,ncov)
    ans <- .C("ARSSAMP",off,beta0,varbeta,as.integer(obsid),as.double(featval),as.integer(nfobs),y,ncov,nobs,betaout,as.double(betamean),as.double(betavar),qcent,xcent,ncent,ninit,xl,xu,XI,
                eta,neval,as.integer(nsamp),as.integer(verboseC))
    betamean <- ans[[11]]
    betavar <- ans[[12]] 
  }
  else{
  #sampling.
  
    for(i in 1:Nsamp){
      ans <- .C("ARSLOGISTICSPARSE",off,beta0,varbeta,as.integer(obsid),as.double(featval),as.integer(nfobs),y,ncov,nobs,betaout,qcent,xcent,ncent,ninit,xl,xu,XI,
                eta,neval,as.integer(verboseC))
      sampv <- ans[[10]]
      #sampv[i,] <- ans[[10]]
      #NEVAL[i,] <- ans[[19]]
      betaout <- as.double(sampv)
      if(verbose)cat("sample=",sampv[1:min(ncov,5)],"\n")
      if(i > burnin){
	j <- i - burnin
        betamean = betamean + (sampv - betamean)/j
        betavar = betavar + (sampv*sampv - betavar)/j
	}



      XItemp <- as.double(ans[[17]])
      if(!any(is.na(XItemp)))XI <- XItemp
      rm(ans)
      #if(any(is.na(apply(sampv,2,mean))))stop("numerical problem")
    }
  }
  betavar = betavar - betamean*betamean
  ret <- list(mean=betamean,var=betavar)
  #samples=sampv[(burnin+1):nsamp,,drop=F],intercept=off[1],neval=NEVAL)
  ret
}
