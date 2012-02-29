### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
###
### Author: Liang Zhang

# doesn't work yet
logexp <- function(eta){
  pos <- eta > 0
  ans <- rep(NA,length(eta))
  ans[pos] <- eta[pos] + log(1 + exp(-eta[pos]))
  ans[!pos] <- log(1 + exp(eta[!pos]))
  ans
}

getg <- function(b,x,y,offset,spline,rho,penalize.firstcolumn){
  eta <- x%*%b+offset
  p = 2*spline/(1+exp(-2*(1-spline)*eta));
  p[eta>=0] = 2*spline - 1 + 2*(1-spline)/(1+exp(-2*spline*eta[eta>=0]));
  s = - sum(y*log(p)+(1-y)*log(1-p));
  s = s + rho*sum(b^2)
  if (!penalize.firstcolumn) s = s - rho*b[1]*b[1]
  s
}

gradg <- function(b,x,y,offset,spline,rho,penalize.firstcolumn){
   eta <- x%*%b+offset
   q <- y*2*(1-spline)/(1+exp(2*(1-spline)*eta)) + (y-1)*4*spline*(1-spline)/(1+exp(2*(1-spline)*eta))/(1-2*spline+exp(-2*(1-spline)*eta));
   q[eta>=0] <- y[eta>=0]*4*spline*(1-spline)/(1+(2*spline-1)*exp(-2*spline*eta[eta>=0]))/(1+exp(2*spline*eta[eta>=0]))+(y[eta>=0]-1)*spline/(1-spline)/(1+exp(-2*spline*eta[eta>=0]));
   grad = -t(x)%*%q
   t = 2*rho*b;
   if (!penalize.firstcolumn) t[1] = 0;
   grad + t
 }

# y is the binary response
# x is the design matrix including intercept
# beta saves the initial values of the coefficients
# offset is the offset of the observation, default NULL (no offset), length = length(y)
# spline is the spline knots for logistic regression, default 0.5, 0<spline<1
# rho is the L2 penalty parameter, the Loss function = - Loglik + rho*sum(beta^2), default 0 (no penalty)
# penalize.firstcolumn indicates whether we want to penalize the first column of beta (usually intercept). default TRUE.
# lower is the lower bound of each coefficient including intercept, If a scalar, lower = rep(lower,length(y))
# upper is the upper bound of each coefficient including intercept. Can also be a scalar
# maxit is the maximum number of iterations for lbfgs
bounded.logistic<-function(y,x,beta,offset=NULL,spline=0.5,rho=0,penalize.firstcolumn=T,lower=-Inf,upper=Inf,maxit=1000)
{
	if (!is.vector(y)) stop("The response y has to be a vector!");
	if (!is.matrix(x)) stop("The design matrix x has to be a matrix!");
	if (!is.vector(beta)) stop("The initial values of beta has to be a vector!");
	if (ncol(x)!=length(beta)) stop("ncol(x)!=length(beta)");
	if (is.null(offset)) {
	   offset = rep(0,length(y));
	} else {
	   if (length(offset)!=length(y)) stop("length(offset)!=length(y)");
	}
	if (rho<0) stop ("rho<0");
	optim(par=beta,fn=getg,gr=gradg,x=x,y=y,offset=offset,spline=spline,rho=rho,penalize.firstcolumn=penalize.firstcolumn,method="L-BFGS-B",lower=lower,upper=upper, control=list(maxit=maxit))
}

