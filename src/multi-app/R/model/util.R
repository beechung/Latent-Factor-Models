### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

##############################################################################
#   Prediction functions
##############################################################################

### Input:
###   response = data.frame(user, app, item,  y, w)
###   param    = list(A, B, b, alpha, beta, var_x, var_y, var_z, var_u)
### Output:
###   mean of response$y based on the input z
### NOTE: response$user > nrow(z[[k]]) are the data of NEW users
predict.y.from.z <- function(response, param, z, add.noise=FALSE){
	warning.any.not.in(c("user", "app", "item"), names(response), msg="response must have: ", stop=TRUE);
	return(predict.obs.from.z(
			user=response$user, app=response$app, index=response$item, 
			B=param$beta, b=param$alpha, z=z,
			w=response$w, var_x=param$var_y, add.noise=add.noise
	));
}

### Input:
###   feature  = data.frame(user, app, index, x, w)
###     param  = list(A, B, b, alpha, beta, var_x, var_y, var_z, var_u)
### Output:
###   mean of feature$x based on the input z
### NOTE: feature$user > nrow(z[[k]]) are the data of NEW users
predict.x.from.z <- function(feature, param, z, add.noise=FALSE){
	warning.any.not.in(c("user", "app", "index"), names(feature), msg="feature must have: ", stop=TRUE);
	return(predict.obs.from.z(
			user=feature$user, app=feature$app, index=feature$index, 
			B=param$B, b=param$b, z=z,
			w=feature$w, var_x=param$var_x, add.noise=add.noise
	));
}

### NOTE: user > nrow(z[[k]]) means NEW users
predict.obs.from.z <- function(user, app, index, B, b, z, w=NULL, var_x=NULL, add.noise=FALSE){
	temp = tapply(seq_len(length(app)), list(app), FUN=c, simplify=FALSE);
	out  = rep(NA, length(user));
	if(length(B) != length(b)) stop("length(B) != length(b)");
	if(length(B) != length(z)) stop("length(B) != length(z)");
	if(add.noise && length(B) != length(var_x)) stop("length(B) != length(var_x)");
	for(ind in seq_len(length(temp))){
		k = as.integer(names(temp)[ind]);
		select = temp[[ind]];
		if(is.null(select)) next;
		i = user[select];
		m = index[select];
		nUsers = max(user);
		if(nUsers > nrow(z[[k]])){
			z_k = rbind(z[[k]], matrix(0, nrow=nUsers-nrow(z[[k]]), ncol=ncol(z[[k]])));
		}else{
			z_k = z[[k]];
		}
		if(length(B[[k]]) == 1 && B[[k]] == 1){
			Bz = z_k[cbind(i,m)] * B[[k]];
		}else{
			Bz = B[[k]][m,,drop=FALSE] * z_k[i,,drop=FALSE];
			Bz = sum_margin(Bz, 1);
		}
		if(is.null(Bz)) stop("Bz is null")
		out[select] = b[[k]][m] + Bz;
		
		if(add.noise){
			if(is.null(w)) var = var_x[k]
			else           var = var_x[k] * w;
			out[select] = out[select] + rnorm(length(select), mean=0, sd=sqrt(var));
		}
	}
	return(out);
}

###
### Prediction based on the Gaussian model
### NOTE: feature and response may contain new users; in this case, notice that
###   (1) The input factor only contains the factors for old users.
###       In particular, consider user i = feature[m,"user"] (the m-th obs).
###       His/her factors are factor$u[i,] and factor$z[k][i,]
###   (2) User IDs of OLD users are 1, ..., nrow(factor$u).
###       User IDs of NEW users are (nrow(factor$u)+1), ..., max(feature$user, response$user).
###
predict.x.and.y <- function(
	feature, response, param, factor
){
	out = list();
	if(!is.null(response)){
		if(is.null(response$w)) weight = 1
		else                    weight = 1/response$w;
		pred.y = predict.y.from.z(response=response, param=param, z=factor$z, add.noise=FALSE);
		temp = predict.response.from.gaussian(pred.y, response$y, 1, param$is.y.logistic, weight);
		# Here, temp$loglik is not right (but not used; so, it's OK)
		out$pred.y = temp$pred.y;
		out$true.y = response$y;
		out$RMSE.y = temp$rmse;
		out$loss.y = temp$test.loss;
	}
	if(!is.null(feature)){
		if(is.null(feature$w)) weight = 1
		else                   weight = 1/feature$w;
		pred.x = predict.x.from.z(feature=feature, param=param, z=factor$z, add.noise=FALSE);
		temp = predict.response.from.gaussian(pred.x, feature$x, 1, param$is.x.logistic, weight);
		# Here, temp$loglik is not right (but not used; so, it's OK)
		out$pred.x = temp$pred.y;
		out$true.x = feature$x;
		out$RMSE.x = temp$rmse;
		out$loss.x = temp$test.loss;
	}
	
	return(out);
}


##############################################################################
#   Likelihood Functions
##############################################################################

###
### Marginal Likelihood
### NOTE: Very computationally inefficient; this is for debugging
###
marginal.loglik.R <- function(feature, response, param, verbose=0){
	
	check.syntax.obs(feature=feature, response=response);
	size = check.syntax.param(param);
	
	obsIndex.f = tapply(seq_len(nrow(feature)),  list(feature$user),  FUN=c, simplify=FALSE);
	obsIndex.r = tapply(seq_len(nrow(response)), list(response$user), FUN=c, simplify=FALSE);
	
	nUsers = max(feature$user, response$user);
	A = NULL;  var_z = NULL;
	for(k in seq_len(size$nApps)){
		if(length(param$A[[k]]) == 1 && param$A[[k]] == 1) A_k = param$A[[k]] * diag(size$nGlobalFactors)
		else                                               A_k = param$A[[k]];
		A = rbind(A, A_k);
		var_z = c(var_z, rep(param$var_z[k], nrow(A_k)));
	}
	AA = A %*% t(A);
	
	ans = 0;
	
	for(i in seq_len(nUsers)){
		this.user = as.character(i); # IMPORTANT: use string to access obsIndex.f and obsIndex.r
		f = feature[ obsIndex.f[[this.user]],];
		r = response[obsIndex.r[[this.user]],];
		if(nrow(f)+nrow(r) == 0) next;
		f.by.app = split(f, f$app);
		r.by.app = split(r, r$app);

		###    log|Sigma| = nGlobalFactors*log(var_u) + log|E| + log|G|
		### d' Sigma^-1 d = L - H G^-1 H'
		###        log|E| = sum_k log|V_ik| + log|D_ik| + log|J_ik|
		###             G = I/var_u + FEF
		###           FEF = sum_k F_ik' D_ik^-1 F_ik - N_ik' J_ik^-1 N_ik
		###          F_ik = C_ik A_ik
		###          N_ik = C_ik' D_ik^-1 F_ik
		###          J_ik = V_ik^-1 + C_ik' D_ik^-1 C_ik
		###             L = sum_k d_ik' D_ik^-1 d_ik - M_ik J_ik^-1 M_ik'
		###          M_ik = d_ik' D_ik^-1 C_ik
		###             H = sum_k M_ik A_ik - M_ik J_ik^-1 N_ik
		log.E=0; FEF=0; L=0; H=0;
		for(k in seq_len(size$nApps)){
			this.app = as.character(k);
			this.f = f.by.app[[this.app]];
			this.r = r.by.app[[this.app]];
			d_xik = NULL; C_xik = NULL; D_xik = NULL;
			d_yik = NULL; C_yik = NULL; D_yik = NULL;
			if(!is.null(this.f)){
				if(length(param$B[[k]]) == 1 && param$B[[k]] == 1) B = param$B[[k]] * diag(size$nLocalFactors[k])
				else                                               B = param$B[[k]];
				d_xik = this.f$x - param$b[[k]][this.f$index];
				C_xik = B[this.f$index,,drop=FALSE];
				if(is.null(this.f$w)) w = rep(1, nrow(this.f))
				else                  w = this.f$w;
				D_xik = w * param$var_x[k];
			} 
			if(!is.null(this.r)){
				d_yik = this.r$y - param$alpha[[k]][this.r$item];
				C_yik = param$beta[[k]][this.r$item,,drop=FALSE];
				if(is.null(this.r$w)) w = rep(1, nrow(this.r))
				else                  w = this.r$w;
				D_yik = w * param$var_y[k];
			}
			d_ik = c(d_xik, d_yik);
			if(length(d_ik) == 0) next;
			D_ik = c(D_xik, D_yik);
			C_ik = rbind(C_xik, C_yik);
			V_ik = rep(param$var_z[k], size$nLocalFactors[k]);
			if(length(param$A[[k]]) == 1 && param$A[[k]] == 1) A_ik = diag(nrow=size$nGlobalFactors)
			else                                               A_ik = param$A[[k]];
			F_ik = C_ik %*% A_ik;
			D.inv.C = array(1/D_ik,dim=dim(C_ik)) * C_ik;
			D.inv.F = array(1/D_ik,dim=dim(F_ik)) * F_ik;
			N_ik = t(D.inv.C) %*% F_ik;
			J_ik = diag(1/V_ik,nrow=length(V_ik)) + t(C_ik) %*% D.inv.C;
			M_ik = t(d_ik) %*% D.inv.C;
			J_ik.inv = solve(J_ik);
			log.E = log.E + sum(log(V_ik)) + sum(log(D_ik)) + log(det(J_ik));
			FEF   = FEF   + t(F_ik) %*% D.inv.F - t(N_ik) %*% J_ik.inv %*% N_ik;
			L     = L     + t(d_ik) %*% (d_ik/D_ik) - M_ik %*% J_ik.inv %*% t(M_ik);
			H     = H     + M_ik %*% A_ik - M_ik %*% J_ik.inv %*% N_ik;
		}
		G = diag(1/param$var_u,nrow=size$nGlobalFactors) + FEF;
		G.inv = solve(G);
		log.det.Sigma = size$nGlobalFactors*log(param$var_u) + log.E + log(det(G));
		d.Sigma.inv.d = L - H %*% G.inv %*% t(H);
		this.loglik = -(1/2) * (log.det.Sigma + d.Sigma.inv.d);
		
		if(verbose >= 5) cat("  User ",i,": ",this.loglik,"\n",sep="");
		
		ans = ans + drop(this.loglik);
	}
	return(ans);
}
# slower implementation for debugging
marginal.loglik.R2 <- function(feature, response, param, verbose=0){
	
	check.syntax.obs(feature=feature, response=response);
	size = check.syntax.param(param);
	
	obsIndex.f = tapply(seq_len(nrow(feature)),  list(feature$user),  FUN=c, simplify=FALSE);
	obsIndex.r = tapply(seq_len(nrow(response)), list(response$user), FUN=c, simplify=FALSE);
	
	nUsers = max(feature$user, response$user);
	A = NULL;  var_z = NULL;
	for(k in seq_len(size$nApps)){
		if(length(param$A[[k]]) == 1 && param$A[[k]] == 1) A_k = param$A[[k]] * diag(size$nGlobalFactors)
		else                                               A_k = param$A[[k]];
		A = rbind(A, A_k);
		var_z = c(var_z, rep(param$var_z[k], nrow(A_k)));
	}
	AA = A %*% t(A);
	
	ans = 0;
	
	for(i in seq_len(nUsers)){
		this.user = as.character(i); # IMPORTANT: use string to access obsIndex.f and obsIndex.r
		f = feature[ obsIndex.f[[this.user]],];
		r = response[obsIndex.r[[this.user]],];
		if(nrow(f)+nrow(r) == 0) next;
		f.by.app = split(f, f$app);
		r.by.app = split(r, r$app);
		d = rep(NA, nrow(f)+nrow(r));
		D = rep(NA, length(d));
		C = matrix(0, nrow=length(d), ncol=nrow(A));
		num = 0; C_col = 0;
		for(k in seq_len(size$nApps)){
			this.app = as.character(k);
			this.f = f.by.app[[this.app]];
			this.r = r.by.app[[this.app]];
			d_xik = NULL; C_xik = NULL; D_xik = NULL;
			d_yik = NULL; C_yik = NULL; D_yik = NULL;
			if(!is.null(this.f)){
				if(length(param$B[[k]]) == 1 && param$B[[k]] == 1) B = param$B[[k]] * diag(size$nLocalFactors[k])
				else                                               B = param$B[[k]];
				d_xik = this.f$x - param$b[[k]][this.f$index];
				C_xik = B[this.f$index,,drop=FALSE];
				if(is.null(this.f$w)) w = rep(1, nrow(this.f))
				else                  w = this.f$w;
				D_xik = w * param$var_x[k];
			} 
			if(!is.null(this.r)){
				d_yik = this.r$y - param$alpha[[k]][this.r$item];
				C_yik = param$beta[[k]][this.r$item,,drop=FALSE];
				if(is.null(this.r$w)) w = rep(1, nrow(this.r))
				else                  w = this.r$w;
				D_yik = w * param$var_y[k];
			}
			d_ik = c(d_xik, d_yik);
			D_ik = c(D_xik, D_yik);
			C_ik = rbind(C_xik, C_yik);
			if(length(d_ik) > 0){
				d[num+(seq_len(length(d_ik)))] = d_ik;
				D[num+(seq_len(length(D_ik)))] = D_ik;
				C[num+(seq_len(nrow(C_ik))), C_col+(seq_len(ncol(C_ik)))] = C_ik;
			}
			
			num = num + length(d_ik);
			C_col = C_col + size$nLocalFactors[k];
		}
		Sigma = param$var_u * C %*% AA %*% t(C) + C %*% diag(var_z,nrow=length(var_z)) %*% t(C) + 
				if(length(D) == 1) D else diag(D);
		det.Sigma = det(Sigma);
		if(is.na(det.Sigma) || abs(det.Sigma) < 1e-16){
			cat("\nWARNING: (i = ",i,"): det(Sigma) = ",det.Sigma,"\n",sep="");
			cat("Sigma=\n"); print(Sigma);
			cat("param$var_u=\n"); print(param$var_u);
			cat("C=\n"); print(C);
			cat("AA=\n"); print(AA);
			cat("var_z=\n"); print(var_z);
			cat("D=\n"); print(D);
			warning("(i = ",i,"): det(Sigma) = ",det.Sigma,"\n",sep="");
		}
		this.loglik = -(1/2) * (log(det.Sigma) + t(d) %*% solve(Sigma) %*% d);
		
		if(verbose >= 5) cat("  User ",i,": ",this.loglik,"\n",sep="");
		
		ans = ans + drop(this.loglik);
	}
	return(ans);
}

###
### Observation log-likelihood (excluding the factor part)
###
obs.loglik <- function(feature, response, param, factor){
	pred.y = predict.y.from.z(response=response, param=param, z=factor$z, add.noise=FALSE);
	pred.x = predict.x.from.z( feature=feature,  param=param, z=factor$z, add.noise=FALSE);
	if(is.null(response$w)) y.w = rep(1, nrow(response))
	else                    y.w = response$w;
	if(is.null(feature$w))  x.w = rep(1, nrow(feature))
	else                    x.w = feature$w;
	y.var = param$var_y[response$app] * y.w;
	x.var = param$var_x[feature$app]  * x.w;
	if(param$is.x.logistic) x = feature$response
	else                    x = feature$x;
	if(param$is.y.logistic) y = response$response
	else                    y = response$y;
	loglik.x = obsLoglik.from.gaussian(pred.x, x, x.var, param$is.x.logistic, weight=1);
	loglik.y = obsLoglik.from.gaussian(pred.y, y, y.var, param$is.y.logistic, weight=1);
	out = list(total=loglik.x+loglik.y, x=loglik.x, y=loglik.y);
	return(out);
}

###
### Complete data log likelihood
###	NOTE: z_{ik} is included even if user i has no data in app k
###
complete.data.loglik <- function(feature, response, param, factor){
	size = check.syntax.all(feature, response, param, factor);
	temp = obs.loglik(feature=feature, response=response, param=param, factor=factor);
	loglik.x = temp$x;
	loglik.y = temp$y;
	loglik.z = 0;
	for(k in seq_len(size$nApps)){
		if(length(param$A[[k]]) == 1 && param$A[[k]] == 1) Au = factor$u * drop(param$A[[k]])
		else                                               Au = factor$u %*% t(param$A[[k]]);
		if(any(dim(Au) != dim(factor$z[[k]]))) stop("dim(Au) != dim(factor$z[[k]])");
		temp = loglik.gaussian(Au, factor$z[[k]], param$var_z[k]);
		loglik.z = loglik.z + temp;
	}
	zero = matrix(0, nrow=nrow(factor$u), ncol=ncol(factor$u));
	loglik.u = loglik.gaussian(zero, factor$u, param$var_u);
	out = list(total=loglik.x+loglik.y+loglik.z+loglik.u,
			   x=loglik.x, y=loglik.y, z=loglik.z, u=loglik.u);
	return(out);
}

##############################################################################
#   Syntax Sanity Checks
##############################################################################

check.syntax.obs <- function(feature, response){
	check_names(feature,  "feature",  required=c("user", "app", "index", "x"), optional=c("w", "response"));
	check_names(response, "response", required=c("user", "app", "item",  "y"), optional=c("w", "response"));
	if(!is.data.frame(feature))  stop("feature is not a data.frame");
	if(!is.data.frame(response)) stop("response is not a data.frame");
}

check.syntax.param <- function(param){
	check_names(param, "param", required=c("A", "B", "b", "alpha", "beta", "var_x", "var_y", "var_z", "var_u"), optional=c("x.xi", "y.xi", "is.x.logistic", "is.y.logistic"));
	nApps = length(param$A);
	if(length(param$B) != nApps) stop("length(B) != length(param$A)");
	if(length(param$b) != nApps) stop("length(b) != length(param$A)");
	if(length(param$alpha) != nApps) stop("length(alpha) != length(param$A)");
	if(length(param$beta)  != nApps) stop("length(beta) != length(param$A)");
	if(length(param$var_x) != nApps) stop("length(var_x) != length(param$A)");
	if(length(param$var_y) != nApps) stop("length(var_y) != length(param$A)");
	if(length(param$var_z) != nApps) stop("length(var_z) != length(param$A)");
	nGlobalFactors = if(length(param$A[[1]]) != 1) ncol(param$A[[1]]) else ncol(param$beta[[1]]);
	nLocalFactors  = rep(NA, nApps);
	nFeatures      = rep(NA, nApps);
	nItems         = rep(NA, nApps);
	for(k in seq_len(nApps)){
		if(length(param$A[[k]]) == 1 && param$A[[k]] == 1){
			nLocalFactors[k] = nGlobalFactors;
		}else{
			nLocalFactors[k] = nrow(param$A[[k]]);
			if(ncol(param$A[[k]]) != nGlobalFactors) stop("ncol(param$A[[k]]) != nGlobalFactors");
		}
		if(length(param$B[[k]]) == 1 && param$B[[k]] == 1){
			nFeatures[k] = nLocalFactors[k];
		}else{
			nFeatures[k] = nrow(param$B[[k]]);
			if(ncol(param$B[[k]]) != nLocalFactors[k]) stop("ncol(param$B[[k]]) != nLocalFactors[k]");
		}
		if(length(param$b[[k]]) != nFeatures[k]) stop("length(b[[k]]) != nFeatures[k]");
		nItems[k] = nrow(param$beta[[k]]);
		if(ncol(param$beta[[k]]) != nLocalFactors[k]) stop("ncol(param$beta[[k]]) != nLocalFactors[k]");
		if(length(param$alpha[[k]]) != nItems[k]) stop("length(alpha[[k]]) != nItems[k]");
	}
	if(length(param$var_u) != 1) stop("length(param$var_u) != 1");
	out = list(nApps=as.integer(nApps), nGlobalFactors=as.integer(nGlobalFactors), nLocalFactors=as.integer(nLocalFactors), nFeatures=as.integer(nFeatures), nItems=as.integer(nItems));
	return(out);
}

check.syntax.factor <- function(factor){
	check_names(factor, "factor", required=c("u", "z"), optional=c("var.x.score", "var.y.score"));
	nUsers = nrow(factor$u);
	nGlobalFactors = ncol(factor$u);
	nApps = length(factor$z);
	nLocalFactors  = rep(NA, nApps);
	for(k in seq_len(nApps)){
		nLocalFactors[k] = ncol(factor$z[[k]]);
		if(nrow(factor$z[[k]]) != nUsers) stop("nrow(factor$z[[k]]) != nUsers");
	}
	out = list(nApps=as.integer(nApps), nGlobalFactors=as.integer(nGlobalFactors), nLocalFactors=as.integer(nLocalFactors), nUsers=as.integer(nUsers));
}

# return the dimensionality of the problem (#applications, #users, #items, ...)
check.syntax.all <- function(feature, response, param, factor, check.indices=FALSE, test.data=FALSE){
	check.syntax.obs(feature=feature, response=response);
	size.p = check.syntax.param(param=param);
	size.f = check.syntax.factor(factor=factor);
	if(size.p$nApps != size.f$nApps) stop("#applications mismatch");
	if(size.p$nGlobalFactors != size.f$nGlobalFactors) stop("#global factors mismatch");
	if(any(size.p$nLocalFactors != size.f$nLocalFactors)) stop("#local factors mismatch");
	size = list(
		nApps=size.f$nApps, nUsers=size.f$nUsers, nItems=size.p$nItems, nFeatures=size.p$nFeatures,
		nGlobalFactors=size.f$nGlobalFactors, nLocalFactors=size.f$nLocalFactors
	);
	if(check.indices){
		if(any(is.na(feature$user))) stop("some feature$user is NA");
		if(any(feature$user < 1)) stop("some feature$user < 1");
		if(!test.data && any(feature$user > size$nUsers)) stop("some feature$user > nUsers");
		if(any(is.na(feature$app))) stop("some feature$app is NA");
		if(any(feature$app < 1)) stop("some feature$app < 1");
		if(any(feature$app > size$nApps)) stop("some feature$app > nApps");
		if(any(is.na(feature$index))) stop("some feature$index is NA");
		if(any(feature$index < 1)) stop("some feature$index < 1");
		if(any(is.na(response$user))) stop("some response$user is NA");
		if(any(response$user < 1)) stop("some response$user < 1");
		if(!test.data && any(response$user > size$nUsers)) stop("some response$user > nUsers");
		if(any(is.na(response$app))) stop("some response$app is NA");
		if(any(response$app < 1)) stop("some response$app < 1");
		if(any(response$app > size$nApps)) stop("some response$app > nApps");
		if(any(is.na(response$item))) stop("some response$item is NA");
		if(any(response$item < 1)) stop("some response$item < 1");
		
		temp = tapply(seq_len(length(feature$app)), list(feature$app), FUN=c, simplify=FALSE);
		for(i in seq_len(length(temp))){
			k = as.integer(names(temp)[i]);
			select = temp[[i]];
			if(any(feature[select,"index"] > size$nFeatures[k])) stop("some feature$index > nFeatures");
		}
		
		temp = tapply(seq_len(length(response$app)), list(response$app), FUN=c, simplify=FALSE);
		for(i in seq_len(length(temp))){
			k = as.integer(names(temp)[i]);
			select = temp[[i]];
			if(any(response[select,"item"] > size$nItems[k])) stop("some response$item > nItems");
		}
	}
	return(size);
}

##############################################################################
#   Output results
##############################################################################

output.results <- function(
	method, feature, response, model, prediction,
	out.dir, out.level, 
	minTestLoss, iter, show.marginal.loglik,
	TimeEM, TimeTest, verbose,
	other=NULL, name="model"
){
	if(!is.null(prediction) && is.null(prediction$loss.x)) stop("prediction$loss.x does not exist");
	if(!is.null(prediction) && is.null(prediction$loss.y)) stop("prediction$loss.y does not exist");
	
	marginal.loglik = NA;
	timeLoglik      = NA;
	timeMgnLoglik   = NA;
	
	if(verbose > 0 || out.level > 0){
		begin.time = proc.time();
		cd.loglik = complete.data.loglik(feature=feature, response=response, param=model$param, factor=model$factor);
		timeLoglik = proc.time() - begin.time;
		if(verbose > 0) 
			cat(sep="",
			"Complete data log likelihood: ",cd.loglik$total,"   (using ",timeLoglik[3]," sec)\n");
		if(verbose >= 2){
			cat(sep="",
			"                    response: ",cd.loglik$y,"\n",
			"                     feature: ",cd.loglik$x,"\n",
			"                           z: ",cd.loglik$z,"\n",
			"                           u: ",cd.loglik$u,"\n");
		}
		if(show.marginal.loglik){
			cat(sep="",
			"     Marginal log likelihood: ");
			begin.time = proc.time();
			marginal.loglik = marginal.loglik.R(feature=feature, response=response, param=model$param);
			timeMgnLoglik = proc.time() - begin.time;
			cat(marginal.loglik,"   (using ",timeMgnLoglik[3]," sec)\n",sep="");
		}
		if(verbose > 0 && !is.null(prediction)){
			cat(sep="",
			"      Test-set response loss: ",prediction$loss.y,"\n",
			"                feature loss: ",prediction$loss.x,"\n");
		}
	}
	
	out = list(marginal.loglik=marginal.loglik);
	
	if(out.level <= 0) return(out);
	
	b.time.write = proc.time();
	
	thisTestLoss   = if(is.null(prediction)) NA else prediction$loss.y;
	thisTestLoss.x = if(is.null(prediction)) NA else prediction$loss.x;
	
	if(iter == 0){
		if(out.level >= 2) save(file=paste(out.dir,"/",name,".0",sep=""), list=c("model"));
	}else{
		file = paste(out.dir,"/",name,".last",sep="");
		if(out.level >= 2){
			file.prev = paste(out.dir,"/",name,".",(iter-1),sep="");
			if(file.exists(file.prev)) file.remove(file.prev);
			file.rename(file, file.prev);
		}
		save(file=file, list=c("model", "prediction", "other"));
		if(!is.null(prediction)){
			if(thisTestLoss == minTestLoss) file.copy(file, paste(out.dir,"/",name,".minTestLoss",sep=""), overwrite=TRUE);
		}
	}
	
	summary = data.frame(
			Method=method, Iter=iter, 
			TestLoss.y=thisTestLoss, TestLoss.x=thisTestLoss.x, MgnLoglik=marginal.loglik, CDLoglik=cd.loglik$total, 
			TimeEM=TimeEM[3], TimeTest=TimeTest[3], TimeLoglik=timeLoglik[3], TimeMgnLoglik=timeMgnLoglik[3]
	);
	file = paste(out.dir,"/summary.txt",sep="");
	if(file.exists(file) && iter>0) write.table(summary, file=file, append=TRUE,  quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	else                            write.table(summary, file=file, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE);
	if(verbose > 0){
		time.used.write = proc.time() - b.time.write;
		cat("Write a model & summary info onto disk (using ",time.used.write[3]," sec)\n",sep="");
	}
	return(out);
}

