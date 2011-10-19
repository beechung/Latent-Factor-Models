### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

rgamma.psize <- function(n, mean, size){
	scale = 1/size;
	shape = mean * size;
	return(rgamma(n, shape=shape, scale=scale));
}

rgamma.mean.var <- function(n, mean, var){
	size = mean/var;
	scale = 1/size;
	shape = mean * size;
	return(rgamma(n, shape=shape, scale=scale));
}

logit <- function(p){
	return(log(p/(1-p)));
}

inv.logit <- function(x){
	return(1/(1+exp(-x)));
}

check_type_size <- function(x, type, size, isNullOK=FALSE, check.NA=TRUE, name=NULL){
	if(is.null(name)) name = "The input"
	if(is.null(x)){
		if(any(size == 0) || isNullOK) return(TRUE)
		else stop(name," is null");
	}
	if(type == "double"){
		if(!is.double(x)) stop(name," should be double");
	}else if(type == "integer" || type == "int"){
		if(!is.integer(x)) stop(name," should be integer");
	}else if(!is.na(type)) stop("Unknown type: ",type,sep="");
	
	if(length(size) > 1 || !all(is.na(size))){
		d = dim(x);
		if(is.null(d)) d = length(x);
		if(length(d) != length(size) || any(d != size)) stop(name," has dimensionality mismatch: (",paste(d,collapse=" x "),") vs (",paste(size,collapse=" x "),")");
	}
	if(check.NA && any(is.na(x))) stop(name," has some elements that are NA");
}

warning.any.not.in <- function(a, b, msg, stop=FALSE, print=FALSE){
	if(any(!(a %in% b))){
		temp = a[!(a %in% b)];
		if(stop) stop(msg,paste(temp, collapse=", "))
		else     warning(msg,paste(temp, collapse=", "),call.=FALSE);
		if(print) cat("\nWARNING: ",msg,paste(temp, collapse=", "),"\n\n",sep="");
	}
}

###
### Functions for the logistic model
### 
### Initialize the observation table
###   Observation table: obs
###	  Add the following to the code before the EM procedure.
###      if(is.logistic) obs = init.obs(obs, target="y");
###   After the call, obs$response will contain the original binary response.
###   Note obs$y will be used to store the Gaussian response values
###   for variational approximation
init.obs.logistic <- function(obs, target="y"){
	if(is.null(obs)) return(NULL);
	obs$response = obs[[target]]; # Make a copy of the original response
	                              # since obs$y may be changed for variational approx.
	# sanity check
	labels = unique(obs$response);
	if(length(labels) != 2) stop("The response is not binary: ",paste(labels[1:min(10,length(labels))],collapse=", "));
	labels = sort(labels);
	if(any(labels != c(0,1)) && any(labels != c(-1,1))) stop("Binary response must be {0,1} or {-1,1}");
	return(obs);
}
init.param.logistic <- function(param, obs, value=1.0, target="y"){
	if(is.null(obs[[target]])) stop("obs$",target," is null");
	param[[paste(target,".xi",sep="")]] = rep(value, nrow(obs));
	return(param);
}
### Generate Gaussain response
### Add the following code in the EM loop before the E-step
### if(is.logistic){
###    obs$y  = get.response.logistic(obs, param, target="y", verbose=verbose);
###    param$var_y = get.var.logistic(obs, param, target="y", verbose=verbose);
### }
get.var.logistic <- function(obs, param, target="y", verbose=0){
	# variational approximation
	if(verbose >= 1) cat("generate gaussian var for logistic\n");
	xi = param[[paste(target,".xi",sep="")]];
	if(length(xi) != nrow(obs)) stop("length(xi) != nObs");
	var_y = 1/(2 * logistic.lambda(xi));
	return(var_y);
}
get.response.logistic <- function(obs, param, target="y", verbose=0){
	if(verbose >= 1) cat("generate gaussian response for logistic\n");
	var_y = get.var.logistic(obs=obs, param=param, target=target, verbose=0);
	if(length(obs$response) != nrow(obs)) stop("length(obs$response) != nObs");
	response = obs$response;
	if(all(response %in% c(0,1))) response[response == 0] = -1;
	if(!all(response %in% c(-1,1))) stop("Binary response must be {-1,1} at this point");
	y = response * var_y / 2;
	return(y);
}
### Update the parameters
### Add the following code in the EM loop after the M-step
### if(is.logistic){
###    param = update.param.logistic(param, mean.score, var.score, target="y");
### }
update.param.logistic <- function(param, mean.score, var.score, target="y"){
	if(length(mean.score) != length(var.score)) stop("length(mean.score) != length(var.score)");
	param[[paste(target,".xi",sep="")]] = sqrt(mean.score^2 + var.score);
	return(param);
}
# input pred.y is based on the Gaussian model
obsLoglik.from.gaussian <- function(pred.y, y, var_y, is.logistic, weight){
	if(length(y) != length(pred.y)) stop("length(y) != length(pred.y)");
	if(!is.null(is.logistic) && is.logistic){
		if(all(y %in% c(0,1))) y[y == 0] = -1;
		if(!all(y %in% c(-1,1))) stop("Binary response must be {-1,1} at this point");
		loglik = sum( -log1p(exp(-y * pred.y) ) );
		attr(loglik, "loss") = -loglik / length(y);
	}else{
		if(is.null(weight)) weight = 1;
		loglik = loglik.gaussian(pred.x=pred.y, x=y, var_x=var_y);
		attr(loglik, "loss") = sqrt(mean( weight*(y - pred.y)^2 ));
	}
	return(loglik);
}
# input pred.y is based on the Gaussian model
# output$pred.y is the input for the Gaussian model
#               is the predicted probability for the Logistic model
predict.response.from.gaussian <- function(pred.y, y, var_y, is.logistic, weight){
	if(!is.null(y)) loglik = obsLoglik.from.gaussian(pred.y=pred.y, y=y, var_y=var_y, is.logistic=is.logistic, weight=weight)
	else            loglik = NA;
	if(!is.null(is.logistic) && is.logistic){
		pred.y = 1/(1+exp(-pred.y));
		if(!is.null(y)){
			if(min(y) == -1) temp.y = 2*pred.y - 1
			else             temp.y = pred.y;
		}
	}else{
		temp.y = pred.y;
	}
	if(!is.null(y)){
		if(is.null(weight)) weight = 1;
		rmse = sqrt(mean( weight*(y - temp.y)^2 ));
		mae  = mean( abs( weight*(y - temp.y) ) );
	}else{
		rmse = NA;
		mae  = NA;
	}
	return(list(pred.y=pred.y, true.y=y, rmse=rmse, mae=mae, loglik=loglik, test.loss=attr(loglik, "loss")));
}
logistic.lambda <- function(xi){
	return(tanh(xi/2) / (4*xi));
}

loglik.gaussian <- function(pred.x, x, var_x){
	if(length(pred.x) != length(x)) stop("Input data lengths mismatch");
	if(length(var_x) == 1){
		loglik = -(1/2) * ( sum((x - pred.x)^2 / var_x) + length(x) * log(var_x) );
	}else if(length(var_x) == length(x)){
		loglik = -(1/2) * sum( (x - pred.x)^2 / var_x + log(var_x) );
	}else if(is.matrix(var_x)){
		if(nrow(var_x) != length(x) || ncol(var_x) != length(x)) stop("dim(var_x) = ",paste(dim(var_x),collapse=" x "));
		diff = x - pred.x;
		loglik = -(1/2) * ( log(det(var_x)) + t(diff) %*% solve(var_x) %*% diff );
	}else stop("length(var_x) is wrong");
	return(loglik);
}

RMSE <- function(x, y, weight=1){
	return(sqrt(mean( weight*((x-y)^2) )));
}

subsample.ROC <- function(perf, nPoints=1000){
	n = length(perf@x.values[[1]]);
	if(n <= nPoints) return(perf);
	size = floor(n / nPoints);
	select = ((1:n) %% size) == 1;
	select[n] = TRUE;
	perf@x.values[[1]] = perf@x.values[[1]][select];
	perf@y.values[[1]] = perf@y.values[[1]][select];
	perf@alpha.values[[1]] = perf@alpha.values[[1]][select];
	return(perf);
}

# Set optional = NULL to disable checking for unused components
# Missing required    => stop
# Missing should.have => warn
check_names <- function(x, display.name, required, should.have=c(), optional=c()){
	names.all = c(required, optional, should.have);
	if(!is.null(optional)) warning.any.not.in(names(x), names.all, paste("The following components should not be in ",display.name,": ", sep=""), stop=TRUE);
	warning.any.not.in(required,  names(x),  paste("You must have the following components in ",display.name,": ",sep=""), stop=TRUE);
	warning.any.not.in(should.have, names(x),paste("The following components are missing in ",display.name,": ",sep=""), stop=FALSE, print=TRUE);
}

###
### Create a deep copy of the input object
###
deepCopy <- function(x){
	if(is.list(x)){
		out = list();
		for(i in seq_len(length(x))){
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
	}else if(is.logical(x)){
		out = x & TRUE;
	}else if(is.null(x)){
		out = NULL;
	}else stop("Type not supported");
	
	return(out);
}

memcpy <- function(dst, src, start_dst=1, start_src=1, len=length(dst)){
	if(len == 0) return(NULL);
	if(is.double(dst)){
		if(!is.double(src)) stop("type mismatch");
		.C("copy_double_array", 
			dst, src, as.integer(length(dst)), as.integer(length(src)), 
			as.integer(start_dst), as.integer(start_src), as.integer(len), DUP=FALSE
		);
	}else if(is.integer(dst)){
		if(!is.integer(src)) stop("type mismatch");
		.C("copy_int_array", 
			dst, src, as.integer(length(dst)), as.integer(length(src)), 
			as.integer(start_dst), as.integer(start_src), as.integer(len), DUP=FALSE
		);
	}else stop("Unknown type");
	return(NULL);
}

rep_matrix <- function(matrix, num){
	ans = array(NA, dim=c(num, nrow(matrix), ncol(matrix)));
	for(i in seq_len(num)){
		ans[i,,] = matrix;
	}
	return(ans);
}

rep_3DA <- function(array_3D, num){
	if(length(dim(array_3D))!=3) stop("length(dim(array_3D))!=3");
	ans = array(NA, dim=c(num, dim(array_3D)));
	for(i in seq_len(num)){
		ans[i,,,] = array_3D;
	}
	return(ans);
}

scalar_times_I <- function(s, nrow){
	ans = array(NA, dim=c(length(s), nrow, nrow));
	for(i in seq_len(length(s))){
		ans[i,,] = s[i] * diag(nrow);
	}
	return(ans);
}

matrix_times_3DA <- function(matrix, array_3D){
	if(length(dim(array_3D))!=3) stop("length(dim(array_3D))!=3");
	dim2 = dim(array_3D)[1]; dim3 = dim(array_3D)[3];
	ans = array(NA, dim(nrow(matrix), dim2, dim3));	
	for(i in seq_len(dim2)){
		ans[,i,] = matrix %*% array_3D[i,,];
	}
	return(ans);
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
	dim1 = dim(x1);
	dim2 = dim(x2);
	if(length(dim1) != length(dim2) || any(dim1 != dim2)){
		cat(prefix,": Different dimensions!\n",sep="");
		return(TRUE);
	}
	if(is.list(x1)){
		if(!is.list(x2)){
			cat(prefix,": Different types! (list vs non-list)\n",sep="");
			return(TRUE);
		}
		name1 = names(x1); if(!is.null(name1)) name1 = sort(name1);
		name2 = names(x2); if(!is.null(name2)) name2 = sort(name2);
		if(is.null(name1) || is.null(name2)){
			if(!is.null(name2) || !is.null(name1)){
				cat(prefix,": One has no names; the other has names!\n",sep="");
				return(TRUE);
			}
			ans = FALSE;
			for(i in seq_len(length(x1))){
				ret = is.diff(x1[[i]], x2[[i]], precision, prefix=paste(prefix,"[[",i,"]]",sep=""));
				ans = ans || ret;
			}
			return(ans);
		}else{
			if(length(name1) != length(name2) || any(name1 != name2)){
				cat(prefix,": Different names!\n",sep="");
				return(TRUE);
			}
			ans = FALSE;
			for(i in seq_len(length(name1))){
				name = name1[i];
				ret = is.diff(x1[[name]], x2[[name]], precision, prefix=paste(prefix,"$",name,sep=""));
				ans = ans || ret;
			}
			return(ans);
		}
	}else{
		indexDiff = (seq_len(length(x1)))[x1 != x2];
		# print(indexDiff);
		if(length(indexDiff) == 0) return(FALSE);
		
		value = cbind(x1[indexDiff], x2[indexDiff]);
		denom = apply(abs(value), 1, max);
		diff  = abs(value[,1]-value[,2]) / denom;
		
		# print(diff);
		
		temp = (seq_len(length(indexDiff)))[diff > precision];
		indexDiff = indexDiff[temp];
		diff      = diff[temp];
		if(length(indexDiff) == 0) return(FALSE);
		
		for(i in seq_len(length(indexDiff))){
			index = indexDiff[i];
			if(length(dim1) == 2){
				row = (index-1) %% dim1[1] + 1;
				col = floor((index-1) / dim1[1]) + 1;
				cat(prefix,"[",row,",",col,"]: ",x1[index]," vs ",x2[index]," (relative diff=",diff[i],")\n",sep="");
			}else if(length(dim1) == 3){
				d1 = (index-1) %% dim1[1] + 1;
				d2 = (floor((index-1) / dim1[1])) %% dim1[2] + 1;
				d3 = floor((index-1) / (dim1[1]*dim1[2])) + 1;
				cat(prefix,"[",d1,",",d2,",",d3,"]: ",x1[index]," vs ",x2[index]," (relative diff=",diff[i],")\n",sep="");
			}else{
				cat(prefix,"[",index,"]: ",x1[index]," vs ",x2[index]," (relative diff=",diff[i],")\n",sep="");
			}
		}
		return(TRUE);
	}
}

###
### Join the columns of the input data frame into a vector of strings
### (one joined string per row)
###
join.columns = function(data, sep="\t"){
	if(!is.data.frame(data)) stop("Input must be a data frame!");
	output = data[,1];
	for(i in 2:ncol(data)){
		temp = data[,i];
		output = paste(output, temp, sep=sep);
	}
	return(output);
}

###
### split.columns(data, "source.column", c("dest.column1", "dest.column2"))
### split(data[,"source.column"]) and put the results in data[,c("dest.column1", "dest.column2")]
###
split.columns = function(data, src, dest, sep="\t", ...){
	if(!is.data.frame(data)) stop("Input must be a data frame!");
	if(!(src %in% names(data))) stop("Column ",src," is not in the input data");
	temp = strsplit(as.character(data[,src]), sep, ...);
	len = sapply(temp, length);
	if(any(len != length(dest))) stop("Not all rows have the same number of columns");
	result = matrix(unlist(temp), nrow=nrow(data), ncol=length(dest), byrow=TRUE);
	for(i in seq_len(length(dest))){
		data[,dest[i]] = result[,i];
	}
	return(data);
}

### Output a table grouped by the group.by columns
### Output schema: (group.by[1], ..., group.by[n], fn[[1]](attr[1]), ..., fn[[k]](attr[k]))
### E.g., table.aggregate(table, c("city", "state"), c("value", "value"), c(mean, var));
### Inefficient if group.by has multiple attributes (because of tapply in aggregate)!!
###
table.aggregate = function(table, group.by, attr, fn, use.joined.column=TRUE, ...){
	
	attrnames = names(table);
	if(!all(group.by %in% attrnames)) stop("Some group.by attributes (",group.by,") are not in the table");
	if(!all(attr %in% attrnames)) stop("Some group.by attributes (",attr,") are not in the table");
	if(length(group.by) == 1) use.joined.column=FALSE;
	
	by = table[,group.by];
	if(use.joined.column) by = join.columns(by);
	if(!is.list(by)) by = list(by);
	
	if(length(fn) == 1){
		output = aggregate(table[,attr], by, fn, ...);
		names(output)[length(by)+(seq_len(length(attr)))] = attr;
	}else{
		if(length(attr) != length(fn)) stop("length(attr) should be equal to length(fn)!");
		output = aggregate(table[,attr[1]], by, fn[[1]], ...);
		names(output)[length(by)+1] = attr[1];
		for(i in 2:length(attr)){
			temp = aggregate(table[,attr[i]], by, fn[[i]], ...);
			names(temp)[length(by)+1] = attr[i];
			output = merge(output, temp);
		}
	}
	
	if(use.joined.column){
		output = split.columns(output, names(output)[1], group.by);
		output = output[,c(group.by, attr)];
		for(i in seq_len(length(group.by))){
			if(     is.factor( table[1,group.by[i]])) output[,i] = factor(output[,i], levels=levels(table[1,group.by[i]]))
			else if(is.integer(table[1,group.by[i]])) output[,i] = as.integer(output[,i])
			else if(is.numeric(table[1,group.by[i]])) output[,i] = as.numeric(output[,i]);
		}
	}else{
		names(output) = c(group.by, attr);
	}
	
	return(output);
}


append.to.file <- function(..., file.name=NULL, quote=FALSE, sep="\t",
		eol="\n", na="NA", dec=".", row.names=FALSE, qmethod="escape")
{
	data = data.frame(...);
	if(is.null(file.name)) stop("file.name cannot be null!");
	if(file.exists(file.name)){
		append = TRUE;
		col.names = FALSE;
	}else{
		append = FALSE;
		col.names = TRUE;
	}
	write.table(data,file=file.name, append=append, col.names=col.names,
			quote=quote,sep=sep,eol=eol,na=na,dec=dec,row.names=row.names,qmethod=qmethod);
}

###
### x is a vector, matrix or 3D array of boolean values
### This function returns the indices of x==TRUE
###
get.index <- function(x){
	if(!is.logical(x)) stop("The input must be logical (TRUE/FALSE)!");
	if(is.vector(x)){
		return((1:length(x))[x]);
	}else if(is.matrix(x)){
		z = cbind(rep(1:nrow(x),times=ncol(x)), rep(1:ncol(x),each=nrow(x)));
		return(z[x,,drop=FALSE]);
	}else if(is.array(x) && length(dim(x)) == 3){
		d = dim(x);
		z = cbind(rep(1:d[1],times=d[2]*d[3]), rep(rep(1:d[2], each=d[1]),times=d[3]), rep(1:d[3],each=d[1]*d[2]));
		return(z[x,,drop=FALSE]);
	}else stop("The function only support vector, matrix and 3D array!");
}

R.square <- function(y.obs, y.pred){
	SS.tot = sum((y.obs - mean(y.obs))^2);
	SS.err = sum((y.obs - y.pred)^2);
	return(1 - SS.err/SS.tot);
}

timing <- function(func, times=1){
	fn <- function(){
		for(i in 1:times){
			func();
		}
	}
	out = system.time(fn());
	return(out);
}


###
### Plot (x[i], y[i]) grouped by the group.by attributes and ordered by the order.by attributes
###
plot.groupBy = function(
		x, y, group.by, order.by, type="b", x.lines=NULL, y.lines=NULL, x.lines.fine=FALSE, attach=FALSE,
		pch.start=NULL, pch=NULL, lty=NULL, col=NULL, cex.eachPoint=NULL, cex=1, ylim=NULL, 
		legend.pos=NULL, legend.bg=NULL, legend.cex=1, legend.box.lty=1, las=1,...
){
	temp = y[!is.na(y) & y != Inf & y!= -Inf];
	if(is.null(ylim)) ylim=c(min(temp), max(temp));
	if(!attach){
		if(is.factor(x)){
			plot(x, rep(NA, length(y)), ylim=ylim, xaxt="n", ...);
			axis(1, at=1:length(levels(x)), labels=levels(x), las=las);
			x.all = 1:nlevels(x);
		}else{
			plot(x, rep(NA, length(y)), ylim=ylim, las=las, ...);
			x.all = unique(x);
		}
		
		if(x.lines.fine){
			axis(1, x.all, labels=FALSE, tck=1, col="#DDDDDD");
		}
		if(!is.null(x.lines)){
			axis(1, x.lines, labels=FALSE, tck=1, col="#BBBBBB");
		}
		if(!is.null(y.lines)){
			axis(2, y.lines, labels=FALSE, tck=1, col="#CCCCCC");
		}
	}
	
	index = tapply(1:length(y), group.by, c, simplify=F);
	
	if(is.null(col)){
		col = rep(c(1,2,3,4,6,8), ceiling(length(index)/6)); 
	}else if(!is.null(names(col))){
		col = col[names(index)];
		if(any(is.na(col))) stop("col should have names for all group.by");
	}else{
		col = rep(col, ceiling(length(index)/length(col)));
	}
	if(is.null(lty)){ 
		lty = rep(1, length(index));
	}else if(!is.null(names(lty))){
		lty = lty[names(index)];
		if(any(is.na(lty))) stop("lty should have names for all group.by");
	}else{
		lty = rep(lty, ceiling(length(index)/length(lty)));
	}
	if(is.null(pch)){
		pch = 1:length(index);
	}else if(!is.null(names(pch))){
		pch = pch[names(index)];
		if(any(is.na(pch))) stop("pch should have names for all group.by");
	}else{
		pch = rep(pch, ceiling(length(index)/length(pch)));
	}
	
	for(k in 1:length(index)){
		select = index[[k]];
		if(is.null(select)) next;
		x.s = x[select];
		y.s = y[select];
		
		if(!is.null(cex.eachPoint)) cex1 = cex.eachPoint[select]
		else cex1 = cex;
		
		if(!is.null(order.by)){
			by  = order.by[select];
			ord = order(by);
			x.s = x.s[ord];
			y.s = y.s[ord];
			if(!is.null(cex.eachPoint)) cex1 = cex1[ord]
		}
		lines(x.s, y.s, type=type, col=col[k], lty=lty[k], pch=pch[k], cex=cex1);
		if(!is.null(pch.start)){
			points(x.s[1], y.s[1], pch=pch.start, col=col[k], cex=cex1[1]);
		}
	}
	if(!is.null(legend.pos)){
		if(is.null(legend.bg))
			legend(legend.pos, legend=names(index), col=col, lty=lty, pch=pch, cex=legend.cex, box.lty=legend.box.lty)
		else
			legend(legend.pos, legend=names(index), col=col, lty=lty, pch=pch, bg=legend.bg, cex=legend.cex, box.lty=legend.box.lty);
	}
}


###
### Convert a matrix into the (row,col,value) format
###
matrix.to.index.value <- function(x){
	if(!is.matrix(x)) stop("The input is not a matrix");
	select = (!is.na(x)) & (x != 0);
	row.index = rep(1:nrow(x), times=ncol(x));
	col.index = rep(1:ncol(x), each =nrow(x));
	out = data.frame(row=row.index[select], col=col.index[select], value=x[select]);
	return(out);
}

###
### Convert (row,col,value) format to a matrix
###
index.value.to.matrix <- function(x, nrow=NULL, ncol=NULL){
	if(!is.data.frame(x)) stop("The input is not a data frame");
	if(!all(c("row", "col", "value") %in% names(x))) stop("The input data frame must have the following columns: row, col, value");
	if(is.null(nrow)) nrow = max(x$row);
	if(is.null(ncol)) ncol = max(x$col);
	out = matrix(0.0, nrow=nrow, ncol=ncol);
	index = cbind(x$row, x$col);
	out[index] = x$value;
	return(out);
}


###
### Pack/unpack list-of-matrices
### Input: list(M1, M2, ...)
### Output:
###		data = a vector that concat all matrices
###     dim  = a vector of the following form
###            #matrices, start(M1), nrow(M1), ncol(M1), start(M2), nrow(M2), ncol(M2), ...
### Note: For sanity check, we append a random number to the end of both data and dim 
###       The two numbers should match
###
### If output != NULL, the result will be filled into the already created output space.
###
pack.list.of.matrices <- function(x, output=NULL){
	if(!is.list(x)) stop("The input is not a list");
	nMatrices = length(x);
	rnd = as.integer(sample.int(n=1e+8, size=1));
	if(!is.null(output)){
		dim  = output$dim;
		if(length(dim) != 2+3*nMatrices) stop("The output$dim has a different size");
		if(!is.integer(dim)) stop("!is.integer(dim)");
	}else{
		dim = as.integer(rep(0, 2+3*nMatrices));
	}
	start = 1;
	memcpy(dim, nMatrices, start_dst=1,           start_src=1, len=1);
	memcpy(dim, rnd,       start_dst=length(dim), start_src=1, len=1);
	for(k in seq_len(nMatrices)){
		if(is.matrix(x[[k]])){      nrow_=nrow(x[[k]]);   ncol_=ncol(x[[k]]);}
		else if(is.vector(x[[k]])){ nrow_=length(x[[k]]); ncol_=1;}
		else stop("x[[",k,"]] is not a matrix or vector");
		if(!is.double(x[[k]])) stop("x[[",k,"]] is not double");
		temp = as.integer(c(start, nrow_, ncol_));
		memcpy(dim, temp, start_dst=((k-1)*3+2), start_src=1, len=3);
		start = start + length(x[[k]]);
	}
	end = start;

	if(!is.null(output)){
		data = output$data;
		if(length(data) != end) stop("The output$data has a different size");
		if(!is.double(data)) stop("!is.double(data)");
	}else{
		data = rep(0.0, end);
	}
	rnd = as.double(rnd);
	memcpy(data, rnd, start_dst=end, start_src=1, len=1);
	start = 1;
	for(k in seq_len(nMatrices)){
		memcpy(data, x[[k]], start_dst=start, start_src=1, len=length(x[[k]]));
		start = start + length(x[[k]]);
	}
	if(is.null(output)) output = list(data=data, dim=dim);
	return(output);
}
unpack.list.of.matrices <- function(x, output=NULL){
	check_names(x, "x", required=c("data", "dim"));
	if(x$data[length(x$data)] != x$dim[length(x$dim)]) stop("sanity check failed");
	nMatrices = x$dim[1];
	if(length(x$dim) != 2+3*nMatrices) stop("length(x$dim) != 2+3*nMatrices");
	if(is.null(output)){
		output = list();
	}else{
		if(length(output) != nMatrices) stop("length(output) != nMatrices");
	}
	size = 0;
	for(k in seq_len(nMatrices)){
		temp = x$dim[((k-1)*3+1) + (1:3)];
		start = temp[1];  nrow = temp[2];  ncol = temp[3];
		if(length(output) < k){
			output[[k]] = matrix(0.0, nrow=nrow, ncol=ncol);
		}else{
			if(is.matrix(output[[k]])){
				check_type_size(output[[k]], "double", c(nrow, ncol), name="output[[k]]");
			}else if(is.vector(output[[k]])){
				if(ncol != 1) stop("output[[k]] is a vector, but ncol = ",ncol);
				check_type_size(output[[k]], "double", nrow*ncol, name="output[[k]]");
			}else stop("output[[k]] is not matrix or a vector");
		}
		memcpy(output[[k]], x$data, start_dst=1, start_src=start, len=nrow*ncol);
		size = size + length(output[[k]]);
	}
	if(size != length(x$data)-1) stop("length(x$data) != the size specified by x$dim: ");
	return(output);
}

###############################################################################
# C Functions (implemented in utils.c)
###############################################################################

###
### Print/get the address of a variable
###
print.address <- function(x){
	if(is.double(x)){
		.C("print_doublePointer", x, DUP=FALSE);
	}else if(is.integer(x)){
		.C("print_intPointer", x, DUP=FALSE);
	}else stop("x is not double or integer")
}
get.address <- function(x){
	if(is.list(x)){
		out = list();
		for(i in seq_len(length(x))){
			out[[i]] = get.address(x[[i]]);
		}
		names(out) = names(x);
		return(out);
	}else if(is.double(x)){
		address = integer(1);
		.C("get_doublePointer", x, address, DUP=FALSE);
		return(address+0);
	}else if(is.integer(x) || is.logical(x)){
		address = integer(1);
		.C("get_intPointer", x, address, DUP=FALSE);
		return(address+0);
	}else stop("x is not double or integer");
}

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
### See generateObsIndex in utils.c
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

indexWithQuantities <- function(x){
	len = sum(x);
	out = integer(len);
	ans = .C("indexWithQuantities",
			out, as.integer(x), as.integer(length(x)),
			DUP=FALSE
	);
	return(out);
}
