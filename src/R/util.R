### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

###
### Replace some rows of x by the rows in y
###   Join by x[,key] and y[,key]
###   then set x[,intersect(names(x),names(y))] = y[,intersect(names(x),names(y))]
###
replace.data.frame <- function(x, y, key){
	if(!all(key %in% names(x))) stop("Join key is not in x");
	if(!all(key %in% names(y))) stop("Join key is not in y");
	if(length(key)==1 && nrow(y) != length(unique(y[,key]))) stop("Join key is not unique in y");
	if(length(key)>=2 && nrow(y) != nrow(unique(y[,key]))) stop("Join key is not unique in y");
	attrs = intersect(names(x),names(y));
	attrs.x1 = c(setdiff(names(x),names(y)), key);
	x1 = x[x[,key] %in% y[,key],attrs.x1];
	x2 = x[!(x[,key] %in% y[,key]),];
	x3 = merge(x1,y)[,names(x)];
	if(nrow(x3) != nrow(x1)) stop("something is wrong");
	out = rbind(x2,x3);
}

###
### Convert a matrix into the (row,col,value) format
###
matrix.to.index.value <- function(x, order.by=NULL){
	if(!is.matrix(x)) stop("The input is not a matrix");
	select = (!is.na(x)) & (x != 0);
	row.index = rep(1:nrow(x), times=ncol(x));
	col.index = rep(1:ncol(x), each =nrow(x));
	out = data.frame(row=row.index[select], col=col.index[select], value=x[select]);
	if(!is.null(order.by)){
		if(order.by == "row"){ out = out[order(out$row, out$col),];} else
		if(order.by == "col"){ out = out[order(out$col, out$row),];}
		else stop("order.by must be either 'row' or 'col'");
	}
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
### Return a dataframe of (by, values[1], ..., values[n])
###
hist.group.by <- function(x, by, values=NULL){
	if(is.null(values)){
		values = unique(x);
	}else{
		if(length(unique(values)) != length(values)) stop("The input values have duplicates!!");
		temp = unique(x);
		if(!all(temp %in% values)) warning("'",paste(temp[!(temp %in% values)], collapse="', '"), "' is/are not in values!!");
		if(!all(values %in% temp)){
			warning("'",paste(values[!(values %in% temp)], collapse="', '"), "' is/are not in x!!");
			values = values[values %in% temp];
		}
	}
	if(is.factor(x)) x = as.character(x);
	data = aggregate(x, by=list(by, x), FUN=length);
	out = data[data[,2] == values[1], c(1,3)];
	names(out) = c("by", as.character(values[1]));
	for(i in 2:length(values)){
		temp = data[data[,2] == values[i], c(1,3)];
		names(temp) = c("by", as.character(values[i]));
		out = merge(out, temp, all=TRUE);
	}
	for(i in 2:(length(values)+1)){
		out[is.na(out[,i]),i] = 0;
	}
	return(out);
}

score.to.rank <- function(score){
	temp = rep(0, length(score));
	temp[order(score, decreasing=T)] = 1:length(score);
	return(temp);
}

editorial.score <- function(data, score=c(excellent=2, good=1, fair=0, bad=-1, abuse=-10)){
	if(!all(c("yuid", "label") %in% names(data))) stop("yuid or label is not a column of the input data");
	data = data[!(data$label %in% c("nj", "")),];
	data$score = score[as.character(data$label)];
	if(any(is.na(data$score))) stop("Some labels are not in the following list: ", paste(names(score),collapse=" ,"));
	temp = data;
	temp$var = data$score;
	temp$nLabels = 1;
	ans = table.aggregate(table=temp, group.by="yuid", attr=c("score", "var", "nLabels"), fn=c(mean, var, length));
	return(ans);
}

###
### Create buzz dataset
###  data.train = data.frame(article_id, comment_time, rating_time, comment_id, voter_id, author_id, rating);
###  data.test  = data.frame(article_id, comment_time, rating_time, comment_id, voter_id, author_id, rating);
###  user.stats = data.frame(yuid, ..., in.U.P, ..., nComments)
###                       in.U.P: Number of unique Users who voted the yuid Positively
###					   nComments: Number of comments that the yuid posted
###  Selection criteria (min.nComments and min.posInDegree) are based on
###  user.stats, instead of the input data
create.buzz.dataset <- function(
	data.train, data.test, user.stats, 
	min.nComments, min.posInDegree,
	out.file=NULL, converTimestamp=TRUE
){
	data = data.train;
	cat("Total number of ratings (in training set): ",nrow(data),"\n",
		"Total number of authors (in training set): ",length(unique(data$author_id)),"\n",
		"Total number of  voters (in training set): ",length(unique(data$voter_id)),"\n", sep="");

	selected.user = user.stat$yuid[user.stat$nComments >= min.nComments & user.stat$in.U.P >= min.posInDegree];
	d.1 = data[data$voter_id %in% selected.user & data$author_id %in% selected.user,];
	out = d.1[d.1$voter_id %in% unique(d.1$author_id) & d.1$author_id %in% unique(d.1$voter_id),];
	
	user.train = unique(c(out$author_id, out$voter_id));
	
	cat("Selected number of ratings (in training set): ",nrow(out),"\n",
		"Selected number of   users (in training set): ",length(user.train),"\n",
		"Selected number of authors (in training set): ",length(unique(out$author_id)),"\n",
		"Selected number of  voters (in training set): ",length(unique(out$voter_id)),"\n",sep="");
	
	if(converTimestamp){
		out$comment_time = str2timestamp(out$comment_time);
		out$rating_time  = str2timestamp(out$rating_time);
	}
	
	out.train = out;
	
	data = data.test;
	cat("Total number of ratings (in test set): ",nrow(data),"\n",
		"Total number of authors (in test set): ",length(unique(data$author_id)),"\n",
		"Total number of  voters (in test set): ",length(unique(data$voter_id)),"\n", sep="");

	out = data[data$voter_id %in% user.train & data$author_id %in% user.train,];
	
	cat("Selected number of ratings (in test set): ",nrow(out),"\n",
		"Selected number of   users (in test set): ",length(unique(c(out.train$author_id, out.train$voter_id))),"\n",
		"Selected number of authors (in test set): ",length(unique(out$author_id)),"\n",
		"Selected number of  voters (in test set): ",length(unique(out$voter_id)),"\n",sep="");
	
	if(converTimestamp){
		out$comment_time = str2timestamp(out$comment_time);
		out$rating_time  = str2timestamp(out$rating_time);
	}
	
	out.test = out;
	out.all = rbind(out.train, out.test);
	
	data = list(train=out.train, test=out.test);
	if(!is.null(out.file)){
		save(data, file=paste(out.file,".RData",sep=""));
		write.table(data$train, file=paste(out.file,".train",sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE);
		write.table(data$test,  file=paste(out.file,".test",sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE);
	}
	
	return(data);
}

format.timestamp = function(ts, format="%Y-%m-%d %H:%M:%S", tz="EST8EDT"){
	t = ISOdatetime(1970,1,1,0,0,0,tz="GMT") + ts;
	format(t, format, tz=tz)
}

str2timestamp = function(str, tz="EST8EDT"){
	time = strptime(str, "%Y-%m-%d %H:%M:%S", tz=tz);
	begin = ISOdatetime(1970,1,1,0,0,0,tz="GMT");
	return(as.integer(difftime(time, begin, units="secs")));
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
		names(output)[length(by)+(1:length(attr))] = attr;
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
		for(i in 1:length(group.by)){
			if(     is.factor( table[1,group.by[i]])) output[,i] = factor(output[,i], levels=levels(table[1,group.by[i]]))
			else if(is.integer(table[1,group.by[i]])) output[,i] = as.integer(output[,i])
			else if(is.numeric(table[1,group.by[i]])) output[,i] = as.numeric(output[,i]);
		}
	}else{
		names(output) = c(group.by, attr);
	}
	
	return(output);
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
	}else if(is.logical(x)){
		out = x & TRUE;
	}else if(is.null(x)){
		out = NULL;
	}else stop("Type not supported");
	
	return(out);
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


check_type_size <- function(x, type, size, isNullOK=FALSE, check.NA=TRUE, name=NULL){
	if(is.null(name)) name = "The input"
	if(is.null(x)){
		if(isNullOK) return(TRUE);
		if(is.list(size)){
			for(i in 1:length(size)) if(any(size[[i]] == 0)) return(TRUE);
			stop(name," is null");
		}
		if(any(size == 0)) return(TRUE)
		else stop(name," is null");
	}
	if(type == "double"){
		if(!is.double(x)) stop(name," should be double");
	}else if(type == "integer" || type == "int"){
		if(!is.integer(x)) stop(name," should be integer");
	}else stop("Unknown type: ",type,sep="");
	
	d = dim(x);
	if(is.null(d)) d = length(x);
	if(is.list(size)){
		dim.correct = FALSE;
		for(i in 1:length(size)){
			if(length(d) == length(size[[i]]) && all(d == size[[i]])){dim.correct = TRUE; break;}
		}
		if(!dim.correct){
			if(length(size) == 1) stop(name," has dimensionality mismatch: (",paste(d,collapse=" x "),") vs (",paste(size[[1]],collapse=" x "),")");
			if(length(size) == 2) stop(name," has dimensionality mismatch: (",paste(d,collapse=" x "),") vs (",paste(size[[1]],collapse=" x "),") or (",paste(size[[2]],collapse=" x "),")");
			if(length(size) == 3) stop(name," has dimensionality mismatch: (",paste(d,collapse=" x "),") vs (",paste(size[[1]],collapse=" x "),") or (",paste(size[[2]],collapse=" x "),") or (",paste(size[[3]],collapse=" x "),")");
			stop(name," has dimensionality mismatch: (",paste(d,collapse=" x "),") vs (",paste(size[[1]],collapse=" x "),") or (",paste(size[[2]],collapse=" x "),") or (",paste(size[[3]],collapse=" x "),") ...");
		}
	}else{
		if(length(d) != length(size) || any(d != size)) stop(name," has dimensionality mismatch: (",paste(d,collapse=" x "),") vs (",paste(size,collapse=" x "),")");
	}
	
	if(check.NA && any(is.na(x))) stop(name," has some elements that are NA");
}

check_size <- function(x, size, isNullOK=FALSE, check.NA=FALSE){
	if(is.null(x)){
		if(all(size == 0) || isNullOK) return(TRUE)
		else stop("The input is null");
	}
	d = dim(x);
	if(is.null(d)) d = length(x);
	if(length(d) != length(size) || any(d != size)) stop("Dimensionality mismatch: (",paste(d,collapse=" x "),") vs (",paste(size,collapse=" x "),")");
	if(check.NA && any(is.na(x))) stop("Some elements are NA");
}

boxplot.quantile <- function(x, y, prob=(0:10)/10, ...){
	if(length(x) != length(y)) stop("length(x) != length(y)");
	if(any(prob > 1) || any(prob < 0)) stop("any(prob > 1) || any(prob < 0)");
	prob = sort(prob);

	ignore = is.na(x) | is.na(y);
	x = x[!ignore]; y = y[!ignore];
	rnd = order(runif(length(x)));
	x = x[rnd];  y = y[rnd];
	ord = order(x);
	x = x[ord]; y = y[ord];
	label = rep("", length(x));
	start = floor(length(x)*prob[1])+1;
	for(i in 2:length(prob)){
		label[start:floor(length(x)*prob[i])] = sprintf("%02d~%02d%%",round(prob[i-1]*100),round(prob[i]*100));
		start = floor(length(x)*prob[i])+1;
	}
	s = label != "";
	y = y[s];  label = label[s];
	cat("Number of cases: ",length(y),"\n",sep="");
	boxplot(y~label, las=2, ...);
	print(quantile(x, prob=prob));
}

length.unique <- function(x){
	return(length(unique(x)));
}

is.same.data.frame <- function(x, y){
	if(!is.data.frame(x)) stop("x is not a data frame");
	if(!is.data.frame(y)) stop("y is not a data frame");
	if(any(dim(x) != dim(y))) return(FALSE);
	if(!is.null(names(x)) && !is.null(names(y)) && any(names(x) != names(y))) return(FALSE);
	eq = (x == y);
	if(any(is.na(eq))){
		same.na = (is.na(x) == is.na(y));
		eq[is.na(eq)] = same.na[is.na(eq)];
	}
	return(all(eq));
}

###
### Per-group AUC / precision at k / RMSE
###	method: auc, prec, rmse
### pred.y: a vector of column names in data
### true.y and by are column names in data
metric.perGroup <- function(
	data, # data.frame that include column names specified in pred.y, true.y and by
	pred.y, true.y, # column names
	by,   # the name of the group-by column
	method, 
	topK=NULL,
	minObs.perGroup=0
){
	if(!(true.y %in% names(data))) stop("true.y=",true.y," is not in names(data)");
	if(!(by %in% names(data))) stop("by=",by," is not in names(data)");
	if(!all(pred.y %in% names(data))) stop("!all(pred.y %in% names(data))");
	if(!(method %in% c("auc", "prec", "rmse"))) stop("Unknown method: ", method);
	if(method == "prec" && is.null(topK)) stop("topK should not be NULL");
	if(method == "prec" && !all(data[,true.y] %in% c(0,1))) stop("!all(data[,true.y] %in% c(0,1))");
	temp   = aggregate(data[,by], list(group=data[,by]), length);
	groups = temp$group[temp$x >= minObs.perGroup];
	index = tapply(1:nrow(data), list(data[,by]), c, simplify=FALSE);
	out = data.frame(group=groups, stringsAsFactors=FALSE);
	for(name in pred.y) out[,name] = NA;
	i = 1;
	for(k in 1:length(index)){
		group.name = names(index)[k];
		this.indices = index[[k]];
		if(length(this.indices) < minObs.perGroup) next;
		d = data[this.indices, c(pred.y, true.y, by)];
		if(method == "auc" && length(unique(d[,true.y])) == 1) next;
		if(any(d[,by] != group.name)) stop("any(d[,by] != group.name)");
		out[i,"group"] = group.name;
		for(name in pred.y){
			if(method == "auc"){
				pred.roc = prediction(d[,name],d[,true.y]);
				out[i,name] = performance(pred.roc,"auc")@y.values[[1]];
			}else if(method == "prec"){
				x = d[,c(name,true.y)];
				x = x[order(x[,name], decreasing=TRUE),];
				out[i,name] = mean(x[1:topK,true.y]);
			}else if(method == "rmse"){
				out[i,name] = sqrt(mean((d[,name] - d[,true.y])^2));
			}else stop("Unknown method: ", method)
		}
		i = i+1;
	}
	return(out[1:(i-1),]);
}

###
### Draw a random sample from a multivariate normal distribtuion
###     This function is implemented to make sure the C implementation generates exactly
###     the same output as the R implementation
###
my_rmvnorm <- function(n, mu, Sigma=NULL, Sigma.inv=NULL, debug=10, tol=1e-8, verbose=0, always.decompose.Sigma.inv=TRUE){
	nDim <- length(mu);
	if (!is.null(Sigma) && !all(dim(Sigma) == c(nDim, nDim))) stop("Sigma is not a square matrix");
	if (!is.null(Sigma.inv) && !all(dim(Sigma.inv) == c(nDim, nDim))) stop("Sigma.inv is not a square matrix");
	name = NULL;
	if(is.null(Sigma.inv)){
		name = dimnames(Sigma);
		if(always.decompose.Sigma.inv){
			Sigma.inv = solve(Sigma);
			eigen = sym_eigen(Sigma.inv);
			eigenValue = 1/eigen$values;
		}else{
			eigen = sym_eigen(Sigma);
			eigenValue = eigen$values;
		}
	}else if(is.null(Sigma)){
		name = dimnames(Sigma.inv);
		eigen <- sym_eigen(Sigma.inv)
		eigenValue <- 1/eigen$values
	}else{
		error("Please specify one of Sigma or Sigma.inv (not both)");
	}
	
	if(debug >= 3){
		if(is.null(Sigma)) Sigma = solve(Sigma.inv);
		if(max(abs(eigen$vectors %*% diag(eigenValue, nDim) %*% t(eigen$vectors) - Sigma)) > tol * abs(eigenValue[1]))
			stop("sym_eigen(Sigma) seems to have some problems!!");
	}
	
	if (!all(eigenValue >= -tol * abs(eigenValue[1]))) stop("'Sigma' is not positive definite")
	
	rnd = matrix(rnorm(nDim * n), n)
	
	if(verbose >= 10){
		cat("   eigen value = ");   print(drop(eigenValue));
		cat("   eigen vector =\n"); print(eigen$vectors, nDim);
	}
	# cat("temp:\n");
	# print(eigen$vectors %*% diag(sqrt(pmax(eigenValue, 0)), nDim));
	# cat("rnd: ");
	# print(drop(rnd));
	
	out = drop(mu) + eigen$vectors %*% diag(sqrt(pmax(eigenValue, 0)), nDim) %*% t(rnd);
	name2 = NULL;
	if(!is.null(names(mu))){
		name2 = names(mu);
	}else if(!is.null(name)){
		name2 = name[[1]];
	}
	dimnames(out) = list(name2, NULL);
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

cross.table <- function(x, y, 
	grand.total=TRUE,   # whether to show grand total
	fraction.by=NULL,   # "row" or "col" sum up to one
	percentage.by=NULL, # "row" or "col" sum up to one
	...
){
	if(!is.null(fraction.by) && !(fraction.by %in% c("row", "col"))) stop("fraction.by must be either 'row' or 'col'");
	if(!is.null(percentage.by) && !(percentage.by %in% c("row", "col"))) stop("percentage.by must be either 'row' or 'col'");
	if(!is.null(percentage.by) && !is.null(fraction.by)) stop("You cannot specify both fraction.by and percentage.by");
	tab = table(x, y, ...);
	by = NULL; if(!is.null(fraction.by)) by = fraction.by;  if(!is.null(percentage.by)) by = percentage.by;
	multiplier = if(!is.null(percentage.by)) 100 else 1;
	
	if(grand.total || !is.null(by)){
		x.total = aggregate(rep(1,length(x)), list(x), length);
		y.total = aggregate(rep(1,length(y)), list(y), length);
		if(!is.null(by)){
			if(     by == "row"){ temp = matrix(x.total[,2], byrow=FALSE, nrow=nrow(tab),ncol=ncol(tab));}
			else if(by == "col"){ temp = matrix(y.total[,2], byrow=TRUE,  nrow=nrow(tab),ncol=ncol(tab));}
			else stop("what??");
			tab = (tab / temp) * multiplier;
		}
	}
	if(grand.total){
		out = matrix(0, nrow=nrow(tab)+1, ncol=ncol(tab)+1);
		attr(out, "dimnames") = list(c(dimnames(tab)[[1]], "TOTAL"), c(dimnames(tab)[[2]], "TOTAL"));
		out[1:nrow(tab),1:ncol(tab)] = tab;
		out[x.total[,1], ncol(out)] = x.total[,2];
		out[nrow(out), y.total[,1]] = y.total[,2];
		total = sum(x.total[,2]);
		out[nrow(out), ncol(out)] = total;
		if(!is.null(by)){
			if(     by == "row"){ out[nrow(out),1:ncol(tab)] = (out[nrow(out),1:ncol(tab)] / total) * multiplier;}
			else if(by == "col"){ out[1:nrow(tab),ncol(out)] = (out[1:nrow(tab),ncol(out)] / total) * multiplier;}
			else stop("what??");
		}
	}else{
		out = matrix(0, nrow=nrow(tab), ncol=ncol(tab));
		attr(out, "dimnames") = list(c(dimnames(tab)[[1]]), c(dimnames(tab)[[2]]));
		out[1:nrow(tab),1:ncol(tab)] = tab;
	}
	return(out);
}

my.print.table <- function(
	x,    # matrix or data.frame
	num.format=NULL, # format for numbers using sprintf
	latex=FALSE,
	show.rownames=TRUE, show.colnames=TRUE,
	justify="right", begin="", sep=" ", end="",
	...
){
	if(!latex && is.null(num.format)){ print(x,...); return();}
	if(!is.data.frame(x) && !is.matrix(x)) stop("input data must be a matrix or data frame");
	row.name = if(show.rownames) rownames(x) else NULL;
	col.name = if(show.colnames) colnames(x) else NULL;
	if(is.matrix(x)) x = data.frame(x, stringsAsFactors=FALSE);
	M = if(is.null(col.name)) nrow(x) else nrow(x)+1;
	N = if(is.null(row.name)) ncol(x) else ncol(x)+1;
	y = matrix("", nrow=M, ncol=N);
	k=0;
	if(!is.null(row.name)){
		if(!is.null(col.name)) row.name = c("", row.name);
		row.name = format(row.name, justify=justify);
		y[,1] = row.name; k=1;
	}
	for(i in 1:ncol(x)){
		if(is.logical(x[,i])){ x[,i] = x[,i]+0; }
		else if(is.integer(x[,i])){ x[,i] = sprintf("%d", x[,i]); }
		else if(is.numeric(x[,i]) && !is.null(num.format)){ x[,i] = sprintf(num.format, x[,i]); }
		else if(is.factor(x[,i])){ x[,i] = as.character(x[,i]); }
		y[,i+k] = format(c(col.name[i], x[,i]), justify=justify);
	}
	if(latex){  sep=" & "; end=" \\\\";}
	for(i in 1:nrow(y)){
		cat(begin,sep=""); cat(y[i,],sep=sep); cat(end,"\n",sep="");
	}
}

logit <- function(x){
	return(log(x) - log(1-x));
}

###
### Allow default values for NA
###
my.merge <- function(x, y, ..., na.as=NULL){
	out = merge(x=x, y=y, ...);
	if(!is.null(na.as)){
		for(i in 1:ncol(out)){
			select = is.na(out[,i]);
			if(any(select)) out[select,i] = na.as;
		}
	}
	return(out);
}

###
### aggregate.char: Make aggregate faster for characters
### (it doesn't really matter)
###
aggregate.char <- function(x, by, FUN, ..., sep="\t"){
	if(!is.data.frame(by)){
		if(!is.list(by) && is.vector(by)) by = list(Group.1=by);
		if(is.list(by)){
			if(is.null(names(by))) names(by) = paste("Group",1:length(by),sep=".");
			by = data.frame(by, stringsAsFactors=FALSE);
		}else stop("by should be a data.frame, list or vector");
	}
	group.names = names(by);
	str.by = by[,1];
	if(ncol(by) > 1){ for(i in 2:ncol(by)){
		str.by = paste(str.by, by[,i], sep=sep);
	}}
	str.keys = unique(str.by);
	by$..GROUP_ID.. = match(str.by, str.keys);
	temp = aggregate(x=x, by=by["..GROUP_ID.."], FUN=FUN, ...);
	out = merge(by, temp, by="..GROUP_ID..")[,c(group.names, "x")];
	return(out);
}

###
### Smooth ratio  (x + mean*size) / (y + size)
### if(mean is null) mean = sum(x) / sum(y)
###
smooth.ratio <- function(
	x, y, mean=NULL, size
){
	if(is.null(mean)) mean = sum(as.double(x)) / sum(as.double(y));
	return( (x + mean*size) / (y + size) );
}
