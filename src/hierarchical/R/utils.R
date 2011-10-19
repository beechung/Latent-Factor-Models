### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

logit <- function(p){
	return(log(p/(1-p)));
}

inv.logit <- function(x){
	return(1/(1+exp(-x)));
}

plot.required.support <- function(
    data, file, conf=0.95, err=c(0.2, 0.1, 0.05, 0.01), type="fixed", relative.err=TRUE
){
    data = data[data$type == type,];
    postscript(file);
    setting = paste(data$mu, data$cv);
    set = split(data, setting);
    for(i in 1:length(set)){
        data_i = set[[i]];
        plot.required.support.single(data_i, conf=conf, param="mean", err=err, relative.err=relative.err);
        plot.required.support.single(data_i, conf=conf, param="CV", err=err, relative.err=relative.err);
    }
    dev.off();
}

###
### data$CB95.mu=x means Pr(|mu.fitted - mu.truth|/mu.truth < x) > 95%
### data$CB95.cv=x means Pr(|cv.fitted - cv.truth| < x) > 95%
###
plot.required.support.single <- function(
    data, # (nObs, nTrials, type, mu, cv, method, CB90.mu, CB95.mu, CB99.mu, CB90.cv, CB95.cv, CB99.cv)
    conf=0.95, param="mean", err=c(0.2, 0.1, 0.05, 0.01), relative.err=TRUE
){
    plot_data = get.required.support.single(
		data=data, conf=conf, param=param, err=err, relative.err=relative.err
	);

	if(is.null(plot_data)){
        plot_data = data.frame(bound="", nTrials=1e8, nObs=1e8);
    }
    plot.groupBy(x=plot_data$nTrials, y=plot_data$nObs, group.by=plot_data$bound, order.by=plot_data$nTrials, 
        legend.pos="bottomleft", log="xy", xlim=c(10,10000), ylim=c(10,10000),
        main=c(paste("Required sample size for ",conf*100,"% bound of ",param, sep=""), 
               paste("Ground truth: mean=",mu,", CV=",cv,sep="")),
        xlab="Number of views", ylab="Number of items");
}

###
### data$CB95.mu=x means Pr(|mu.fitted - mu.truth|/mu.truth < x) > 95%
### data$CB95.cv=x means Pr(|cv.fitted - cv.truth| < x) > 95%
###
get.required.support.single <- function(
	data, # (nObs, nTrials, type, mu, cv, method, CB90.mu, CB95.mu, CB99.mu, CB90.cv, CB95.cv, CB99.cv)
	conf=0.95, param="mean", err=c(0.2, 0.1, 0.05, 0.01), relative.err=TRUE
){
	if(nrow(data) == 0) stop("nrow(data) == 0");
	mu = data$mu[1];  cv = data$cv[1];  method = data$method[1];
	if(any(data$mu != mu) || any(data$cv != cv) || any(data$type != data$type[1]) || any(data$method != method)) stop("Multiple setting fixed together");
	out = NULL;
	if(param == "mean"){
		attr = paste("CB",conf*100,".mu",sep="");
		if(relative.err) bound = data[,attr]
		else             bound = data[,attr] * mu;
	}else if(param == "CV"){
		attr = paste("CB",conf*100,".cv",sep="");
		if(relative.err) bound = data[,attr] / cv
		else             bound = data[,attr];
	}
	for(e in err){
		temp  = data[bound <= e, c("nObs", "nTrials")];
		if(nrow(temp) == 0) next;
		temp2 = aggregate(temp$nObs, list(temp$nTrials), min);
		if(relative.err) err.type = "relative"
		else             err.type = "absolute";
		temp3 = data.frame(
					method=method, mu=mu, cv=cv, error=e, param=param, err.type=err.type, 
					nTrials=temp2[,1], nObs=temp2[,2], sim.type=data$type[1]
			    );
		temp3 = temp3[order(temp3$nTrials),];
		if(nrow(temp3) > 1){
			for(i in nrow(temp3):2){
				if(temp3$nObs[i] > temp3$nObs[i-1]) temp3$nObs[i-1] = temp3$nObs[i];
			}
		}
		out = rbind(out, temp3);
	}
	return(out);
}

get.required.support <- function(
	data, conf=0.95, err=c(0.2, 0.1, 0.05, 0.01), 
	type="fixed", param=c("mean", "CV"), relative.err=TRUE
){
	data = data[data$type == type,];
	setting = paste(data$mu, data$cv);
	set = split(data, setting);
	out = NULL;
	for(i in 1:length(set)){
		data_i = set[[i]];
		for(j in 1:length(param)){
			temp = get.required.support.single(data_i, conf=conf, param=param[j], err=err, relative.err=relative.err);
			out = rbind(out, temp);
		}
	}
	return(out);
}

###
###  INPUT: dag = data.frame(nodeID, parentID)
###
get.all.leaves <- function(dag){
	temp = unique(dag$nodeID);
	temp = temp[!(temp %in% unique(dag$parentID))];
	return(temp);
}
