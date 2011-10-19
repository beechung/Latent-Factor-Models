### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
### 
### Author: Bee-Chung Chen

check_type_size <- function(x, type, size, isNullOK=FALSE, check.NA=TRUE){
	if(is.null(x)){
		if(all(size == 0) || isNullOK) return(TRUE)
		else stop("The input is null");
	}
	if(type == "double"){
		if(!is.double(x)) stop("The input should be double");
	}else if(type == "integer" || type == "int"){
		if(!is.integer(x)) stop("The input should be integer");
	}else stop("Unknown type: ",type,sep="");
	d = dim(x);
	if(is.null(d)) d = length(x);
	if(length(d) != length(size) || any(d != size)) stop("Dimensionality mismatch: (",paste(d,collapse=" x "),") vs (",paste(size,collapse=" x "),")");
	if(check.NA && any(is.na(x))) stop("Some elements are NA");
}

### Fit a Beta-Binomial DAG model
###
### Input: obs.data = data.frame(pNodeID, tID, nPositives, nTrials)
###					  pNodeID: ID of a population node
###					  tID: ID of a treatment
###        pop.dag  = data.frame(nodeID, parentID)
###		   tre.dag  = data.frame(nodeID, parentID)
###		   tre.map  = data.frame(tID, tNodeID)
###
### NODE ID ORDERING:
###   * Root should have ID = 0.
###   * Population nodes and treatment nodes are labeled separately; i.e., both start from 0.
###   * For any given node, every one of its ancestors should must ID < the ID of the given node.
###
### INPUT OBSERVATION DATA ORDERING:
###   * Observation data table: data_(pNodeID, tID, nPositives, nTrials)
###   * This table must be ordered by pNodeID
###
### Output:
###	  param[pnode+1,tnode+1,] = (mean, size) of the Beta distribution for (pnode, tnode)
###   status[pnode+1,tnode+1] = 1: Moment estimate
###								2: Posterior mode
###								0: Low support (no estimate)
###							   -1: Invalid (the node was not reached)
###							   -2: Should compute posterior mode, but only the moment estimate is ready
###							   -3: Posterior mode computation failed
### Options:
###		option=0: moment estimates only
###		option=1: postierior mode (use the likelihood criterion for selecting immediate parents)
###		option=2: postierior mode (use fixed priors for all nodes, specified by root_priorMean and root_priorSize)
###
fit.BetaBinomialDAG.C <- function(
	obs.data, # observation data: (pNodeID, tID, nPositives, nTrials)
	pop.dag,  # population DAG: (nodeID, parentID)
	tre.dag,  # treatment  DAG: (nodeID, parentID)
	tre.map,  # treatment ID to node mapping: (tID, tNodeID)
	option=1, # -1: naive estimates; 0: moment estimates only;  1: also estimate postierior modes
	ME_nRefines=0,  # number of refinements for the moment estimator
	root_priorMean, # prior for the mean of Beta at root
	root_priorSize, # prior for the size of Beta at root
	var_logitMean,  # variance of logit(mu) from its parent node
	var_logSize,    # variance of log(size) from its parent node
	ME_threshold_numObs,    # Use moment estimates for notes with
	ME_threshold_numTrials, # #obs >= ME_threshold_numObs and avg(#trial) >= ME_threshold_numTrials
	nTrials_ForSmall=0, nObs_ForLarge=1e6, epsilon=1e-10, stepSize=1, 
	minMean=1e-10, maxMean=1-1e-10, minSize=1e-10, maxSize=1e10,
	minNumObs=3, minTotalTrials=10, minTotalPostives=0,
	maxIter1=20, nLnsrchStep1=0, maxIter2=100, nLnsrchStep2=10,
	verbose=0, debug=0
){
	data.length = as.integer(nrow(obs.data));
	check_type_size(obs.data$pNodeID, "int", data.length);
	check_type_size(obs.data$tID, "int", data.length);
	check_type_size(obs.data$nPositives, "int", data.length);
	check_type_size(obs.data$nTrials, "int", data.length);
	
	nPNodes = as.integer(max(pop.dag$nodeID)+1);
	nTNodes = as.integer(max(tre.dag$nodeID)+1);
	nTIDs   = as.integer(max(tre.map$tID)+1);
	
	num_pDAG_edges = as.integer(nrow(pop.dag));
	check_type_size(pop.dag$nodeID, "int", num_pDAG_edges);
	check_type_size(pop.dag$parentID, "int", num_pDAG_edges);
	
	num_tDAG_edges = as.integer(nrow(tre.dag));
	check_type_size(tre.dag$nodeID, "int", num_tDAG_edges);
	check_type_size(tre.dag$parentID, "int", num_tDAG_edges);
	
	num_tMap_edges = as.integer(nrow(tre.map));
	check_type_size(tre.map$tID, "int", num_tMap_edges);
	check_type_size(tre.map$tNodeID, "int", num_tMap_edges);
	
	out = list(
		param  = array(0.0, dim=c(nPNodes, nTNodes, 2)),
		status = array(as.integer(0), dim=c(nPNodes, nTNodes))
	);
	
	.C("fit_BetaBinomialDAG",
		# OUTPUT
		out$param, out$status,
		# INPUT observation data
		obs.data$pNodeID, obs.data$tID, obs.data$nPositives, obs.data$nTrials, data.length,
		# INPUT population & treatment DAG
		pop.dag$nodeID, pop.dag$parentID, nPNodes, num_pDAG_edges,
		tre.dag$nodeID, tre.dag$parentID, nTNodes, num_tDAG_edges,
		# INPUT mapping from treatment IDs to nodes
		tre.map$tID, tre.map$tNodeID, nTIDs, num_tMap_edges,
		# INPUT for the estimation procedures
		as.integer(ME_nRefines),  as.integer(option), 
		as.double(root_priorMean), as.double(root_priorSize), as.double(var_logitMean), as.double(var_logSize),
		as.integer(ME_threshold_numObs), as.double(ME_threshold_numTrials),
		as.integer(nTrials_ForSmall), as.integer(nObs_ForLarge), as.double(epsilon),
		as.double(stepSize), as.integer(maxIter1), as.integer(nLnsrchStep1), as.integer(maxIter2), as.integer(nLnsrchStep2),
		as.double(minMean), as.double(maxMean), as.double(minSize), as.double(maxSize),
		as.integer(minNumObs), as.integer(minTotalTrials), as.integer(minTotalPostives),
		# INPUT other options
		as.integer(verbose), as.integer(debug),
		DUP=FALSE
	);
	
	pnodes = paste("P",0:(nPNodes-1),sep="");
	tnodes = paste("T",0:(nTNodes-1),sep="")
	dimnames(out$param) = list(pnodes, tnodes, c("mean", "size"));
	dimnames(out$status) = list(pnodes, tnodes);
	
	out$numLowSupportNodes = sum(out$status == 0);
	out$numFailedNodes     = sum(out$status < 0);
	
	return(out);
}

###
### Compute the log likelihood of the observations
###		option=0: Compute log likelihood for the leaves of the treatment DAG
###		option=1: Compute log likelihood for all the nodes
###		option=2: Compute log likelihood for all the frontiers (of the treatment DAG),
###				  which are leaves or have at least one child with low support (status=0)
###     avg.loglik=TRUE:  Compute E[log-likelihood] over observations for each node
###                       Then, compute the E[log-likelihood] over nodes
###     avg.loglik=FALSE: sum(log-likelihood) over all observations
### * To ignore high-level nodes in the population DAG, please remove them from obs.data
###   before calling this function
loglik.BetaBinomialDAG <- function(
	model,    # list(param, status)
	obs.data, # observation data: (pNodeID, tID, nPositives, nTrials)
	pop.dag,  # population DAG: (nodeID, parentID)
	tre.dag,  # treatment  DAG: (nodeID, parentID)
	tre.map,  # treatment ID to node mapping: (tID, tNodeID)
	option=1, # 0: leaves only, 1: the whole DAG, 2: frontiers only
	output.nodeLoglik=TRUE,
	nTrials_ForSmall=0, nObs_ForLarge=1e6, 
	minMean=1e-10, maxMean=1-1e-10, minSize=1e-10, maxSize=1e10,
	avg.loglik=TRUE, # whether to compute the average log likelihood per observation
	verbose=0, debug=0
){
	data.length = as.integer(nrow(obs.data));
	check_type_size(obs.data$pNodeID, "int", data.length);
	check_type_size(obs.data$tID, "int", data.length);
	check_type_size(obs.data$nPositives, "int", data.length);
	check_type_size(obs.data$nTrials, "int", data.length);
	
	nPNodes = as.integer(max(pop.dag$nodeID)+1);
	nTNodes = as.integer(max(tre.dag$nodeID)+1);
	nTIDs   = as.integer(max(tre.map$tID)+1);
	
	num_pDAG_edges = as.integer(nrow(pop.dag));
	check_type_size(pop.dag$nodeID, "int", num_pDAG_edges);
	check_type_size(pop.dag$parentID, "int", num_pDAG_edges);
	
	num_tDAG_edges = as.integer(nrow(tre.dag));
	check_type_size(tre.dag$nodeID, "int", num_tDAG_edges);
	check_type_size(tre.dag$parentID, "int", num_tDAG_edges);
	
	num_tMap_edges = as.integer(nrow(tre.map));
	check_type_size(tre.map$tID, "int", num_tMap_edges);
	check_type_size(tre.map$tNodeID, "int", num_tMap_edges);
	
	check_type_size(model$param, "double", c(nPNodes, nTNodes, 2));
	check_type_size(model$status, "int", c(nPNodes, nTNodes));

	nodeLoglik = NULL;  output_nodeLoglik = 0;
	if(output.nodeLoglik){
		nodeLoglik = array(0.0, dim=c(nPNodes, nTNodes));
		output_nodeLoglik = 1;
	}
	avgloglik = as.integer(0);
	if(avg.loglik){
		avgloglik = as.integer(1);
	}
	loglik = c(0.0);
		
	.C("loglik_BetaBinomialDAG",
			# OUTPUT
			loglik, nodeLoglik,
			# INPUT model
			model$param, model$status,
			# INPUT observation data
			obs.data$pNodeID, obs.data$tID, obs.data$nPositives, obs.data$nTrials, data.length,
			# INPUT population & treatment DAG
			pop.dag$nodeID, pop.dag$parentID, nPNodes, num_pDAG_edges,
			tre.dag$nodeID, tre.dag$parentID, nTNodes, num_tDAG_edges,
			# INPUT mapping from treatment IDs to nodes
			tre.map$tID, tre.map$tNodeID, nTIDs, num_tMap_edges,
			# OPTIONS
			as.integer(option), as.integer(output_nodeLoglik),
			as.integer(nTrials_ForSmall), as.integer(nObs_ForLarge),
			as.double(minMean), as.double(maxMean), as.double(minSize), as.double(maxSize),
			avgloglik,
			# INPUT other options
			as.integer(verbose), as.integer(debug),
			DUP=FALSE
	);
	
	out = list(loglik=loglik);
	
	if(output.nodeLoglik){
		out$node.loglik = nodeLoglik;
		if(option == 1){
			
			if(avg.loglik) out$avg.entire.loglik = loglik
			else           out$sum.entire.loglik = loglik;
			
			# leaf loglik
			selected = array(FALSE, dim=c(nPNodes, nTNodes));
			pLeafIDs = get.all.leaves(pop.dag);
			tLeafIDs = get.all.leaves(tre.dag);
			selected[pLeafIDs+1,tLeafIDs+1] = TRUE;
			selected[nodeLoglik == 0] = FALSE;
			
			if(avg.loglik) out$avg.leaf.loglik = mean(nodeLoglik[selected])
			else           out$sum.leaf.loglik =  sum(nodeLoglik[selected]);
			
			out$nLeaves.valid = sum(selected);
			out$nLeaves = length(pLeafIDs) * length(tLeafIDs);
			# frontier loglik
			selected = array(FALSE, dim=c(nPNodes, nTNodes));
			indices = get.frontierIndices(status=model$status, pop.dag=pop.dag, tre.dag=tre.dag);
			selected[indices] = TRUE;
			selected[nodeLoglik == 0] = FALSE;

			if(avg.loglik) out$avg.frontier.loglik = mean(nodeLoglik[selected])
			else           out$sum.frontier.loglik =  sum(nodeLoglik[selected]);
			
			out$nFrontiers.valid = sum(selected);
			out$nFrontiers = nrow(indices);
		}
	}
	
	return(out);
}

CI.Coverage.BetaBinomialDAG <- function(
		model,    # list(param, status)
		prob,     # (prob*100)% confidence interval will be constructed
		obs.data, # observation data: (pNodeID, tID, nPositives, nTrials)
		pop.dag,  # population DAG: (nodeID, parentID)
		tre.dag,  # treatment  DAG: (nodeID, parentID)
		tre.map,  # treatment ID to node mapping: (tID, tNodeID)
		option=1, # 0: leaves only, 1: the whole DAG, 2: frontiers only
		output.nodeStats=TRUE,
		minMean=1e-10, maxMean=1-1e-10, minSize=1e-10, maxSize=1e10,
		verbose=0, debug=0
){
	check_type_size(prob, "double", 1);
	
	data.length = as.integer(nrow(obs.data));
	check_type_size(obs.data$pNodeID, "int", data.length);
	check_type_size(obs.data$tID, "int", data.length);
	check_type_size(obs.data$nPositives, "int", data.length);
	check_type_size(obs.data$nTrials, "int", data.length);
	
	nPNodes = as.integer(max(pop.dag$nodeID)+1);
	nTNodes = as.integer(max(tre.dag$nodeID)+1);
	nTIDs   = as.integer(max(tre.map$tID)+1);
	
	num_pDAG_edges = as.integer(nrow(pop.dag));
	check_type_size(pop.dag$nodeID, "int", num_pDAG_edges);
	check_type_size(pop.dag$parentID, "int", num_pDAG_edges);
	
	num_tDAG_edges = as.integer(nrow(tre.dag));
	check_type_size(tre.dag$nodeID, "int", num_tDAG_edges);
	check_type_size(tre.dag$parentID, "int", num_tDAG_edges);
	
	num_tMap_edges = as.integer(nrow(tre.map));
	check_type_size(tre.map$tID, "int", num_tMap_edges);
	check_type_size(tre.map$tNodeID, "int", num_tMap_edges);
	
	check_type_size(model$param, "double", c(nPNodes, nTNodes, 2));
	check_type_size(model$status, "int", c(nPNodes, nTNodes));
	
	node.nObs = NULL;  node.nObsInCI = NULL;  output_nodeStats = 0;
	if(output.nodeStats){
		node.nObs = array(as.integer(0), dim=c(nPNodes, nTNodes));
		node.nObsInCI = array(as.double(0), dim=c(nPNodes, nTNodes));
		output_nodeStats = 1;
	}
	
	select = NULL;  useSelect = 0;
	if(option %in% c(0, 2)){
		useSelect = 1;
		select = array(as.integer(0), dim=c(nPNodes, nTNodes));
		if(option == 0){
			pLeafIDs = get.all.leaves(pop.dag);
			tLeafIDs = get.all.leaves(tre.dag);
			select[pLeafIDs+1,tLeafIDs+1] = as.integer(1);
		}else if(option == 2){
			frontier.indices = get.frontierIndices(status=model$status, pop.dag=pop.dag, tre.dag=tre.dag);
			select[frontier.indices] = as.integer(1);
		}else stop("error");
	}

	coverage = c(0.0);  nFailed = as.integer(0);
	
	.C("confIntCoverage_BetaBinomialDAG",
			# OUTPUT
			coverage, node.nObsInCI, node.nObs, nFailed,
			# INPUT model
			model$param, model$status, prob, select,
			# INPUT observation data
			obs.data$pNodeID, obs.data$tID, obs.data$nPositives, obs.data$nTrials, data.length,
			# INPUT population & treatment DAG
			pop.dag$nodeID, pop.dag$parentID, nPNodes, num_pDAG_edges,
			tre.dag$nodeID, tre.dag$parentID, nTNodes, num_tDAG_edges,
			# INPUT mapping from treatment IDs to nodes
			tre.map$tID, tre.map$tNodeID, nTIDs, num_tMap_edges,
			# OPTIONS
			as.integer(useSelect), as.integer(output_nodeStats),
			as.double(minMean), as.double(maxMean), as.double(minSize), as.double(maxSize),
			# INPUT other options
			as.integer(verbose), as.integer(debug),
			DUP=FALSE
	);
	
	out = list(coverage=coverage+0, nFailed=nFailed+0);
	
	if(output.nodeStats){
		
		out$node.nObs = node.nObs;
		out$node.nObsInCI = node.nObsInCI;

		if(option %in% c(1,2)){
			
			# leaf coverage
			selected = array(FALSE, dim=c(nPNodes, nTNodes));
			pLeafIDs = get.all.leaves(pop.dag);
			tLeafIDs = get.all.leaves(tre.dag);
			selected[pLeafIDs+1,tLeafIDs+1] = TRUE;
			selected[node.nObs == 0] = FALSE;
			
			out$leaf.coverage = sum(node.nObsInCI[selected]) / sum(node.nObs[selected]);
			
			out$nLeaves.valid = sum(selected);
			out$nLeaves = length(pLeafIDs) * length(tLeafIDs);
		}
		if(option == 1){
			# frontier coverage
			selected = array(FALSE, dim=c(nPNodes, nTNodes));
			indices = get.frontierIndices(status=model$status, pop.dag=pop.dag, tre.dag=tre.dag);
			selected[indices] = TRUE;
			selected[node.nObs == 0] = FALSE;
			
			out$frontier.coverage = sum(node.nObsInCI[selected]) / sum(node.nObs[selected]);
			
			out$nFrontiers.valid = sum(selected);
			out$nFrontiers = nrow(indices);
		}else{
			out$frontier.coverage = coverage+0;
		}
	}
	
	return(out);
}

cubeAggregate.pDAG <- function(
	obs.data, # observation data: (pNodeID, tID, nPositives, nTrials)
	pop.dag,  # population DAG: (nodeID, parentID)
	verbose=0, debug=0
){
	# check
	pLeaves = get.all.leaves(pop.dag);
	if(!all(unique(obs.data$pNodeID) %in% pLeaves)) stop("Some observations are not at the leaf level");
	
	obs.data = obs.data[order(obs.data$tID),];
	
	data.length = as.integer(nrow(obs.data));
	check_type_size(obs.data$pNodeID, "int", data.length);
	check_type_size(obs.data$tID, "int", data.length);
	check_type_size(obs.data$nPositives, "int", data.length);
	check_type_size(obs.data$nTrials, "int", data.length);
	
	nPNodes = as.integer(max(pop.dag$nodeID)+1);
	
	num_pDAG_edges = as.integer(nrow(pop.dag));
	check_type_size(pop.dag$nodeID, "int", num_pDAG_edges);
	check_type_size(pop.dag$parentID, "int", num_pDAG_edges);
	
	nObs = c(as.integer(0));
	
	.C("cubeAggregate_pDAG_get_nObs",
		# OUTPUT
		nObs,
		# INPUT observation data
		obs.data$pNodeID, obs.data$tID, data.length,
		# INPUT population DAG
		pop.dag$nodeID, pop.dag$parentID, nPNodes, num_pDAG_edges,
		# INPUT other options
		as.integer(verbose), as.integer(debug),
		DUP=FALSE
	);
	
	out = data.frame(pNodeID=integer(nObs), tID=integer(nObs), nPositives=integer(nObs), nTrials=integer(nObs));
	
	.C("cubeAggregate_pDAG",
		# OUTPUT
		out$pNodeID, out$tID, out$nPositives, out$nTrials,
		# INPUT
		as.integer(nrow(out)),
		# INPUT observation data
		obs.data$pNodeID, obs.data$tID, obs.data$nPositives, obs.data$nTrials, data.length,
		# INPUT population DAG
		pop.dag$nodeID, pop.dag$parentID, nPNodes, num_pDAG_edges,
		# INPUT other options
		as.integer(verbose), as.integer(debug),
		DUP=FALSE
	);
	
	out = out[order(out$pNodeID, out$tID),];
	
	return(out);
}

get.frontierIndices <- function(
	status,   # status[pnode,tnode] = 0 means low support
	pop.dag,  # population DAG: (nodeID, parentID)
	tre.dag,  # treatment  DAG: (nodeID, parentID)
	verbose=0, debug=0
){
	nPNodes = as.integer(nrow(status));
	nTNodes = as.integer(ncol(status));
	check_type_size(status, "int", c(nPNodes, nTNodes));
	
	num_pDAG_edges = as.integer(nrow(pop.dag));
	check_type_size(pop.dag$nodeID, "int", num_pDAG_edges);
	check_type_size(pop.dag$parentID, "int", num_pDAG_edges);
	
	num_tDAG_edges = as.integer(nrow(tre.dag));
	check_type_size(tre.dag$nodeID, "int", num_tDAG_edges);
	check_type_size(tre.dag$parentID, "int", num_tDAG_edges);
	
	num = c(as.integer(0));
	.C("get_frontierIndices",
		NULL, num, status,
		# INPUT population & treatment DAG
		pop.dag$nodeID, pop.dag$parentID, nPNodes, num_pDAG_edges,
		tre.dag$nodeID, tre.dag$parentID, nTNodes, num_tDAG_edges,
		# INPUT other options
		as.integer(verbose), as.integer(debug),
		DUP=FALSE
	);
	indices = array(integer(1), dim=c(num,2));
	.C("get_frontierIndices",
			indices, num, status,
			# INPUT population & treatment DAG
			pop.dag$nodeID, pop.dag$parentID, nPNodes, num_pDAG_edges,
			tre.dag$nodeID, tre.dag$parentID, nTNodes, num_tDAG_edges,
			# INPUT other options
			as.integer(verbose), as.integer(debug),
			DUP=FALSE
	);
	return(indices+as.integer(1));
}

fit.BetaBin <- function(
	pos, size, mu0, gamma0, var_mu, var_gamma,
	threshold=0, nSamplesForLarge=1e6, epsilon=1e-10, stepSize=1, 
	maxIter1=20, nLnsrchStep1=0, maxIter2=100, nLnsrchStep2=10, 
	minMean=1e-10, maxMean=1-1e-10, minSize=1e-10, maxSize=1e10,
	verbose=0, debug=0
){
	mu = double(1);    gamma = double(1);  code = integer(1);  nIter = integer(1);
	mu.ME = double(1); gamma.ME = double(1);
	len = length(pos);
	check_type_size(pos,"int", len);
	check_type_size(size,"int", len);
	check_type_size(len,"int", 1);
		
	ans = .C("fit_BetaBinomial_mean_size",
			# OUTPUT
			mu, gamma, code, nIter, mu.ME, gamma.ME,
			# INPUT
			pos, size, len, as.integer(threshold), as.integer(nSamplesForLarge),
			as.double(mu0), as.double(var_mu), as.double(gamma0), as.double(var_gamma),
			as.double(epsilon), as.double(stepSize), 
			as.integer(maxIter1), as.integer(nLnsrchStep1), as.integer(maxIter2), as.integer(nLnsrchStep2),
			as.double(minMean), as.double(maxMean), as.double(minSize), as.double(maxSize),
			as.integer(verbose), as.integer(debug),
			DUP=FALSE
	);
	
	out = c(mu=mu, gamma=gamma);
	attr(out, "code") = code;
	attr(out, "nIter") = nIter;
	attr(out, "MomEst") = c(mu=mu.ME, gamma=gamma.ME);
	return(out);
}

gen.cube.DAG <- function(tree.list){
	cube = NULL;
	for(d in 1:length(tree.list)){
		temp = tree.list[[d]];
		temp$dim = as.integer(d-1);
		cube = rbind(cube, temp);
	}
	
	check_type_size(cube$dim, "int", nrow(cube));
	check_type_size(cube$nodeID, "int", nrow(cube));
	check_type_size(cube$parentID, "int", nrow(cube));
	nDagNodes = integer(1);  nDagEdges = integer(1);
	
	ans = .C("get_numDagEdges_from_cube",
		nDagNodes, nDagEdges,
		cube$dim, cube$nodeID, cube$parentID, 
		as.integer(nrow(cube)), as.integer(length(tree.list)),
		DUP=FALSE
	);
	
	dag = data.frame(nodeID=integer(nDagEdges), parentID=integer(nDagEdges));
	map = list(nodeID=integer(nDagNodes), dimVector=matrix(as.integer(0), nrow=nDagNodes, ncol=length(tree.list)));
	
	ans = .C("gen_DAG_from_cube",
		map$nodeID, map$dimVector, dag$nodeID, dag$parentID,
		cube$dim, cube$nodeID, cube$parentID, as.integer(nrow(cube)), 
		as.integer(length(tree.list)), nDagNodes, nDagEdges,
		DUP=FALSE
	);
	dimnames(map$dimVector) = list(map$nodeID, NULL);
	out = list(dag=dag, map=map);
	return(out);
}


fit.BetaBin.old <- function(
	pos, size, mu0, gamma0, var_mu, var_gamma,
	threshold=0, nSamplesForLarge=1e6,
	epsilon=1e-10, stepSize=1, maxIter=20, nLnsrchSteps=10, transform=TRUE, 
	verbose=0, debug=0
){
	mu = double(1);  gamma = double(1);  code = integer(1);  nIter = integer(1);
	len = length(pos);
	check_type_size(pos,"int", len);
	check_type_size(size,"int", len);
	check_type_size(len,"int", 1);

	if(transform) trans = as.integer(1)
	else          trans = as.integer(0);
	
	ans = .C("fit_BetaBinomial_mean_size",
		# OUTPUT
		mu, gamma, code, nIter,
		# INPUT
		pos, size, len, as.integer(threshold), as.integer(nSamplesForLarge),
		as.double(mu0), as.double(var_mu), as.double(gamma0), as.double(var_gamma),
		as.double(epsilon), as.double(stepSize), as.integer(maxIter), as.integer(nLnsrchSteps), trans,
		as.integer(verbose), as.integer(debug),
		DUP=FALSE
	);
	
	out = c(mu=mu, gamma=gamma);
	attr(out, "code") = code;
	attr(out, "nIter") = nIter;
	return(out);
}
