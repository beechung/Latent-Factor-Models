### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
###
### Author: Liang Zhang
fit.bst <- function(
	code.dir = "", # The top-level directory of where code get installed, "" if you are in that directory
	obs.train, # The training response data
	obs.test = NULL, # The testing response data
	x_obs.train = NULL, # The data of training observation features
	x_obs.test = NULL, # The data of testing observation features
	x_src = NULL, # The data of context features for source nodes
	x_dst = NULL, # The data of context features for destination nodes
	x_ctx = NULL, # The data of context features for edges
	out.dir = "", # The directory of output files
	model.name = "model", #The name of the model, can be any string or a vector of strings
	nFactors, # Number of factors, can be any positive integer or vector of positive intergers with length=length(model.name)
	nIter = 20, # Number of EM iterations
	nSamplesPerIter = 200, # Number of Gibbs samples per E step, can be a vector of numbers with length=nIter
	is.logistic = FALSE, # whether to use the logistic model for binary rating, default Gaussian, can be a vector of booleans with length=length(model.name)
	src.dst.same = FALSE, # Whether src_id and dst_id are the same
	control = fit.bst.control(), # A list of control params
	...
) {
  library(Matrix);
  if (code.dir!="") code.dir = sprintf("%s/",code.dir);
  #if (out.dir!="") out.dir = sprintf("%s/",out.dir);
  
  if (floor(nIter)!=nIter || nIter<=0 || length(nIter)>1) stop("When calling fit.bst, nIter must be a positive integer scalar!");
  if (floor(nSamplesPerIter)!=nSamplesPerIter || nSamplesPerIter<=0 || length(nSamplesPerIter)>1) stop("When calling fit.bst, nSamplesPerIter must be a positive integer scalar!");

  # Load all the required libraries and source code  
  if (class(try(load.code(code.dir)))=="try-error") stop("When calling fit.bst, code.dir is not specified correctly. code.dir must be the path of the top-level directory of this package (which contains src/ as a sub-directory and file LICENSE)."); 
  
  # Make sure all the data have required columns
  if (is.null(obs.train$src_id) || is.null(obs.train$dst_id) || is.null(obs.train$y)) stop("When calling fit.bst, obs.train must have the following three columns: src_id, dst_id, y");
  if (!is.null(obs.test)) {
    if (is.null(obs.test$src_id) || is.null(obs.test$dst_id) || is.null(obs.test$y)) stop("When calling fit.bst, obs.test must have the following three columns: src_id, dst_id, y");
  }

  if (is.null(x_obs.train) && !is.null(x_obs.test)) stop("When calling fit.bst, x_obs.train does not exist while x_obs.test is used!  Either set both to NULL or supply both observation feature tables.");
  if (is.null(x_obs.test) && !is.null(x_obs.train) && !is.null(obs.test)) stop("When calling fit.bst, x_obs.test does not exist while x_obs.train is used!  Either set both to NULL or supply both observation feature tables.");
  if (!is.null(x_obs.train) && !is.null(x_obs.test)) {
    if (ncol(x_obs.train)!=ncol(x_obs.test)) stop("When calling fit.bst, x_obs.train and x_obs.test have different numbers of columns!  The number of features for training and test cases should be exactly the same!");
  }
  # Index data: Put the input data into the right form
  # Convert IDs into numeric indices and
  # Convert some data frames into matrices
  # Index training data
  if (length(src.dst.same)!=1) stop("When calling fit.bst, src.dst.same must be a scalar boolean!");
  if (src.dst.same!=0 && src.dst.same!=1) stop("When calling fit.bst, src.dst.name must be boolean!");
  data.train = indexData(
                obs=obs.train, src.dst.same=src.dst.same, rm.self.link=control$rm.self.link,
                x_obs=x_obs.train, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx,
                add.intercept=control$add.intercept
		);
  # Index test data
  if (!is.null(obs.test)) {
    data.test = indexTestData(
              data.train=data.train, obs=obs.test,
    		  x_obs=x_obs.test, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx
		);
  } else {
    data.test = NULL;
  }
  
  # Model Settings
  if (is.null(control$has.gamma)) {
     control$has.gamma = FALSE;
     if (is.null(x_src) && is.null(x_dst) && !is.null(x_ctx)) control$has.gamma = TRUE;
  }
  if (length(nFactors)==1) nFactors = rep(nFactors, length(model.name));
  if (length(control$has.gamma)==1) control$has.gamma = rep(control$has.gamma,length(model.name));
  if (length(is.logistic)==1) is.logistic = rep(is.logistic,length(model.name));
  if (length(model.name)!=length(nFactors)) stop("When calling fit.bst, model.name and nFactors must have the same length.");
  if (length(model.name)!=length(control$has.gamma)) stop("When calling fit.bst, model.name and control$has.gamma must have the same length.");
  if (length(model.name)!=length(is.logistic)) stop("When calling fit.bst, model.name and is.logistic must have the same length");
  for (i in 1:length(is.logistic))
  {
	if (is.logistic[i]!=0 && is.logistic[i]!=1) stop("When calling fit.bst, is.logistic must be boolean!");
	if (is.logistic[i] && length(which(obs.train$y!=0 & obs.train$y!=1))>0) stop("When calling fit.bst, logistic link function should not be used for non-binary training data! Please set is.logistic=FALSE");
	if (is.logistic[i] && length(which(obs.test$y!=0 & obs.test$y!=1))>0) stop("When calling fit.bst, logistic link function should not be used for non-binary test data! Please set is.logistic=FALSE");
  }

  setting = data.frame(
                name          = model.name,
                nFactors      = nFactors, # number of interaction factors
                has.u         = rep(!src.dst.same,length(model.name)), # whether to use u_i' v_j or v_i' v_j
                has.gamma     = control$has.gamma, 
                nLocalFactors = rep(0,length(model.name)), # just set to 0
                is.logistic   = is.logistic  # whether to use the logistic model for binary rating
  );
  if (is.null(control$nBurnin)) nBurnin = floor(nSamplesPerIter*0.1) else nBurnin=control$nBurnin;

  if (!is.null(control$reg.algo)) {
     if (control$reg.algo=="GLMNet") {
     	source(sprintf("%ssrc/R/model/GLMNet.R",code.dir));
	reg.algo = GLMNet;
     } 
     if (control$reg.algo=="RandomForest") {
     	source(sprintf("%ssrc/R/model/RandomForest.R",code.dir));
	reg.algo = RandomForest;
     }
  } else {
    reg.algo = NULL;
  }
  init.params = control$init.params;
  ans = run.multicontext(
		  		data.train=data.train, data.test=data.test,
                setting=setting,    # Model setting
                nSamples=nSamplesPerIter,   # Number of samples drawn in each E-step: could be a vector of size nIter.
                nBurnIn=nBurnin,     # Number of burn-in draws before take samples for the E-step: could be a vector of size nIter.
                nIter=nIter,       # Number of EM iterations
                approx.interaction=TRUE, # In prediction, predict E[uv] as E[u]E[v].
                reg.algo=reg.algo,     # The regression algorithm to be used in the M-step (NULL => linear regression)
                reg.control=control$reg.control,  # The control paramter for reg.algo
                # initialization parameters
                var_alpha=init.params$var_alpha, var_beta=init.params$var_beta, var_gamma=init.params$var_gamma,
                var_v=init.params$var_v, var_u=init.params$var_u, var_w=init.params$var_w, var_y=init.params$var_y,
                relative.to.var_y=init.params$relative.to.var_y, var_alpha_global=init.params$var_alpha_global, var_beta_global=init.params$var_beta_global,
                # others
                IDs=data.test$IDs,
                out.level=1,      # out.level=1: Save the factor & parameter values to out.dir/model.last and out.dir/model.minTestLoss
                out.dir=out.dir,  # out.level=2: Save the factor & parameter values of each iteration i to out.dir/model.i
                out.overwrite=TRUE,  # whether to overwrite the output directory if it exists
                debug=0,      # Set to 0 to disable internal sanity checking; Set to 100 for most detailed sanity checking
                verbose=1,    # Set to 0 to disable console output; Set to 100 to print everything to the console
                verbose.M=2,
                rm.factors.without.obs.in.loglik=TRUE,
                ridge.lambda=1, # Add diag(lambda) to X'X in linear regression
                zero.mean=rep(0,0), # zero.mean["alpha"] = TRUE  ->  g = 0, etc.
                fix.var=NULL,       # fix.var[["u"]] = n -> var_u = n (NULL -> default, list() -> fix no var)
                max.nObs.for.b=NULL,# maximum number of observations to be used to fit b
                rnd.seed.init=control$random.seed, rnd.seed.fit=control$random.seed+1
  );
  # Do prediction
  if (!is.null(data.test)) {
     pred.y = list();
     for (i in 1:length(model.name)) {
     	 load(sprintf("%s_%s/model.last",out.dir,model.name[i]));
	 pred.model = predict.multicontext(
	 	    model=list(factor=factor, param=param),
		    obs=data.test$obs, feature=data.test$feature, is.logistic=is.logistic[i]
	 );
	 d = data.frame(y=obs.test$y,pred_y=pred.model$pred.y);
	 write.table(d,sprintf("%s_%s/prediction",out.dir,model.name[i]),row.names=F,col.names=T,quote=F,sep="\t");
	 pred.y = c(pred.y, list(pred.model$pred.y));
     }
     names(pred.y) = model.name;
     ans$pred.y = pred.y;
  }
  ans
}

load.code <- function(code.dir)
{
  dyn.load(sprintf("%slib/c_funcs.so",code.dir));
  source(sprintf("%ssrc/R/c_funcs.R",code.dir));
  source(sprintf("%ssrc/R/util.R",code.dir));
  source(sprintf("%ssrc/R/model/util.R",code.dir));
  source(sprintf("%ssrc/R/model/multicontext_model_utils.R",code.dir));
  source(sprintf("%ssrc/R/model/multicontext_model_MStep.R",code.dir));
  source(sprintf("%ssrc/R/model/multicontext_model_EM.R",code.dir));
}

fit.bst.control <- function (
	rm.self.link = FALSE, # Allow Self link?
	add.intercept = TRUE, # Whether to add intercept to each feature matrix
	has.gamma = NULL, # Whether to include context main effect into the model, can be a vector, but length must be equal to the number of model names
	reg.algo=NULL,     # The regression algorithm to be used in the M-step (NULL => linear regression), "GLMNet", or "RandomForest"
	reg.control=NULL,  # The control paramter for reg.algo
	nBurnin = NULL, # Default is 10% of the nSamplesPerIter
	init.params = list(var_alpha=1, var_beta=1, var_gamma=1,
                var_u=1, var_v=1, var_w=1, var_y=NULL,
                relative.to.var_y=FALSE, var_alpha_global=1, var_beta_global=1), # Initial params for all variance components
	random.seed = 0, # The random seed
	...
) {
  if (length(rm.self.link)!=1) stop("When calling fit.bst, rm.self.link in control=fit.bst.control(..., rm.self.link, ...) must be a scalar boolean!");
  if (length(add.intercept)!=1) stop("When calling fit.bst, add.intercept in control=fit.bst.control(..., add.intercept, ...) must be a scalar boolean!");
  if (rm.self.link!=0 && rm.self.link!=1) stop("When calling fit.bst, rm.self.link in control=fit.bst.control(..., rm.self.link, ...) must be boolean!");
  if (add.intercept!=0 && add.intercept!=1) stop("When calling fit.bst, add.intercept in control=fit.bst.control(..., add.intercept, ...) must be boolean!");
  if (!is.null(has.gamma)) {
    for (i in 1:length(has.gamma)) 
    {
    	if (has.gamma[i]!=0 && has.gamma[i]!=1) stop("When calling fit.bst, has.gamma in control=fit.bst.control(..., has.gamma, ...) must be boolean!");
    }
  }
  if (!is.null(reg.algo)) {
     if (reg.algo!="GLMNet" && reg.algo!="RandomForest") stop("When calling fit.bst, reg.algo in control=fit.bst.control(..., reg.algo, ...) must be NULL, 'GLMNet', or 'RandomForest'.  Make sure they are strings if not NULL.");
  }
  if (!is.null(nBurnin)) {
     if (nBurnin<0 || floor(nBurnin)!=nBurnin || length(nBurnin)>1) stop("When calling fit.bst, nBurnin in control=fit.bst.control(..., nBurnin, ...) must be a positive integer");
  }
  list(rm.self.link=rm.self.link,add.intercept=add.intercept, has.gamma=has.gamma, reg.algo=reg.algo, reg.control=reg.control, nBurnin=nBurnin, init.params=init.params, random.seed=random.seed)
}

predict.bst <- function(
	model.file,
	obs.test, # The testing response data
	x_obs.test = NULL, # The data of testing observation features
	x_src = NULL, # The data of context features for source nodes
	x_dst = NULL, # The data of context features for destination nodes
	x_ctx = NULL, # The data of context features for edges	
	code.dir = "" # The top-level directory of where code get installed, "" if you are in that directory
){  
	library(Matrix);
	if (code.dir!="") code.dir = sprintf("%s/",code.dir);
	# Load all the required libraries and source code  
	if (class(try(load.code(code.dir)))=="try-error") stop("When calling predict.bst, code.dir is not specified correctly. code.dir must be the path of the top-level directory of this package (which contains src/ as a sub-directory and file LICENSE)."); 
	
	if(!file.exists(model.file)) stop("When calling predict.bst, the file specified by model.file='",model.file,"' does not exist.  Please specify an existing model file.");
	if("factor" %in% ls()) rm(factor);
	if("param"  %in% ls()) rm(param);
	if("data.train" %in% ls()) rm(data.train);
	load(model.file);
	if(!all(c("factor", "param", "data.train") %in% ls())) stop("When calling predict.bst, some problem with model.file='",model.file,"': The file is not a model file or is corrupted.");
	data.test = indexTestData(
			data.train=data.train, obs=obs.test,
			x_obs=x_obs.test, x_src=x_src, x_dst=x_dst, x_ctx=x_ctx
	);
	
	pred = predict.multicontext(
			model=list(factor=factor, param=param), 
			obs=data.test$obs, feature=data.test$feature, is.logistic=param$is.logistic
	);

	return(pred);
}

