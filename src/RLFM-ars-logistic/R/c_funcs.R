### Copyright (c) 2011, Yahoo! Inc.  All rights reserved.
### Copyrights licensed under the New BSD License. See the accompanying LICENSE file for terms.
###
### Author: Liang Zhang

MC_EStep_logistic_arscid.C <- function(
    nSamples, nBurnin,
    user, item, y, x, w, z,
    alpha, beta, u, v,
    ars_XI_alpha, ars_XI_beta,
    ars_XI_u, ars_XI_v,
    b, g0, G, d0, D,
    var_alpha, var_beta, var_u, var_v=rep(1,nFactors),
    ars_ninit, ars_qcent, ars_xl, ars_xu, ars_alpha,
    debug=0,
    beta.int=FALSE, center=FALSE,  main.effects=FALSE,
    check.NA=FALSE
){
    if(debug >= 1) check.input.logistic(user, item, y, x, w, z, alpha, beta, u, v, b, g0, G, d0, D, var_alpha, var_beta, var_u, var_v);
    nObs     = as.integer(length(y));
    nUsers   = as.integer(length(alpha));
    nItems   = as.integer(length(beta));
    nFactors = as.integer(ncol(u));

    xb  = as.double(x %*% b);
    g0w = as.double(w %*% g0);
    if(beta.int==F) d0z = as.double(z %*% d0) else d0z = as.double(cbind(1,z) %*% d0);
    Gw  = as.double(w %*% G);
    Dz  = as.double(z %*% D);

    output = list(
        o.mean       = vector(mode="double",length=nObs),
        alpha.mean   = vector(mode="double",length=nUsers),
        alpha.sumvar = vector(mode="double",length=1),
        beta.mean    = vector(mode="double",length=nItems),
        beta.sumvar  = vector(mode="double",length=1),
        u.mean       = matrix(0.0, nrow=nUsers, ncol=nFactors),
        u.sumvar     = vector(mode="double",length=nFactors),
        v.mean       = matrix(0.0, nrow=nItems, ncol=nFactors),
        v.sumvar     = vector(mode="double",length=nFactors),
        ars_XI_alpha = as.double(ars_XI_alpha),
        ars_XI_beta  = as.double(ars_XI_beta),
        ars_XI_u = as.double(ars_XI_u),
        ars_XI_v = as.double(ars_XI_v)
    );
    if(!is.integer(user)) stop("!is.integer(user)");
    if(!is.integer(item)) stop("!is.integer(item)");
    if(!is.double(y)) stop("!is.double(y)");
    if(!is.double(xb)) stop("!is.double(xb)");
    if(!is.double(g0w)) stop("!is.double(g0w)");
    if(!is.double(d0z)) { cat(d0z,"\n"); stop("!is.double(d0z)"); }
    if(!is.double(Gw)) stop("!is.double(Gw)");
    if(!is.double(Dz)) stop("!is.double(Dz)");
    if(!is.double(alpha)) stop("!is.double(alpha)");
    if(!is.double(beta)) stop("!is.double(beta)");
    if(!is.double(u)) stop("!is.double(u)");
    if(!is.double(v)) stop("!is.double(v)");
    if(!is.double(var_alpha)) stop("!is.double(var_alpha)");
    if(!is.double(var_beta)) stop("!is.double(var_beta)");
    if(!is.double(var_u)) stop("!is.double(var_u)");
    if(!is.double(var_v)) stop("!is.double(var_v)");
    if(!is.logical(main.effects)) stop("!is.logical(main.effects)");
    if(!is.logical(center)) stop("!is.logical(center)");
    if(!is.logical(beta.int)) stop("is.logical(beta.int)");

    if(check.NA){
        if(any(is.na(y))) stop("is.na(y)");
        if(any(is.na(xb))) stop("is.na(xb)");
        if(any(is.na(g0w))) stop("is.na(g0w)");
        if(any(is.na(d0z))) stop("is.na(d0z)");
        if(any(is.na(Gw))) stop("is.na(Gw)");
        if(any(is.na(Dz))) stop("is.na(Dz)");
        if(any(is.na(alpha))) stop("is.na(alpha)");
        if(any(is.na(beta))) stop("is.na(beta)");
        if(any(is.na(u))) stop("is.na(u)");
        if(any(is.na(v))) stop("is.na(v)");
        if(any(is.na(var_alpha))) stop("is.na(var_alpha)");
        if(any(is.na(var_beta))) stop("is.na(var_beta)");
        if(any(is.na(var_u))) stop("is.na(var_u)");
        if(any(is.na(var_v))) stop("is.na(var_v)");
        if(any(is.na(ars_XI_alpha))) stop("is.na(ars_XI_alpha)");
        if(any(is.na(ars_XI_beta))) stop("is.na(ars_XI_beta)");
        if(any(is.na(ars_XI_u))) stop("is.na(ars_XI_u)");
        if(any(is.na(ars_XI_v))) stop("is.na(ars_XI_v)");
    }
        
    .Call("MCEM_EStep_arscid_Call",
        # OUTPUT
        output$o.mean,
        output$alpha.mean,   output$alpha.sumvar,
        output$beta.mean,    output$beta.sumvar,
        output$u.mean,       output$u.sumvar,
        output$v.mean,       output$v.sumvar,
        output$ars_XI_alpha, output$ars_XI_beta,
        output$ars_XI_u,     output$ars_XI_v,
        # INPUT
        as.integer(nSamples), as.integer(nBurnin),
        user, item, y, xb, g0w, d0z, Gw, Dz,
        alpha, beta, u, v,
        var_alpha, var_beta, var_u, var_v,
        nObs, nUsers, nItems, nFactors,
        as.integer(ars_ninit), as.double(ars_qcent),
        as.double(ars_xl), as.double(ars_xu), as.double(ars_alpha),
        # OTHER
        as.integer(debug),
        as.integer(main.effects),
        as.integer(beta.int),
        as.integer(center)
    );
        
    #browser()
    return(output);
}