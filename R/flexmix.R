#
#  Copyright (C) 2004-2005 Friedrich Leisch and Bettina Gruen
#  $Id: flexmix.R 1943 2005-12-21 10:36:16Z leisch $
#

setClass("FLXcontrol",
         representation(iter.max="numeric",
                        minprior="numeric",
                        tolerance="numeric",
                        verbose="numeric",
                        classify="character",
                        nrep="numeric"),
         prototype(iter.max=200,
                   minprior=0.05,
                   tolerance=10e-7,
                   verbose=0,
                   classify="auto",
                   nrep=1))

setAs("list", "FLXcontrol",
function(from, to){
    z = list2object(from, to)
    z@classify = match.arg(z@classify,
                           c("auto", "weighted", "hard", "random"))
    z
})

setAs("NULL", "FLXcontrol",
function(from, to){
    new(to)
})


###**********************************************************

setClass("FLXmodel",
         representation(fit="function",
                        weighted="logical",
                        name="character",
                        formula="formula",
                        fullformula="formula",
                        x="matrix",
                        y="matrix",
                        preproc.x="function",
                        preproc.y="function"),
         prototype(formula=.~.,
                   fullformula=.~.,
                   preproc.x = function(x) x,
                   preproc.y = function(x) x))

setMethod("show", "FLXmodel",
function(object){
    cat("FlexMix model of type", object@name,"\n\nformula: ")
    print(object@formula)
    cat("Weighted likelihood possible:", object@weighted,"\n\n")
    if(nrow(object@x)>0){
        cat("Regressors:\n")
        print(summary(object@x))
    }
    if(nrow(object@y)>0){
        cat("Response:\n")
        print(summary(object@y))
    }
    cat("\n")
})

setClass("FLXcomponent",
         representation(df="numeric",
                        logLik="function",
                        parameters="list",
                        predict="function"))

setMethod("show", "FLXcomponent",
function(object){
    if(length(object@parameters)>0)
        print(object@parameters)
})
    


###**********************************************************

setClass("FLXdist",
         representation(model="ANY",
                        prior="ANY",
                        components="list",
                        formula="formula",
                        call="call",
                        k="integer"),
         validity=function(object) {
           (object@k == length(object@prior))
         },
         prototype(formula=.~.))

setClass("flexmix",
         representation(posterior="ANY",
                        iter="numeric",
                        cluster="integer",
                        logLik="numeric",
                        df="numeric",
                        control="FLXcontrol",
                        group="factor",
                        size="integer",
                        converged="logical"),
         prototype(group=(factor(integer(0))),
                   formula=.~.),
         contains="FLXdist")

setMethod("show", "flexmix",
function(object){
    cat("\nCall:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\nCluster sizes:\n")
#    print(as.table(TAB))
    print(object@size)
    cat("\n")
    if(!object@converged) cat("no ")
    cat("convergence after", object@iter, "iterations\n")
})


###**********************************************************

setClass("summary.flexmix",
         representation(call="call",
                        AIC="numeric",
                        BIC="numeric",
                        logLik="logLik",
                        comptab="ANY"))

setGeneric("summary")

setMethod("summary", "flexmix",
function(object, eps=1e-4, ...){    
    z <- new("summary.flexmix",
             call = object@call,
             AIC = AIC(object),
             BIC = BIC(object),
             logLik = logLik(object))

    TAB <- data.frame(prior=object@prior,
                      size=object@size)
    rownames(TAB) <- paste("Comp.", 1:nrow(TAB), sep="")
    TAB[["post>0"]] <- colSums(object@posterior$scaled > eps)
    TAB[["ratio"]] <- TAB[["size"]]/TAB[["post>0"]]
    
    z@comptab = TAB
    z
    
})

setMethod("show", "summary.flexmix",
function(object){
    cat("\nCall:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\n")
    print(object@comptab, digits=3)
    cat("\n")
    print(object@logLik)
    cat("AIC:", object@AIC, "  BIC:", object@BIC, "\n")
    cat("\n")    
})


###**********************************************************


setGeneric("flexmix",
           function(formula, data=list(), k=NULL, cluster=NULL,
                    model=NULL, control=NULL)
           standardGeneric("flexmix"))


setMethod("flexmix",
          signature(formula = "formula", model="missing"),
function(formula, data=list(), k=NULL, cluster=NULL,
         model=NULL, control=NULL)
{
    mycall = match.call()
    z <- flexmix(formula=formula, data=data, k=k, cluster=cluster,
                 model=list(FLXglm()), control=control)
    z@call <- mycall
    z
})


setMethod("flexmix",
          signature(formula = "formula", model="list"),
function(formula, data=list(), k=NULL, cluster=NULL,
         model=NULL, control=NULL)
{
    mycall = match.call()
    control = as(control, "FLXcontrol")
    if(is(model, "FLXmodel")) model = list(model)

    groups <- .group(formula, data)
    model <- lapply(model, .model, groups$formula, data)
    postunscaled <- initPosteriors(k, cluster, nrow(model[[1]]@x))
    
    z = FLXfit(model=model, control=control, 
      postunscaled=postunscaled, group=groups$group)
    z@formula = formula
    z@call = mycall
    z
})

setMethod("flexmix",
          signature(formula = "formula", model="FLXmodel"),
function(formula, data=list(), k=NULL, cluster=NULL,
         model=NULL, control=NULL)
{
    mycall = match.call()
    z <- flexmix(formula=formula, data=data, k=k, cluster=cluster,
                 model=list(model), control=control)
    z@call <- mycall
    z
})

###**********************************************************
setGeneric("FLXfit", function(model, control, postunscaled=NULL, group) standardGeneric("FLXfit"))

setMethod("FLXfit", signature(model="list"),
function(model, control, postunscaled=NULL, group)
{
    N = nrow(postunscaled)
    k = ncol(postunscaled)
    
    control <- allweighted(model, control)
    if(control@verbose>0)
      cat("Classification:", control@classify, "\n")
    
    if(length(group)>0){
        groupfirst <- groupFirst(group)
        postunscaled <- groupPosteriors(postunscaled, group)
    }
    else{
        groupfirst <- rep(TRUE, N)
    }

    postscaled <- postunscaled/rowSums(postunscaled)    
    prior <- colMeans(postscaled[groupfirst,,drop=FALSE])

    llh <- -Inf
    converged <- FALSE
        
    components <- list()
    for(K in 1:k) components[[K]] <- list()

    for(iter in 1:control@iter.max){
        postunscaled <- matrix(0, nrow=N, ncol=k)

        if(control@classify != "weighted")
            ok = .FLXgetOK(postscaled, control)
        
        for(m in 1:k){
            for(n in 1:length(model)){
                
                if(control@classify == "weighted"){
                    components[[m]][[n]] <-
                        model[[n]]@fit(model[[n]]@x, model[[n]]@y,
                                        postscaled[,m])
                }
                else{
                    components[[m]][[n]] <-
                        model[[n]]@fit(model[[n]]@x, model[[n]]@y,
                                        ok[,m])
                }
                postunscaled[,m] <- postunscaled[,m] +
                    components[[m]][[n]]@logLik(model[[n]]@x,
                                                  model[[n]]@y)
            }
        }

        if(length(group)>0)
            postunscaled <- groupPosteriors(postunscaled, group)

        postunscaled <- exp(postunscaled)
        ##<FIXME>: wenn eine beobachtung in allen Komonenten extrem
        ## kleine postunscaled-werte hat, ist exp(-postunscaled)
        ## numerisch Null, und damit postscaled Inf
        postunscaled[postunscaled<.Machine$double.eps] <- .Machine$double.eps
        ##</FIXME>
        
        for(m in 1:k)
            postunscaled[,m] <- prior[m] * postunscaled[,m]

        postscaled <- postunscaled/rowSums(postunscaled)

        prior <- colMeans(postscaled[groupfirst,,drop=FALSE])
        
        llh.old <- llh
        llh <- sum(log(rowSums(postunscaled[groupfirst,,drop=FALSE])))
        if(is.na(llh))
            stop(paste(formatC(iter, width=4),
                       "Log-likelihood:", llh))
        if(abs(llh-llh.old)/(abs(llh)+0.1) < control@tolerance){
            if(control@verbose>0){
                printIter(iter, llh)
                cat("converged\n")
            }
            converged <- TRUE
            break
        }
        if(control@verbose && (iter%%control@verbose==0))
            printIter(iter, llh)
        if(any(prior < control@minprior)){
            nok <- which(prior < control@minprior)
            if(control@verbose>0)
                cat("*** Removing",length(nok), "component(s) ***\n")
            prior <- prior[-nok]
            prior <- prior/sum(prior)
            postscaled <- postscaled[,-nok,drop=FALSE]
            components <- components[-nok]
            k <- length(prior)
        }
    }

    df <- k-1   # for the prior probabilities
    for(m in 1:k){
        df <- df+sum(sapply(components[[m]],
                            function(x) x@df))
    }

    names(components) <- paste("Comp", 1:k, sep=".")
    cluster <- max.col(postscaled)
    size <-  tabulate(cluster, nbins=k)
    names(size) <- 1:k
    
    retval <- new("flexmix", model=model, prior=prior,
                  posterior=list(scaled=postscaled,
                                 unscaled=postunscaled),
                  iter=iter, cluster=cluster, size = size,
                  logLik=llh, components=components,
                  control=control, df=df, group=group, k=as(k, "integer"),
                  converged=converged)
                     
    retval
})

.FLXgetOK = function(p, control){

    n = ncol(p)
    N = 1:n
    if(control@classify=="weighted")
        return(matrix(TRUE, nrow=nrow(p), ncol=n))
    else{
        z = matrix(FALSE,  nrow=nrow(p), ncol=n)

        if(control@classify=="hard")
            m = max.col(p)
        else if(control@classify=="random")
            m = apply(p, 1, function(x){sample(N, size=1, prob=x)})
        else
            stop("Unknown classification method")

        z[cbind(1:nrow(p), m)] = TRUE
    }
    z   
}    

.group <- function(formula, data) {
  group = factor(integer(0))
  lf = length(formula)
  formula1 = formula
  if(length(formula[[lf]])>1 &&
     deparse(formula[[lf]][[1]]) == "|"){
    group = as.factor(eval(formula[[lf]][[3]], data))
    formula1[[lf]] = formula[[lf]][[2]]
  }
  return(list(group=group, formula=formula1))
}

setGeneric(".model", function(model, formula, data, rhs=TRUE) standardGeneric(".model"))

setMethod(".model", signature(model="FLXmodel"),
function(model, formula, data, rhs=TRUE)
{
  if(is.null(model@formula))
    model@formula = formula
  
  ## model@fullformula = update.formula(formula, model@formula)
  ## <FIXME>: ist das der richtige weg, wenn ein punkt in beiden
  ## formeln ist?
  model@fullformula = update(terms(formula, data=data), model@formula)
  ## </FIXME>
  
  if (rhs) {
    mf <- model.frame(model@fullformula, data=data)
    model@x = model.matrix(attr(mf, "terms"), data=mf)
    model@y = as.matrix(model.response(mf))
  }
  else {
    mt1 <- terms(model@fullformula, data=data)
    mf <- model.frame(delete.response(mt1), data=data)
    mt <- attr(mf, "terms")
    ## <FIXME>: warum war das da???
    ## attr(mt, "intercept") <- attr(mt1, "intercept")
    ## </FIXME>
    model@x <- model.matrix(mt, data=mf)
  }
  model@x <- model@preproc.x(model@x)
  model@y <- model@preproc.y(model@y)
  model
})

## groupfirst: for grouped observation we need to be able to use
## the posterior of each group, but for computational simplicity
## post(un)scaled has N rows (with mutiple identical rows for each
## group). postscaled[groupfirst,] extracts posteriors of each
## group ordered as the appear in the data set.
groupFirst <- function(x)
{
    x <- as.factor(x)
    z <- rep(FALSE, length(x))
    for(g in levels(x)){
        gok <- x==g
        if(any(gok)){
            z[min(which(gok))] <- TRUE
        }
    }
    z
}

## if we have a group variable, set the posterior to the product
## of all density values for that group (=sum in logarithm)
groupPosteriors <- function(x, group)
{    
    for(g in levels(group)){
        gok <- group==g
        if(any(gok)){
            x[gok,] <- matrix(colSums(x[gok,,drop=FALSE]),
                              nrow=sum(gok), ncol=ncol(x), byrow=TRUE)
        }
    }
    x
}

allweighted <- function(model, control) {
  allweighted = all(sapply(model, function(x) x@weighted))
  if(allweighted){
    if(control@classify=="auto")
      control@classify="weighted"
  }
  else{
    if(control@classify %in% c("auto", "weighted"))
      control@classify="hard"
  }
  control
}

initPosteriors <- function(k, cluster, N) {
  if(is(cluster, "matrix")){
    postunscaled <- cluster
    if (!is.null(k)) if (k != ncol(postunscaled)) stop("specified k does not match the number of columns of cluster")
  }
  else{
    if(is.null(cluster)){
      if(is.null(k))
        stop("either k or cluster must be specified")
      else
        cluster <- sample(1:k, size=N, replace=TRUE)
    }
    else{
      cluster <- as(cluster, "integer")
      if (!is.null(k)) if (k != max(cluster)) stop("specified k does not match the values in cluster")
      k <- max(cluster)
    }
    postunscaled <- matrix(0.1, nrow=N, ncol=k)
    for(K in 1:k){
      postunscaled[cluster==K, K] <- 0.9
    }
  }
  postunscaled
}


###**********************************************************


setGeneric("predict")

setMethod("predict", signature(object="FLXdist"),
function(object, newdata=list(), ...){
    x = list()
    for(m in 1:length(object@model)) {
      comp <- lapply(object@components, "[[", m)
      x[[m]] <- predict(object@model[[m]], newdata, comp, ...)
    }
    z <- list()
    for (k in 1:object@k) {
      z[[k]] <- do.call("cbind", lapply(x, "[[", k))
    }
    names(z) <- paste("Comp", 1:object@k, sep=".")
    z
})

###**********************************************************

setGeneric("parameters",
           function(object, ...) standardGeneric("parameters"))

setMethod("parameters", signature(object="FLXdist"),
function(object, component=1, model=1)
{
    object@components[[component]][[model]]@parameters
})
    
setGeneric("posterior",
           function(object, newdata, ...) standardGeneric("posterior"))

setMethod("posterior", signature(object="flexmix", newdata="missing"),
function(object)
{
    object@posterior$scaled
})

setMethod("posterior", signature(object="FLXdist", newdata="data.frame"),
          function(object, newdata, unscaled=FALSE,...) {
            postunscaled <- matrix(0, nrow = nrow(newdata), ncol = object@k)
            for (m in 1:length(object@model)) {
              comp <- lapply(object@components, "[[", m)
              postunscaled <- postunscaled + posterior(object@model[[m]], newdata, comp, 
                                                      ...)
            }
            if("group" %in% slotNames(object) && length(object@group)>0)
              postunscaled <- groupPosteriors(postunscaled, object@group)
            for(m in 1:object@k)
              postunscaled[,m] <- object@prior[m] * exp(postunscaled[,m])
            if (unscaled) return(postunscaled)
            else return(postunscaled/rowSums(postunscaled))
})            

setMethod("posterior", signature(object="FLXmodel", newdata="data.frame"),
          function(object, newdata, components, ...) {
            mf <- model.frame(terms(object@fullformula, data=newdata), data = newdata)
            x <- model.matrix(attr(mf, "terms"), data = mf)
            y <- as.matrix(model.response(mf))
            sapply(components, function(z) z@logLik(x, y))
})
    
setGeneric("cluster",
           function(object, ...) standardGeneric("cluster"))

setMethod("cluster", signature(object="flexmix"),
function(object)
{
    object@cluster
})
    

###**********************************************************



stepFlexmix <- function(..., K=NULL, nrep=3,
                        compare=c("logLik", "BIC", "AIC"), verbose=TRUE)
{

    compare <- match.arg(compare)
    COMPFUN <- get(compare, mode="function")
    MYCALL <- match.call()
    
    bestFlexmix <- function(...)
    {
        z = new("flexmix", logLik=-Inf)
        for(m in 1:nrep){
            if(verbose) cat(" *")
            x = flexmix(...)
            if(compare=="logLik"){
                if(logLik(x) > logLik(z))
                    z = x
            }
            else{
                if(COMPFUN(x) < COMPFUN(z))
                    z = x
            }
        }
        z
    }

    if(is.null(K)){
        z = bestFlexmix(...)
        z@call <- MYCALL
        if(verbose) cat("\n")
        return(z)
    }
    
    K = as.integer(K)
    if(length(K)==0)
        return(list())
    
    z = list()
    for(n in 1:length(K)){
        if(verbose) cat(K[n], ":")
        z[[as.character(K[n])]] = bestFlexmix(..., k=K[n])
        z[[as.character(K[n])]]@call <- MYCALL
        if(verbose) cat("\n")
    }
    if(length(K)==1)
        return(z[[1]])
    
    z
}

            
