#
#  FlexMix: Flexible mixture modelling in R
#  Copyright (C) 2004 Friedrich Leisch
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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
    z = .list2object(from, to)
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
                        y="matrix"),
         prototype(formula=.~.,
                   fullformula=.~.))

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

setClass("flexmix",
         representation(model="ANY",
                        prior="ANY",
                        posterior="ANY",
                        iter="numeric",
                        cluster="integer",
                        logLik="numeric",
                        df="numeric",
                        components="list",
                        formula="formula",
                        control="FLXcontrol",
                        call="call",
                        group="factor",
                        k="integer",
                        size="integer",
                        converged="logical"),
         prototype(group=(factor(integer(0))),
                   formula=.~.))


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

    group = factor(integer(0))
    lf = length(formula)
    formula1 = formula
    if(length(formula[[lf]])>1 &&
       deparse(formula[[lf]][[1]]) == "|"){
        group = as.factor(eval(formula[[lf]][[3]], data))
        formula1[[lf]] = formula[[lf]][[2]]
    }
    
    for(n in seq(1,length(model))){
        if(is.null(model[[n]]@formula))
            model[[n]]@formula = formula1
        
        model[[n]]@fullformula = update.formula(formula1, model[[n]]@formula)
        mf <- model.frame(model[[n]]@fullformula, data=data)
        model[[n]]@x = model.matrix(attr(mf, "terms"), data=mf)
        model[[n]]@y = as.matrix(model.response(mf))
    }
    
    z = FLXfit(model=model, control=control, k=k,
               cluster=cluster, group=group)
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


FLXfit <- function(k=NULL, cluster=NULL, model, control, group)
{
    N = nrow(model[[1]]@x)
        
    if(is.null(cluster)){
        if(is.null(k))
            stop("either k or cluster must be specified")
        else
            cluster <- sample(1:k, size=N, replace=TRUE)
    }
    else{
        cluster <- as(cluster, "integer")
        k <- max(cluster)
    }

    allweighted = all(sapply(model, function(x) x@weighted))
    if(allweighted){
        if(control@classify=="auto")
            control@classify="weighted"
    }
    else{
        if(control@classify %in% c("auto", "weighted"))
            control@classify="hard"
    }

    if(control@verbose>0)
        cat("Classification:", control@classify, "\n")
        

    prior <- rep(1/k, k)
    postscaled <- matrix(0.1, nrow=N, ncol=k)
    components <- list()
    for(K in 1:k){
        components[[K]] <- list()
        postscaled[cluster==K, K] <- 0.9
    }
    llh <- -Inf
    converged <- FALSE
    
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

        
        ## if we have a group variable, set the posterior to the product
        ## of all density values for that group (=sum in logarithm)
        if(length(group)>0){
            for(g in levels(group)){
                gok <- group==g
                if(any(gok)){
                    postunscaled[gok,] <-
                        matrix(colSums(postunscaled[gok,,drop=FALSE]),
                               nrow=sum(gok), ncol=k, byrow=TRUE)
                }
            }
        }

        for(m in 1:k)
            postunscaled[,m] <- prior[m] * exp(postunscaled[,m])
        postscaled <- postunscaled/rowSums(postunscaled)

        prior <- colMeans(postscaled)
        
        llh.old <- llh
        llh <- sum(log(rowSums(postunscaled)))
        if(is.na(llh))
            stop(paste(formatC(iter, width=4),
                       "Log-likelihood:", llh))
        if(abs(llh-llh.old)/(abs(llh)+0.1) < control@tolerance){
            if(control@verbose>0){
                .FLXprintLogLik(iter, llh)
                cat("converged\n")
            }
            converged <- TRUE
            break
        }
        if(control@verbose && (iter%%control@verbose==0))
            .FLXprintLogLik(iter, llh)
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
    cluster <- apply(postscaled, 1, which.max)
    size <-  tabulate(cluster)
    names(size) <- 1:k
    
    retval <- new("flexmix", model=model, prior=prior,
                  posterior=list(scaled=postscaled,
                                 unscaled=postunscaled),
                  iter=iter, cluster=cluster, size = size,
                  logLik=llh, components=components,
                  control=control, df=df, group=group, k=as(k, "integer"),
                  converged=converged)
                     
    retval
}

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

.FLXprintLogLik = function(iter, logLik, label="Log-likelihood")
    cat(formatC(iter, width=4),
        label, ":", formatC(logLik, width=12, format="f"),"\n")
    

###**********************************************************


setGeneric("predict")

setMethod("predict", signature(object="flexmix"),
function(object, newdata=list(), ...){

    K = length(object@components)
    N = length(object@model)
    z = list(length=K)
    for(n in 1:N){
        f <- update.formula(object@formula,
                            object@model[[n]]@formula)
        mt1 <- terms(f)
        mf <- model.frame(delete.response(mt1), data=newdata)
        mt <- attr(mf, "terms")
        attr(mt, "intercept") <- attr(mt1, "intercept")
        x <- model.matrix(mt, data=mf)

        for(k in 1:K){
            if(n==1)
                z[[k]] = object@components[[k]][[n]]@predict(x)
            else
                z[[k]] = cbind(z[[k]], object@components[[k]][[n]]@predict(x))
        }
    }
    names(z) <- paste("Comp", 1:K, sep=".")
    z
})

###**********************************************************

setGeneric("parameters",
           function(object, ...) standardGeneric("parameters"))

setMethod("parameters", signature(object="flexmix"),
function(object, component=1, model=1)
{
    object@components[[component]][[model]]@parameters
})
    
setGeneric("posterior",
           function(object, ...) standardGeneric("posterior"))

setMethod("posterior", signature(object="flexmix"),
function(object)
{
    object@posterior$scaled
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

            
