.list2object = function(from, to){
    n = names(from)
    s = slotNames(to)
    p = pmatch(n, s)
    if(any(is.na(p)))
        stop(paste("\nInvalid slot name(s) for class",
                   to, ":", paste(n[is.na(p)], collapse=" ")))
    names(from) = s[p]
    do.call("new", c(from, Class=to))
}

###**********************************************************

setClass("FLXcontrol",
         representation(iter.max="numeric",
                        minprior="numeric",
                        tolerance="numeric",
                        verbose="numeric",
                        classify="character"),
         prototype(iter.max=200,
                   minprior=0.05,
                   tolerance=10e-7,
                   verbose=10,
                   classify="auto"))

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
                        x="matrix",
                        y="matrix"))

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
                        call="call"))


setMethod("show", "flexmix",
function(object){
    cat("\nCall:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\nCluster sizes:")
    print(table(object@cluster))
    cat("\nComponent priors:\n")
    print(object@prior)
    cat("\n")
    print(logLik(object))
    cat("AIC:", AIC(object), "  BIC:", BIC(object), "\n")
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
    
    for(n in seq(1,length(model))){
        if(is.null(model[[n]]@formula))
            model[[n]]@formula = formula
        
        famform = update.formula(formula, model[[n]]@formula)

        mf <- model.frame(famform, data=data)
        model[[n]]@x = model.matrix(attr(mf, "terms"), data=mf)
        model[[n]]@y = as.matrix(model.response(mf))
    }
    
    z = FLXfit(model=model, control=control, k=k, cluster=cluster)
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


FLXfit <- function(k=NULL, cluster=NULL, model, control)
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
    retval <- new("flexmix", model=model, prior=prior,
                   posterior=list(scaled=postscaled,
                                  unscaled=postunscaled),
                   iter=iter,
                   cluster=apply(postscaled,1,which.max),
                   logLik=llh, components=components, formula=.~.,
                  control=control, df=df)

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

.FLXprintLogLik = function(iter, logLik)
    cat(formatC(iter, width=4),
        "Log-likelihood:", formatC(logLik, width=12, format="f"),"\n")
    

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
    

###**********************************************************



