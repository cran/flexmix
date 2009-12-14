#
#  Copyright (C) 2004-2008 Friedrich Leisch and Bettina Gruen
#  $Id: stepFlexmix.R 4411 2009-09-23 15:03:19Z gruen $
#

setClass("stepFlexmix",
         representation(models="list",
                        k="integer",
                        nrep="integer",
                        logLiks="matrix",
                        call="call"))



stepFlexmix <- function(..., k=NULL, nrep=3, verbose=TRUE, drop=TRUE,
                        unique=FALSE)
{
    MYCALL <- match.call()
    MYCALL1 <- MYCALL
    
    bestFlexmix <- function(...)
    {
        z = new("flexmix", logLik=-Inf)
        logLiks = rep(NA, length.out = nrep)
        for(m in 1:nrep){
            if(verbose) cat(" *")
            x = try(flexmix(...))
            if (!is(x, "try-error")) {
              logLiks[m] <- logLik(x)
              if(logLik(x) > logLik(z))
                z = x
            }
        }
        return(list(z = z, logLiks = logLiks))
    }

    z = list()
    if(is.null(k)){
        RET = bestFlexmix(...)
        z[[1]] <- RET$z
        logLiks <- as.matrix(RET$logLiks)
        z[[1]]@call <- MYCALL
        z[[1]]@control@nrep <- nrep
        names(z) <- as.character(z[[1]]@k)
        if(verbose) cat("\n")
    }
    else{
        k = as.integer(k)
        logLiks <- matrix(nrow = length(k), ncol = nrep)
        for(n in 1:length(k)){
            ns <- as.character(k[n])
            if(verbose) cat(k[n], ":")
            RET <- bestFlexmix(..., k=k[n])
            z[[ns]] = RET$z
            logLiks[n,] <- RET$logLiks
            MYCALL1[["k"]] <- as.numeric(k[n])
            z[[ns]]@call <- MYCALL1
            z[[ns]]@control@nrep <- nrep
            if(verbose) cat("\n")
        }
    }
    
    z <- z[is.finite(sapply(z, logLik))]
    if (!length(z)) stop("no convergence to a suitable mixture")
    
    if(drop & (length(z)==1)){
        return(z[[1]])
    }
    else{
        z <- return(new("stepFlexmix",
                        models=z,
                        k=as.integer(names(z)),
                        nrep=as.integer(nrep),
                        logLiks=logLiks,
                        call=MYCALL))
        if(unique) z <- unique(z)
        return(z)
    }
}

###**********************************************************

setMethod("unique", "stepFlexmix",
function(x, incomparables=FALSE, ...)
{
    z <- list()
    K <- sapply(x@models, function(x) x@k)
    logLiks <- x@logLiks
    
    for(k in sort(unique(K))){
        n <- which(k==K)
        if(length(n)>1){
            l <- sapply(x@models[n], logLik)
            z[as.character(k)] <- x@models[n][which.max(l)]
            is.na(logLiks[n[-which.max(l)],]) <- TRUE
          }
        else
            z[as.character(k)] <- x@models[n]
    }
    logLiks <- na.omit(logLiks)
    rownames(logLiks) <- names(z)
    attr(logLiks, "na.action") <- NULL
    mycall <- x@call
    mycall["unique"] <- TRUE
    
    return(new("stepFlexmix",
               models=z,
               k=as.integer(names(z)),
               nrep=x@nrep,
               logLiks=logLiks,
               call=mycall))
})

              

###**********************************************************

setMethod("show", "stepFlexmix",
function(object)
{
    cat("\nCall:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\n")
    
    z <- data.frame(iter = sapply(object@models, function(x) x@iter),
                    converged = sapply(object@models, function(x) x@converged),
                    k = sapply(object@models, function(x) x@k),
                    k0 = sapply(object@models, function(x) x@k0),
                    logLik = sapply(object@models, function(x) logLik(x)),
                    AIC = AIC(object),
                    BIC = BIC(object),
                    ICL = ICL(object))
    
    print(z, na.string="")
})

setMethod("AIC", "stepFlexmix",
function(object, ..., k = 2)
{
   sapply(object@models, function(x) AIC(x))
})

setMethod("BIC", "stepFlexmix",
function(object, ...)
{
   sapply(object@models, function(x) BIC(x, ...))
})

setMethod("ICL", "stepFlexmix",
function(object, ...)
{
   sapply(object@models, function(x) ICL(x, ...))
})

setMethod("EIC", "stepFlexmix",
function(object, ...)
{
   sapply(object@models, function(x) EIC(x, ...))
})

###**********************************************************

setMethod("getModel", "stepFlexmix",
function(object, which="BIC")
{
    if(which=="AIC")
        which <- which.min(sapply(object@models, function(x) AIC(x)))
    
    if(which=="BIC")
        which <- which.min(sapply(object@models, function(x) BIC(x)))
    
    if(which=="ICL")
        which <- which.min(sapply(object@models, function(x) ICL(x)))
    
    object@models[[which]]
}
)

###**********************************************************

setMethod("plot", signature(x="stepFlexmix", y="missing"),
function(x, y, what=c("AIC", "BIC", "ICL"), xlab=NULL, ylab=NULL,
         legend="topright", ...)
{
    X <- x@k
    Y <- NULL
    for(w in what){
        Y <- cbind(Y, do.call(w, list(object=x)))
    }
            
    if(is.null(xlab))
        xlab <- "number of components"
    
    if(is.null(ylab)){
        if(length(what)==1)
            ylab <- what
        else
            ylab <- ""
    }
    
    matplot(X, Y, xlab=xlab, ylab=ylab, type="b", lty=1,
            pch=1:length(what), ...)

    if(legend!=FALSE && length(what)>1)
        legend(x=legend, legend=what,
               pch=1:length(what),
               col=1:length(what))

    for(n in 1:ncol(Y)){
        m <- which.min(Y[,n])
        points(X[m], Y[m,n], pch=16, cex=1.5, col=n)
    }
        
    
})
