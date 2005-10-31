#
#  Copyright (C) 2004-2005 Friedrich Leisch
#  $Id: models.R 1664 2005-06-13 06:11:03Z leisch $
#

setClass("FLXglmmodel",
         representation(family="character",
                        refit="function"),
         contains="FLXmodel")

FLXglm <- function(formula=.~.,
                   family=c("gaussian", "binomial", "poisson", "Gamma"),
                   offset=NULL)
{
    family <- match.arg(family)

    glmrefit <- function(x, y, w)
        glm.fit(x, y, weights=w, offset=offset,
                family=get(family, mode="function")())
                
    z <- new("FLXglmmodel", weighted=TRUE, formula=formula,
             name=paste("FLXglm", family, sep=":"),
             family=family, refit=glmrefit)
    
    if(family=="gaussian"){

        z@refit <- function(x, y, w) lm.wfit(x, y, w=w, offset=offset)
        
        z@fit <- function(x, y, w){
            fit <- lm.wfit(x, y, w=w, offset=offset)
            sigma <- sqrt(sum(fit$weights * fit$residuals^2 /
                              mean(fit$weights))/ (nrow(x)-fit$rank))
            fit = fit[c("coefficients")]

            predict <- function(x, ...) {
                dotarg = list(...)
                if("offset" %in% names(dotarg)) offset <- dotarg$offset
                p <- x%*%coef(fit)
                if (!is.null(offset)) p <-  p + offset
                p
            }

            logLik <- function(x, y)
                dnorm(y, mean=predict(x), sd=sigma, log=TRUE)
            
            new("FLXcomponent",
                parameters=list(coef=coef(fit), sigma=sigma),
                logLik=logLik,
                predict=predict,
                df=ncol(x)+1)
        }
    }
    else if(family=="binomial"){
        z@fit <- function(x, y, w){
            fit <- glm.fit(x, y, weights=w, family=binomial(), offset=offset)
            fit = fit[c("coefficients","family")]

            predict <- function(x, ...) {
                dotarg = list(...)
                if("offset" %in% names(dotarg)) offset <- dotarg$offset
                p <- x%*%coef(fit)
                if (!is.null(offset)) p <- p + offset
                fit$family$linkinv(p)
            }

            logLik <- function(x, y){
                dbinom(y[,1], size=rowSums(y), prob=predict(x), log=TRUE)
            }
            
            new("FLXcomponent",
                parameters=list(coef=coef(fit)),
                logLik=logLik,
                predict=predict,
                df=ncol(x))
        }
    }
    else if(family=="poisson"){
        z@fit <- function(x, y, w){
            fit <- glm.fit(x, y, weights=w, family=poisson(), offset=offset)
            fit = fit[c("coefficients","family")]
            rm(w)
            
            predict <- function(x, ...) {
                dotarg = list(...)
                if("offset" %in% names(dotarg)) offset <- dotarg$offset
                p <- x%*%coef(fit)
                if (!is.null(offset)) p <- p + offset
                fit$family$linkinv(p)
            }

            logLik <- function(x, y){
                dpois(y, lambda=predict(x), log=TRUE)
            }

            new("FLXcomponent",
                parameters=list(coef=coef(fit)),
                logLik=logLik,
                predict=predict,
                df=ncol(x))
        }
    }
    else if(family=="Gamma"){
        z@fit <- function(x, y, w){
            fit <- glm.fit(x, y, weights=w, family=Gamma(), offset=offset)
            shape <- sum(fit$prior.weights)/fit$deviance
            fit = fit[c("coefficients","family")]
            rm(w)

            predict <- function(x, ...) {
                dotarg = list(...)
                if("offset" %in% names(dotarg)) offset <- dotarg$offset
                p <- x%*%coef(fit)
                if (!is.null(offset)) p <- p + offset
                fit$family$linkinv(p)
            }

            logLik <- function(x, y){
                p = fit$family$linkinv(x%*%coef(fit))
                dgamma(y, shape = shape, scale=p/shape, log=TRUE)
            }
            new("FLXcomponent", predict=predict,
                parameters=list(coef=coef(fit), shape=shape),
                logLik=logLik, df=ncol(x)+1)
        }
    }
    else stop(paste("Unknown family", family))
    z
}

###**********************************************************

setClass("FLXrefit",
         representation(model="list",
                        call="call"))

setClass("FLXrefitglm",
         representation(fitted="list"))

setGeneric("refit", function(object, ...) standardGeneric("refit"))


setMethod("refit", signature(object="flexmix"),
function(object, model=1, ...)
{
    z = new("FLXrefit",
            call=sys.call(-1))

    z@model <- refit(object@model[[model]],
                     weights=object@posterior$scaled)

    names(z@model) <- paste("Comp", 1:object@k, sep=".")
    z
})

setMethod("refit", signature(object="FLXglmmodel"),
function(object, weights, ...)
{
  z <- list()
  for (k in 1:ncol(weights)) {
    z[[k]] = new("FLXrefitglm")

    z[[k]]@fitted =
        object@refit(object@x,
                     object@y,
                     weights[,k])
  }
  z
})


setMethod("show", signature(object="FLXrefit"),
function(object)
{
    cat("\nCall:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\nNumber of components:", length(object@model), "\n\n")
})

setMethod("summary", signature(object="FLXrefitglm"),
function(object)
{
    printCoefmat(coef(summary.glm(object@fitted)), signif.stars=FALSE)
})

setMethod("summary", signature(object="FLXrefit"),
function(object)
{
    cat("\nCall:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\n")
    
    k = length(object@model)
    for(n in 1:k){
        cat("Component", n, ":\n")
        summary(object@model[[n]])
        cat(ifelse(k==n, "\n\n", "\n-------------\n"))
    }
})

###**********************************************************

setGeneric("fitted")

setMethod("fitted", signature(object="flexmix"),
function(object, drop=TRUE, ...)
{
    x<- list()
    for(m in 1:length(object@model)) {
      comp <- lapply(object@components, "[[", m)
      x[[m]] <- fitted(object@model[[m]], comp, ...)
    }
    z <- list()
    for (k in 1:object@k) {
      z[[k]] <- do.call("cbind", lapply(x, "[[", k))
    }
    names(z) <- paste("Comp", 1:object@k, sep=".")
    if(drop && all(lapply(z, ncol)==1)){
        z <- sapply(z, unlist)
    }
    z
})

setMethod("fitted", signature(object="FLXmodel"),
function(object, components, ...) {
  lapply(components, function(z) z@predict(object@x))
})

setMethod("fitted", signature(object="FLXrefitglm"),
function(object, ...)
{
    object@fitted$fitted.values
})
        
setMethod("fitted", signature(object="FLXrefit"),
function(object, ...)
{
    sapply(object@model, fitted)
})

setMethod("predict", signature(object="FLXglmmodel"), function(object, newdata, components, ...) {
  mt1 <- terms(object@fullformula)
  mf <- model.frame(delete.response(mt1), data=newdata)
  mt <- attr(mf, "terms")
  attr(mt, "intercept") <- attr(mt1, "intercept")
  x <- model.matrix(mt, data=mf)
  
  z <- list()
  for(k in 1:length(components)){
    z[[k]] <- components[[k]]@predict(x, ...)
  }
  z
})
                           
###**********************************************************

FLXmclust <- function(formula=.~., diagonal=TRUE)
{
    z <- new("FLXmodel", weighted=TRUE, formula=formula,
             name="model-based Gaussian clustering")

    require(mvtnorm)

    z@fit <- function(x, y, w){
        
        para <- cov.wt(y, wt=w)[c("center","cov")]
        df <- (3*ncol(y) + ncol(y)^2)/2
        if(diagonal){
            para$cov <- diag(diag(para$cov))
            df <- 2*ncol(y)
        }
        
        predict <- function(x, ...){
            matrix(para$center, nrow=nrow(y), ncol=length(para$center),
                   byrow=TRUE)
        }
        
        logLik <- function(x, y){
            dmvnorm(y, mean=para$center, sigma=para$cov, log=TRUE)
        }
            
        new("FLXcomponent", parameters=para, df=df,
            logLik=logLik, predict=predict)
    }
    z
}



###**********************************************************
