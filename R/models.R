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

setClass("FLXglmmodel",
         representation(family="character"),
         contains="FLXmodel")

FLXglm <- function(formula=.~.,
                   family=c("gaussian", "binomial", "poisson", "Gamma"))
{
    family <- match.arg(family)
    z <- new("FLXglmmodel", weighted=TRUE, formula=formula,
             name=paste("FLXglm", family, sep=":"),
             family=family)

    if(family=="gaussian"){
        z@fit <- function(x, y, w){
            fit <- lm.wfit(x, y, w=w)
            sigma <- sqrt(sum(fit$weights * fit$residuals^2 /
                              mean(fit$weights))/ fit$df.residual)
            fit = fit[c("coefficients")]

            predict <- function(x)
                x%*%coef(fit)
            
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
            fit <- glm.fit(x, y, weights=w, family=binomial())
            fit = fit[c("coefficients","family")]

            predict <- function(x)
                fit$family$linkinv(x%*%coef(fit))
            
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
            fit <- glm.fit(x, y, weights=w, family=poisson())
            fit = fit[c("coefficients","family")]
            rm(w)
            predict <- function(x)
                fit$family$linkinv(x%*%coef(fit))
            
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
            fit <- glm.fit(x, y, weights=w, family=Gamma())
            shape <- sum(fit$prior.weights)/fit$deviance
            fit = fit[c("coefficients","family")]
            rm(w)
            logLik <- function(x, y){
                p = fit$family$linkinv(x%*%coef(fit))
                dgamma(y, shape = shape, scale=p/shape, log=TRUE)
            }
            new("FLXcomponent",
                parameters=list(coef=coef(fit), shape=shape),
                logLik=logLik, df=ncol(x)+1)
        }
    }
    else
        error(paste("Unknown family", family))
    
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
    k = length(object@prior)

    z = new("FLXrefit",
            call=sys.call(-1))

    for(n in 1:k){
        z@model[[n]] <- refit(object@model[[model]],
                              weights=object@posterior$scaled[,n])
    }
    names(z@model) <- paste("Comp", 1:k, sep=".")
    z
})
        
setMethod("refit", signature(object="FLXglmmodel"),
function(object, weights, ...)
{
    z = new("FLXrefitglm")
    
    z@fitted =
        glm.fit(object@x,
                object@y,
                family=get(object@family)(),
                weights=weights)
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
        
        predict <- function(x){
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



