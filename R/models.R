#
#  Copyright (C) 2004-2006 Friedrich Leisch
#  $Id: models.R 3467 2007-04-27 16:05:26Z gruen $
#

FLXMRglm <- function(formula=.~.,
                     family=c("gaussian", "binomial", "poisson", "Gamma"),
                     offset=NULL)
{
    family <- match.arg(family)
    glmrefit <- function(x, y, w)
        glm.fit(x, y, weights=w, offset=offset,
                family=get(family, mode="function")())
                
    z <- new("FLXMRglm", weighted=TRUE, formula=formula,
             name=paste("FLXMRglm", family, sep=":"),
             family=family, refit=glmrefit)
    
    if(family=="gaussian"){
      z@defineComponent <- expression({
        predict <- function(x, ...) {
          dotarg = list(...)
          if("offset" %in% names(dotarg)) offset <- dotarg$offset
          p <- x%*%coef
          if (!is.null(offset)) p <-  p + offset
          p
        }

        logLik <- function(x, y)
          dnorm(y, mean=predict(x), sd=sigma, log=TRUE)

        new("FLXcomponent",
            parameters=list(coef=coef, sigma=sigma),
            logLik=logLik, predict=predict,
            df=df)
      })

      z@fit <- function(x, y, w){
        fit <- lm.wfit(x, y, w=w, offset=offset)
        sigma <- sqrt(sum(fit$weights * fit$residuals^2 /
                          mean(fit$weights))/ (nrow(x)-fit$rank))
        fit = fit[c("coefficients")]
        coef <- coef(fit)
        df <- ncol(x)+1
        eval(z@defineComponent)
      }
    }
    else if(family=="binomial"){
      z@preproc.y <- function(x){
        if (ncol(x) != 2)
          stop("for the binomial family, y must be a 2 column matrix\n",
               "where col 1 is no. successes and col 2 is no. failures")
        x
      }
      
      z@defineComponent <- expression({
        predict <- function(x, ...) {
          dotarg = list(...)
          if("offset" %in% names(dotarg)) offset <- dotarg$offset
          p <- x%*%coef
          if (!is.null(offset)) p <- p + offset
          get(family, mode = "function")()$linkinv(p)
        }
        logLik <- function(x, y)
          dbinom(y[,1], size=rowSums(y), prob=predict(x), log=TRUE)

        new("FLXcomponent",
            parameters=list(coef=coef),
            logLik=logLik, predict=predict,
            df=df)
      })

      z@fit <- function(x, y, w){
        fit <- glm.fit(x, y, weights=w, family=binomial(), offset=offset)
        fit = fit[c("coefficients","family")]
        coef <- coef(fit)
        df <- ncol(x)
        eval(z@defineComponent)
      }
    }
    else if(family=="poisson"){
      z@defineComponent <- expression({
        predict <- function(x, ...) {
          dotarg = list(...)
          if("offset" %in% names(dotarg)) offset <- dotarg$offset
          p <- x%*%coef
          if (!is.null(offset)) p <- p + offset
          get(family, mode = "function")()$linkinv(p)
        }
        logLik <- function(x, y)
          dpois(y, lambda=predict(x), log=TRUE)
        
        new("FLXcomponent",
            parameters=list(coef=coef),
            logLik=logLik, predict=predict,
            df=df)
      })
          
      z@fit <- function(x, y, w){
        fit <- glm.fit(x, y, weights=w, family=poisson(), offset=offset)
        fit = fit[c("coefficients","family")]
        rm(w)
        coef <- coef(fit)
        df <- ncol(x)
        eval(z@defineComponent)
      }
    }
    else if(family=="Gamma"){
      z@defineComponent <- expression({
        predict <- function(x, ...) {
          dotarg = list(...)
          if("offset" %in% names(dotarg)) offset <- dotarg$offset
          p <- x%*%coef
          if (!is.null(offset)) p <- p + offset
          get(family, mode = "function")()$linkinv(p)
        }
        logLik <- function(x, y)
          dgamma(y, shape = shape, scale=predict(x)/shape, log=TRUE)
        
        new("FLXcomponent", 
            parameters = list(coef, shape = shape),
            predict = predict, logLik = logLik,
            df = df)
      })

      z@fit <- function(x, y, w){
        fit <- glm.fit(x, y, weights=w, family=Gamma(), offset=offset)
        shape <- sum(fit$prior.weights)/fit$deviance
        fit = fit[c("coefficients","family")]
        coef <- coef(fit)
        df <- ncol(x)+1
        eval(z@defineComponent)
      }
    }
    else stop(paste("Unknown family", family))
    z
}

###**********************************************************

FLXMCmvnorm <- function(formula=.~., diagonal=TRUE)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             dist = "mvnorm", name="model-based Gaussian clustering")

    require(mvtnorm)
    z@defineComponent <- expression({
      logLik <- function(x, y)
        dmvnorm(y, mean=center, sigma=cov, log=TRUE)
    
      predict <-  function(x, ...)
        matrix(center, nrow=nrow(x), ncol=length(center),
               byrow=TRUE)
      new("FLXcomponent", parameters=list(center = center, cov = cov),
          df=df, logLik=logLik, predict=predict)
    })
    
    z@fit <- function(x, y, w){
      para <- cov.wt(y, wt=w)[c("center","cov")]
      df <- (3*ncol(y) + ncol(y)^2)/2
      if(diagonal){
        para$cov <- diag(diag(para$cov))
        df <- 2*ncol(y)
      }
      with(para, eval(z@defineComponent))
    }
    z
}


###**********************************************************

FLXMCmvbinary <- function(formula=.~., truncated = FALSE) {
  if (truncated) return(MCmvbinary_truncated())
  else return(MCmvbinary())
}

MCmvbinary <- function(formula=.~.)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             dist = "mvbinary", name="model-based binary clustering")

    ## make sure that y is binary
    z@preproc.y <- function(x){
        x <- as.matrix(x)
        storage.mode(x) <- "logical"
        storage.mode(x) <- "integer"
        x
    }
    z@defineComponent <- expression({
      predict <- function(x, ...){
        matrix(center, nrow=nrow(x), ncol=length(center),
               byrow=TRUE)
      }
        
      logLik <- function(x, y){
        p <- matrix(center, nrow=nrow(x), ncol=length(center),
                    byrow=TRUE)
        rowSums(log(y*p+(1-y)*(1-p)))
      }
            
      new("FLXcomponent", parameters=list(center=center), df=df,
          logLik=logLik, predict=predict)
    })

    z@fit <- function(x, y, w){
      center <- colSums(w*y)/sum(w)
      df <- ncol(y)
      eval(z@defineComponent)
    }
    
    z
}


###**********************************************************

binary_truncated <- function(y, w, maxit = 200, epsilon = .Machine$double.eps) {
  r_k <- colSums(y*w)/sum(w)
  r_0 <- 0
  llh.old <- -Inf
  for (i in 1:maxit) {
    p <- r_k/(1+r_0)
    llh <- sum((r_k*log(p))[r_k > 0])+ sum(((1 - r_k + r_0) * log(1-p))[(1-r_k+r_0) > 0])
    if (abs(llh - llh.old)/(abs(llh) + 0.1) < epsilon) break    
    llh.old <- llh
    prod_p <- prod(1-p)
    r_0 <- prod_p/(1-prod_p)
  }
  p
}

MCmvbinary_truncated <- function(formula=.~.)
{
    z <- MCmvbinary(formula=formula)
    z@defineComponent <- expression({
      predict <- function(x, ...) {
        matrix(center, nrow = nrow(x), ncol = length(center), 
               byrow = TRUE)
      }
      logLik <- function(x, y) {
        p <- matrix(center, nrow = nrow(x), ncol = length(center), 
                    byrow = TRUE)
        rowSums(log(y * p + (1 - y) * (1 - p))) - log(1 - prod(1-center))
      }
      new("FLXcomponent", parameters = list(center = center), df = df, 
          logLik = logLik, predict = predict)
    })
    z@fit <- function(x, y, w){
      center <- binary_truncated(y, w)
      df <- ncol(y)
      eval(z@defineComponent)
    }
    
    z
}
