FLXglm <- function(formula=.~.,
                   family=c("gaussian", "binomial", "poisson", "Gamma"))
{
    family <- match.arg(family)
    z <- new("FLXmodel", weighted=TRUE, formula=formula,
             name=paste("FLXglm", family, sep=":"))

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

FLXmclust <- function(formula=.~., diagonal=TRUE)
{
    require(mvtnorm)
    
    z <- new("FLXmodel", weighted=TRUE, formula=formula,
             name="model-based clustering")

    z@fit <- function(x, y, w){

        para <- cov.wt(y, wt=w)[c("center","cov")]
        df <- 2*ncol(x)
        if(diagonal){
            para$cov <- diag(diag(para$cov))
            df <- ncol(x) + ncol(x)^2
        }
        
        predict <- function(x){
            matrix(para$center, nrow=nrow(x), ncol=length(para$center),
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

