FLXMCmvpois <- function(formula=.~.)
{
  z <- new("FLXMC", weighted=TRUE, formula=formula,
                dist="mvpois", name="model-based Poisson clustering")
  
  require(stats)
  z@preproc.y <- function(x){
    x <- as.matrix(x)
    storage.mode(x) <- "integer"
    x
  }

  z@defineComponent <- expression({
    logLik <- function(x, y){
      colSums(dpois(t(y), lambda, log=TRUE))
    }

    predict <- function(x, ...){
      matrix(lambda, nrow = nrow(x), ncol=length(lambda),
             byrow=TRUE)
    }
    
    new("FLXcomponent", parameters=list(lambda=lambda), df=df,
        logLik=logLik, predict=predict)
  })
  z@fit <- function(x, y, w){
    lambda <- colSums(w*y)/sum(w)
    df <- 1
    eval(z@defineComponent)
  }
  z 
}

