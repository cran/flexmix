mymclust <- function (formula = .~., diagonal = TRUE) 
{    
    require("mvtnorm")
    
    retval <- new("FLXmodel", weighted = TRUE, formula = formula, 
                  name = "my model-based clustering")
    
    retval@fit <- function(x, y, w) {
        
        para <- cov.wt(y, wt = w)[c("center", "cov")]
        df <- (3 * ncol(y) + ncol(y)^2)/2
        
        if (diagonal) {
            para$cov <- diag(diag(para$cov))
            df <- 2 * ncol(y)
        }

        logLik <- function(x, y) {
            dmvnorm(y, mean = para$center, sigma = para$cov, 
                log = TRUE)
        }

        predict <- function(x) {
            matrix(para$center, nrow = nrow(y),
                   ncol = length(para$center), byrow = TRUE)
        }

        new("FLXcomponent", parameters = para, df = df,
            logLik = logLik, predict = predict)
    }
    retval
}
