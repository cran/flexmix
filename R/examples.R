#
#  Copyright (C) 2004-2005 Friedrich Leisch
#  $Id: examples.R 1896 2005-10-31 12:37:07Z leisch $
#

ExNPreg = function(n)
{
    if(n %% 2 != 0) stop("n must be even")
    
    x <- runif(2*n, 0, 10)
    mp <- exp(c(2-0.2*x[1:n], 1+0.1*x[(n+1):(2*n)]))
    mb <- binomial()$linkinv(c(x[1:n]-5, 5-x[(n+1):(2*n)]))

    data.frame(x=x,
               yn=c(5*x[1:n], 40-(x[(n+1):(2*n)]-5)^2)+3*rnorm(n),
               yp=rpois(2*n, mp),
               yb=rbinom(2*n, size=1, prob=mb),
               class = rep(1:2, c(n,n)),
               id1 = factor(rep(1:n, rep(2, n))),
               id2 = factor(rep(1:(n/2), rep(4, n/2))))
}


    
ExNclus = function(n=100)
{
    if(n %% 2 != 0) stop("n must be even")

    require("mvtnorm")
    
    rbind(rmvnorm(n, mean=rep(0,2)),
          rmvnorm(n, mean=c(8,0), sigma=diag(1:2)),
          rmvnorm(1.5*n, mean=c(-2,6), sigma=diag(2:1)),
          rmvnorm(2*n, mean=c(4,4), sigma=matrix(c(1,.9,.9,1), 2)))
}

    

    



    
