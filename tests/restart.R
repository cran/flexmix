library("flexmix")
data(Nclus)
mycont=new("FLXcontrol", iter.max=1)

set.seed(123)
ex0 <- flexmix(Nclus ~ 1, k = 4, model = FLXMCmvnorm())

set.seed(123)
ex1 <- flexmix(Nclus ~ 1, k = 4, model = FLXMCmvnorm(), control=mycont)
ex2 <- flexmix(Nclus ~ 1, cluster=posterior(ex1), model = FLXMCmvnorm())

stopifnot(all.equal(ex0@size, ex2@size))
stopifnot(ex0@iter-1==ex2@iter)

ex3a <- flexmix(Nclus ~ 1, cluster=cluster(ex1), model = FLXMCmvnorm())
ex3b <- flexmix(Nclus ~ 1, cluster=cluster(ex1), model = FLXMCmvnorm())

stopifnot(all.equal(ex3a, ex3b))

###**********************************************************

## for one cluster a grouping variable should have no effect on the model
## fit:

data(NPreg)

ex4a <- flexmix(yp ~ x | id1, data = NPreg, k = 1,
                model = FLXMRglm(family = "poisson"))

ex4b <- flexmix(yp ~ x, data = NPreg, k = 1,
                 model = FLXMRglm(family = "poisson"))

stopifnot(all.equal(logLik(ex4a)[1],logLik(ex4b)[1]))

###**********************************************************

## fit with an observation which has a very small likelihood for each of the components
## -> log likelihood would be equal to -Inf

data(NPreg)
NPregNoise <- data.frame(x = c(rep(NPreg$x, 50), 5),
                         yn = c(rep(NPreg$yn, 50), 400))
ex5 <- flexmix(yn ~ x+I(x^2), data = NPregNoise, k = 2)
posterior(ex5, newdata = data.frame(x = 5, yn = 400))


