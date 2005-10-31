library("flexmix")
data(Nclus)
mycont=new("FLXcontrol", iter.max=1)

set.seed(123)
ex0 <- flexmix(Nclus ~ 1, k = 4, model = FLXmclust())

set.seed(123)
ex1 <- flexmix(Nclus ~ 1, k = 4, model = FLXmclust(), control=mycont)
ex2 <- flexmix(Nclus ~ 1, cluster=posterior(ex1), model = FLXmclust())

stopifnot(all.equal(ex0@size, ex2@size))
stopifnot(ex0@iter-1==ex2@iter)

ex3a <- flexmix(Nclus ~ 1, cluster=cluster(ex1), model = FLXmclust())
ex3b <- flexmix(Nclus ~ 1, cluster=cluster(ex1), model = FLXmclust())

stopifnot(all.equal(ex3a, ex3b))

###**********************************************************

## for one cluster a grouping variable should have no effect on the model
## fit:

data(NPreg)

ex4a <- flexmix(yp ~ x | id1, data = NPreg, k = 1,
                model = FLXglm(family = "poisson"))

ex4b <- flexmix(yp ~ x, data = NPreg, k = 1,
                 model = FLXglm(family = "poisson"))

stopifnot(all.equal(logLik(ex4a),logLik(ex4b)))

###**********************************************************

