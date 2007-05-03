library(flexmix)

set.seed(123)

klassen <- rep(1:2, c(800, 200))

z<-rep(0,1000)
z[1:800] <- rbinom(800, 3, prob=0.9)
z[801:1000] <- rbinom(200, 3, prob=0.1)
z <- cbind(z, 3-z)

ex1 <- flexmix(z~1, k=2, model= FLXMRglm(family="binomial"))
table(wahr=klassen, flexmix=cluster(ex1))

y <- runif(1000)
z<-rep(0,1000)
z[1:800] <- rbinom(800, 3, prob=y[1:800])
z[801:1000] <- rbinom(200, 3, prob=1-y[801:1000])
z <- cbind(z, 3-z)

ex2 <- flexmix(z~y, k=2, model= FLXMRglm(family="binomial"))
table(wahr=klassen, flexmix=cluster(ex2))


