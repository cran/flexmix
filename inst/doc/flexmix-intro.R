### R code from vignette source 'flexmix-intro.Rnw'

###################################################
### code chunk number 1: flexmix-intro.Rnw:30-36
###################################################
suppressWarnings(RNGversion("3.5.0"))
set.seed(1504)
options(width=70, prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE)
grDevices::ps.options(family="Times")
library("flexmix")
data("NPreg")


###################################################
### code chunk number 2: flexmix-intro.Rnw:39-40
###################################################
run <- all(sapply(c("ellipse", "mvtnorm"), requireNamespace, quietly = TRUE))


###################################################
### code chunk number 3: flexmix-intro.Rnw:326-330
###################################################
library("flexmix")
data("NPreg")
m1 <- flexmix(yn ~ x + I(x^2), data = NPreg, k = 2)
m1


###################################################
### code chunk number 4: flexmix-intro.Rnw:333-334
###################################################
parameters(m1, component = 1)


###################################################
### code chunk number 5: flexmix-intro.Rnw:337-338
###################################################
parameters(m1, component = 2)


###################################################
### code chunk number 6: flexmix-intro.Rnw:343-344
###################################################
table(NPreg$class, clusters(m1))


###################################################
### code chunk number 7: flexmix-intro.Rnw:347-348
###################################################
summary(m1)


###################################################
### code chunk number 8: flexmix-intro.Rnw:364-367
###################################################
par(mfrow=c(1,2))
plot(yn~x, col=class, pch=class, data=NPreg)
plot(yp~x, col=class, pch=class, data=NPreg)


###################################################
### code chunk number 9: flexmix-intro.Rnw:385-386
###################################################
print(plot(m1))


###################################################
### code chunk number 10: flexmix-intro.Rnw:406-408
###################################################
rm1 <- refit(m1)
summary(rm1)


###################################################
### code chunk number 11: flexmix-intro.Rnw:429-430
###################################################
options(width=55)


###################################################
### code chunk number 12: flexmix-intro.Rnw:432-435
###################################################
m2 <- flexmix(yp ~ x, data = NPreg, k = 2, 
  model = FLXMRglm(family = "poisson"))
summary(m2)


###################################################
### code chunk number 13: flexmix-intro.Rnw:437-438
###################################################
options(width=65)


###################################################
### code chunk number 14: flexmix-intro.Rnw:442-443
###################################################
print(plot(m2))


###################################################
### code chunk number 15: flexmix-intro.Rnw:486-489
###################################################
m3 <- flexmix(~ x, data = NPreg, k = 2,
  model=list(FLXMRglm(yn ~ . + I(x^2)), 
    FLXMRglm(yp ~ ., family = "poisson")))


###################################################
### code chunk number 16: flexmix-intro.Rnw:504-505
###################################################
print(plot(m3))


###################################################
### code chunk number 17: flexmix-intro.Rnw:534-536
###################################################
m4 <- flexmix(yn ~ x + I(x^2) | id2, data = NPreg, k = 2)
summary(m4)


###################################################
### code chunk number 18: flexmix-intro.Rnw:552-554
###################################################
m5 <- flexmix(yn ~ x + I(x^2), data = NPreg, k = 2,
  control = list(iter.max = 15, verbose = 3, classify = "hard"))


###################################################
### code chunk number 19: flexmix-intro.Rnw:571-575
###################################################
m6 <- flexmix(yp ~ x + I(x^2), data = NPreg, k = 4,
  control = list(minprior = 0.2))

m6  


###################################################
### code chunk number 20: flexmix-intro.Rnw:585-588
###################################################
m7 <- stepFlexmix(yp ~ x + I(x^2), data = NPreg,
  control = list(verbose = 0), k = 1:5, nrep = 5)



###################################################
### code chunk number 21: flexmix-intro.Rnw:594-595
###################################################
getModel(m7, "BIC")


###################################################
### code chunk number 22: flexmix-intro.Rnw:730-735
###################################################
library("flexmix")
set.seed(1504)
options(width=60)
grDevices::ps.options(family="Times")
source("mymclust.R")


###################################################
### code chunk number 23: chunk1 (eval = FALSE)
###################################################
## data("Nclus")
## m1 <- flexmix(Nclus ~ 1, k = 4, model = mymclust())
## summary(m1)


###################################################
### code chunk number 24: flexmix-intro.Rnw:746-749
###################################################
if (run) {
data("Nclus")
m1 <- flexmix(Nclus ~ 1, k = 4, model = mymclust())
summary(m1)
}


###################################################
### code chunk number 25: chunk2 (eval = FALSE)
###################################################
## m2 <- flexmix(Nclus ~ 1, k = 4, model = mymclust(diagonal = FALSE))
## summary(m2)


###################################################
### code chunk number 26: flexmix-intro.Rnw:764-767
###################################################
if (run) {
m2 <- flexmix(Nclus ~ 1, k = 4, model = mymclust(diagonal = FALSE))
summary(m2)
}


###################################################
### code chunk number 27: flexmix-intro.Rnw:773-783
###################################################
par(mfrow=1:2)
if (run) {  
  plotEll(m1, Nclus)
  plotEll(m2, Nclus)
} else {
  plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
  box()
  plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
  box()
}


###################################################
### code chunk number 28: flexmix-intro.Rnw:822-826
###################################################
SI <- sessionInfo()
pkgs <- paste(sapply(c(SI$otherPkgs, SI$loadedOnly), function(x) 
                     paste("\\\\pkg{", x$Package, "} ", 
                           x$Version, sep = "")), collapse = ", ")
