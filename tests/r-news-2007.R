###################################################
### Changes:
###   o use upper right corner for the legend
###   o construct components plot by hand
###   o cluster -> clusters
###################################################

###################################################
### chunk number 1: 
###################################################
library("flexmix") 
WIDTH=32
options(width=WIDTH, digits=6, prompt = "R> ")
data("whiskey")


###################################################
### chunk number 2: 
###################################################
library("flexmix") 
data("whiskey")
set.seed(1802)


###################################################
### chunk number 3: whiskey
###################################################
Key <- list(corner = c(1,0.9), col = 1:2, rectangle = TRUE, 
                                      text = list(lab=levels(whiskey_brands$Type)))
Col <- as.numeric(factor(whiskey_brands$Type))
print(barchart(rev(colSums(whiskey$Incidence*whiskey$Freq)/sum(whiskey$Freq)),
               origin = 0, col = rev(Col), xlab= "Probability", key = Key))


###################################################
### chunk number 4: 
###################################################
wh_mix <- stepFlexmix(Incidence ~ 1, 
  weights = ~ Freq, data = whiskey, 
  model = FLXMCmvbinary(truncated = TRUE),
  control = list(minprior = 0.005),
  k = 1:7, nrep = 3)


###################################################
### chunk number 5: 
###################################################
wh_best <- getModel(wh_mix, "BIC")
DF <- data.frame(Probability = as.vector(parameters(wh_best)),
                 Components = rep(paste("Comp.", 1:wh_best@k, sep = ""), each = nrow(parameters(wh_best))),
                 Brand = gsub("center.", "", rownames(parameters(wh_best))))
print(barchart(factor(Brand, rev(whiskey_brands$Brand)) ~ Probability | Components, data = DF, origin = 0,
               layout = c(5, 1), col = Col, key = Key,
               xlim = c(-0.1, 1.1)))

###################################################
### chunk number 6: 
###################################################
BIC(wh_mix)
wh_best <- getModel(wh_mix, "BIC")
wh_best


###################################################
### chunk number 7: 
###################################################
prior(wh_best)
parameters(wh_best, component=4:5)[1:2,]


###################################################
### chunk number 8: 
###################################################
purchased <- colSums(parameters(wh_best))
names(purchased) <- 1:length(purchased)
ordering <- order(purchased, decreasing = TRUE)


###################################################
### chunk number 9: 
###################################################
options(width=25, digits=6)
data("patent")


###################################################
### chunk number 10: patent
###################################################
oldpar <- par(mar = c(4, 4, 1, 1) + 0.1)
plot(Patents ~ lgRD, data = patent, pch = 1, col = 1)


###################################################
### chunk number 11: 
###################################################
data("patent")
pat_mix <- flexmix(Patents ~ lgRD, 
  k = 3, data = patent, 
  model = FLXMRglm(family = "poisson"),
  concomitant = FLXPmultinom(~RDS))


###################################################
### chunk number 12: 
###################################################
oldpar <- par(mar = c(4, 4, 1, 1) + 0.1)
plot(Patents ~ lgRD, data = patent, pch = clusters(pat_mix), 
     col = clusters(pat_mix)+1)
lgRDv <- seq(-3, 5, by = 0.05)
newdata <- data.frame(lgRD = lgRDv)
patent.pred <- predict(pat_mix, newdata)
sapply(1:3, function(i) lines(lgRDv, patent.pred[[i]], col = i+1))


###################################################
### chunk number 13: 
###################################################
plot(pat_mix, mark = 3)


###################################################
### chunk number 14: 
###################################################
print(plot(pat_mix, mark = 3))


###################################################
### chunk number 15: 
###################################################
options(show.signif.stars = TRUE, width=80)


###################################################
### chunk number 16: 
###################################################
refit(pat_mix)


###################################################
### chunk number 17: 
###################################################
options(width=WIDTH)


###################################################
### chunk number 18: 
###################################################
plot(refit(pat_mix), bycluster = FALSE)


###################################################
### chunk number 19: 
###################################################
print(plot(refit(pat_mix), bycluster = FALSE))


###################################################
### chunk number 20: 
###################################################
Model_2 <- FLXMRglmfix(family = "poisson",
  nested = list(k = c(1,2), 
    formula = ~lgRD))
Post_1 <- posterior(pat_mix)[,c(2,1,3)]
pat_mix2 <- flexmix(Patents ~ 1,
  concomitant = FLXPmultinom(~RDS), 
  data = patent, cluster = Post_1, 
  model = Model_2)
c(M_1 = BIC(pat_mix), M_2 = BIC(pat_mix2))


