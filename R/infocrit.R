#
#  Copyright (C) 2004-2005 Friedrich Leisch
#  $Id: infocrit.R 2609 2006-05-15 11:32:16Z gruen $
#

setMethod("logLik", signature(object="flexmix"),
function(object, ...){
    z <- object@logLik
    attr(z, "df") <- object@df
    attr(z, "nobs") <- nrow(object@posterior$scaled)
    class(z) <- "logLik"
    z
})

setMethod("AIC", signature(object="flexmix"),
function(object, ..., k=2){
    -2 * object@logLik + object@df * k
})

setMethod("BIC", signature(object="flexmix"),
function(object, ...){
    -2 * object@logLik + object@df * log(nrow(object@posterior$scaled))
})

setMethod("ICL", signature(object="flexmix"),
function(object, ...){
    -2 * clogLik(object) + object@df * log(nrow(object@posterior$scaled))
})

setMethod("clogLik", signature(object="flexmix"),
function(object, ...){
    n <- nrow(object@posterior$unscaled)
    sum(log(object@posterior$unscaled[1:n + (cluster(object) - 1)*n]))
})
