#
#  Copyright (C) 2004-2005 Friedrich Leisch
#  $Id: infocrit.R 1664 2005-06-13 06:11:03Z leisch $
#

setGeneric("logLik")
setGeneric("AIC")
setGeneric("BIC")

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
