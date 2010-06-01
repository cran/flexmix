#
#  Copyright (C) 2004-2008 Friedrich Leisch and Bettina Gruen
#  $Id: infocrit.R 4556 2010-05-14 13:20:36Z gruen $
#

setGeneric("nobs", function(object, ...) standardGeneric("nobs"))

setMethod("nobs", signature(object="flexmix"),
function(object, ...) {          
  if (is.null(object@weights)) nrow(object@posterior$scaled) else  sum(object@weights)
})

setMethod("logLik", signature(object="flexmix"),
function(object, ...){
    z <- object@logLik
    attr(z, "df") <- object@df
    attr(z, "nobs") <- nobs(object)
    class(z) <- "logLik"
    z
})

setMethod("AIC", signature(object="flexmix"),
function(object, ..., k=2){
    -2 * object@logLik + object@df * k
})

setMethod("BIC", signature(object="flexmix"),
function(object, ...){
    -2 * object@logLik + object@df * log(nobs(object))
})

setMethod("ICL", signature(object="flexmix"),
function(object, ...){
    -2 * clogLik(object) + object@df * log(nobs(object))
})

setMethod("clogLik", signature(object="flexmix"),
function(object, ...){
    first <- if (length(object@group)) groupFirst(object@group) else TRUE
    post <- object@posterior$unscaled[first,,drop=FALSE]
    n <- nrow(post)
    sum(log(post[seq_len(n) + (clusters(object)[first] - 1)*n]))
})

setMethod("EIC", signature(object="flexmix"),
function(object, ...) {
    first <- if (length(object@group)) groupFirst(object@group) else TRUE
    post <- object@posterior$scaled[first,,drop=FALSE]
    n <- nrow(post)
    lpost <- log(post)
    if (any(is.infinite(lpost))) lpost[is.infinite(lpost)] <- -10^3
    1 + sum(post * lpost)/(n * log(object@k))
})
