#
#  Copyright (C) 2004-2008 Friedrich Leisch and Bettina Gruen
#  $Id: infocrit.R 3913 2008-03-13 15:13:55Z gruen $
#

setGeneric("nobs", function(object, ...) standardGeneric("nobs"))

setMethod("nobs", signature(object="flexmix"),
function(object, ...) {          
  if (is.null(object@weights)) {
    n <- if (length(object@group)) sum(groupFirst(object@group)) else nrow(object@posterior$scaled)
  }
  else {
    n <- if (length(object@group)) sum(groupFirst(object@group) * object@weights) else sum(object@weights)
  }
  n
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
    sum(log(post[1:n + (cluster(object)[first] - 1)*n]))
})
