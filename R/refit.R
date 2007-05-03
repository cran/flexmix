#
#  Copyright (C) 2004-2006 Friedrich Leisch
#  $Id: refit.R 2600 2006-05-11 15:25:49Z leisch $
#

setMethod("refit", signature(object="flexmix"),
function(object, model=1, which = c("model", "concomitant"), summary=TRUE,...)
{
    z = new("FLXR",
            call=sys.call(-1), k = object@k)
    which <- match.arg(which)
    if (which == "model") {
      z@refit <- refit(object@model[[model]],
                       weights=object@posterior$scaled)
      names(z@refit) <- paste("Comp", 1:object@k, sep=".")
      if (summary) z@refit <- lapply(z@refit, summary)
    }
    else {
      z@refit <- refit(object@concomitant, object@posterior$scaled, object@group, w = object@weights)
      if (summary) z@refit <- summary(z@refit)
    }
    z
})

setMethod("refit", signature(object="FLXM"),
function(object, weights, ...)
{
  z <- list()
  for (k in 1:ncol(weights)) {
    z[[k]] = new("FLXRM")

    z[[k]]@fitted =
        object@fit(object@x,
                     object@y,
                     weights[,k])@parameters
  }
  z
})

setMethod("refit", signature(object="FLXMRglm"),
function(object, weights, ...)
{
  z <- list()
  for (k in 1:ncol(weights)) {
    z[[k]] = new("FLXRMRglm")

    z[[k]]@fitted =
        object@refit(object@x,
                     object@y,
                     weights[,k])
  }
  z
})

setMethod("summary", signature(object="FLXR"),
function(object) {
  ## <fixme> for R 2.5.0
  object@call <- match.call()
  if (!length(object@refit) | typeof(object@refit) == "S4") {
    if (is.null(object@refit@summary)) object@refit <- summary(object@refit)
  }
  else if (any(sapply(object@refit, function(x) is.null(x@summary)))) {
    object@refit <- lapply(object@refit, summary)
  }
  object
})

setMethod("summary", signature(object="FLXRM"), function(object) {
  object@summary <- unlist(object@fitted)
  object
})

setMethod("summary", signature(object="FLXRMRglm"),
function(object)
{
  object@summary <- new("Coefmat", coef(summary.glm(object@fitted)))
  object
})

###**********************************************************

setMethod("fitted", signature(object="flexmix"),
function(object, drop=TRUE, aggregate = FALSE, ...)
{
    x<- list()
    for(m in 1:length(object@model)) {
      comp <- lapply(object@components, "[[", m)
      x[[m]] <- fitted(object@model[[m]], comp, ...)
    }
    if (aggregate) {
      z <- sapply(x, function(z) do.call("cbind", z) %*% object@prior, simplify = drop)
    }
    else {
      z <- list()
      for (k in 1:object@k) {
        z[[k]] <- do.call("cbind", lapply(x, "[[", k))
      }
      names(z) <- paste("Comp", 1:object@k, sep=".")
      if(drop && all(lapply(z, ncol)==1)){
        z <- sapply(z, unlist)
      }
    }
    z
})

setMethod("fitted", signature(object="FLXM"),
function(object, components, ...) {
  lapply(components, function(z) z@predict(object@x))
})

setMethod("fitted", signature(object="FLXRMRglm"),
function(object, ...)
{
    object@fitted$fitted.values
})
        
setMethod("fitted", signature(object="FLXR"),
function(object, ...)
{
  sapply(object@refit, fitted)
})

setMethod("predict", signature(object="FLXM"), function(object, newdata, components, ...)
{
  object <- FLXgetModelmatrix(object, newdata, formula = object@fullformula, lhs = FALSE) 
  z <- list()
  for(k in 1:length(components)){
    z[[k]] <- components[[k]]@predict(object@x, ...)
  }
  z
})
                           
###**********************************************************


