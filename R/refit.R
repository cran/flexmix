#
#  Copyright (C) 2004-2006 Friedrich Leisch
#  $Id: refit.R 2600 2006-05-11 15:25:49Z leisch $
#

setMethod("refit", signature(object="flexmix", newdata="missing"),
function(object, newdata, model=1, which = c("model", "concomitant"), summary=TRUE,...)
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
      z@refit <- refit(object@concomitant, posterior = object@posterior$scaled,
                       group = object@group, w = object@weights)
      if (summary) z@refit <- summary(z@refit)
    }
    z
})

setMethod("refit", signature(object="FLXM", newdata="missing"),
function(object, newdata, weights, ...)
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

setMethod("refit", signature(object="FLXMRglm", newdata="missing"),
function(object, newdata, weights, ...)
{
  z <- list()
  for (k in 1:ncol(weights)) {
    fit <- object@refit(object@x,
                        object@y,
                        weights[,k])
    fit <- c(fit,
             list(formula = object@fullformula,
                  terms = object@terms,
                  contrasts = object@contrasts,
                  xlevels = object@xlevels))
    class(fit) <- c("glm", "lm")
    z[[k]] <- new("FLXRMRglm", fitted = fit)
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
  dispersion <- NULL
  if (is.null(object@fitted$dispersion)) {
    if (object@fitted$df.residual > 0) 
      dispersion <- sum((object@fitted$weights *
                         object@fitted$residuals^2)[object@fitted$weights > 0])/(object@fitted$df.residual *
                                                                                 mean(object@fitted$weights))
  }
  object@summary <- new("Coefmat", coef(summary.glm(object@fitted, dispersion = dispersion)))
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
      z <- lapply(x, function(z) matrix(rowSums(matrix(sapply(1:object@k, function(K) z[[K]] * object@prior[K]), ncol = object@k)),
                                        nrow = nrow(z[[1]])))
      if(drop && all(lapply(z, ncol)==1)){
        z <- sapply(z, unlist)
      }
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

setMethod("Lapply", signature(object="FLXR"), function(object, FUN, component = TRUE, ...) {
  X <- object@refit
  if (!is.vector(X) || is.object(X)) 
    X <- as.list(X)
  lapply(X[component], function(x) Lapply(x, FUN, ...))
})

setMethod("Lapply", signature(object="FLXRM"), function(object, FUN, ...) {
  FUN <- match.fun(FUN)
  FUN(object@fitted, ...)
})

###*********************************************************

setMethod("refit", signature(object="flexmix", newdata="listOrdata.frame"),
function(object, newdata, model=1, which = c("model", "concomitant"), summary=TRUE,...)
{
    z = new("FLXR",
            call=sys.call(-1), k = object@k)
    which <- match.arg(which)
    if (which == "model") {
      z@refit <- refit(object@model[[model]], newdata = newdata,
                       weights=posterior(object, newdata = newdata, ...), ...)
      names(z@refit) <- paste("Comp", 1:object@k, sep=".")
      if (summary) z@refit <- lapply(z@refit, summary)
    }
    else {
      z@refit <- refit(object@concomitant, newdata, object@posterior$scaled, object@group, w = object@weights)
      if (summary) z@refit <- summary(z@refit)
    }
    z
})

setMethod("refit", signature(object="FLXMRglm", newdata="listOrdata.frame"),
function(object, newdata, weights, ...)
{
  z <- list()
  w <- weights
  for (k in 1:ncol(w)) {
    newdata$weights <- weights <- w[,k]
    fit <-  weighted.glm(formula = object@fullformula, data = newdata,
                         family = object@family, weights = weights, ...)
    z[[k]] <- new("FLXRMRglm", fitted = fit)
  }
  z
})

weighted.glm <- function(weights, ...) {
  fit <- eval(as.call(c(as.symbol("glm"), c(list(...), list(weights = weights, x = TRUE)))))
  fit$df.null <- sum(weights) + fit$df.null - fit$df.residual - fit$rank
  fit$df.residual <- sum(weights) - fit$rank
  fit$method <- "weighted.glm.fit"
  fit
}

weighted.glm.fit <- function(x, y, weights, offset = NULL, family = "gaussian", ...) {
  if (!is.function(family) & !is(family, "family"))
    family <- get(family, mode = "function", envir = parent.frame())
  fit <- c(glm.fit(x, y, weights = weights, offset=offset,
                   family=family),
           list(call = sys.call(), offset = offset,
                control = eval(formals(glm.fit)$control),            
                method = "weighted.glm.fit"))
  fit$df.null <- sum(weights) + fit$df.null - fit$df.residual - fit$rank
  fit$df.residual <- sum(weights) - fit$rank
  fit$x <- x
  fit
}
