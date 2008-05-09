#
#  Copyright (C) 2004-2008 Friedrich Leisch and Bettina Gruen
#  $Id: refit.R 3969 2008-05-05 07:12:05Z gruen $
#
###*********************************************************

setMethod("getParameters", signature(object="FLXdist"),
function(object, model) {
  if (missing(model)) model <- seq_along(object@model)
  coefficients <- unlist(lapply(model, function(m) {
    Model <- unlist(getParameters(object@model[[m]], lapply(object@components, "[[", m)))
    names(Model) <- paste("model.", m, "_", names(Model), sep = "")
    Model
  }))
  c(coefficients, getParameters(object@concomitant))
})

setMethod("getParameters", signature(object="FLXM"),
function(object, components, ...) {
  lapply(components, function(x) unlist(slot(x, "parameters")))
})

setMethod("getParameters", signature(object="FLXMRglm"),
function(object, components, ...) {
  if (object@family == "Gamma") stop("not implemented yet")
  parms <- lapply(components, function(x) unlist(slot(x, "parameters")))
  Design <- FLXgetDesign(object, components)
  if (object@family == "gaussian") {
    parms <- lapply(parms, function(x) {
      x["sigma"] <- log(x["sigma"])
      x
    })
    colnames(Design) <- gsub("sigma$", "log(sigma)", colnames(Design))
  }
  parms_unique <- vector(length = ncol(Design))
  names(parms_unique) <- colnames(Design)
  for (k in seq_along(parms)) 
    parms_unique[as.logical(Design[k,])] <- parms[[k]]
  parms_unique
})

setMethod("getParameters", signature(object="FLXP"),
function(object, ...) {
  if (length(object@coef) == 1) return(NULL) 
  alpha <- log(object@coef[-1]) - log(object@coef[1])
  names(alpha) <- paste("concomitant", paste("Comp", 2:length(object@coef), "alpha", sep = "."), sep = "_")
  return(alpha)
})  

setMethod("getParameters", signature(object="FLXPmultinom"),
function(object, ...) {
  coefficients <- object@coef[,-1,drop=FALSE]
  if (ncol(coefficients) > 0) {
    Names <- paste("Comp", rep(2:(ncol(coefficients)+1), each = nrow(coefficients)), rownames(coefficients), sep = ".")
    coefficients <- as.vector(coefficients)
    names(coefficients) <- paste("concomitant", Names, sep = "_")
    return(coefficients)
  }else return(NULL)
})  

setMethod("VarianceCovariance", signature(object="flexmix"),
function(object, model = TRUE, gradient, ...) {
  if (object@control@classify != "weighted") stop("Only for weighted ML estimation possible.")
  if (length(getParameters(object)) != object@df) stop("not implemented yet for restricted parameters.")
  if (missing(gradient)) gradient <- FLXgradlogLikfun(object)
  fit <- optim(fn = FLXlogLikfun(object), par = getParameters(object), gr = gradient,
               hessian = TRUE, method = "BFGS", control = list(fnscale = -1))
  list(coef = fit$par, vcov = -solve(as.matrix(fit$hessian)))
})  

setMethod("Likfun_comp", signature(object="flexmix"),
function(object) {
  function(parms) {
    postunscaled <- matrix(0, nrow = FLXgetObs(object@model[[1]]), ncol = object@k)
    for (m in seq_along(object@model))
      postunscaled <- postunscaled + FLXdeterminePostunscaled(object@model[[m]], lapply(object@components, "[[", m))
    if(length(object@group)>0)
      postunscaled <- groupPosteriors(postunscaled, object@group)
    exp(postunscaled)
  }
})

setMethod("FLXlogLikfun", signature(object="flexmix"),
function(object, ...) function(parms) {
  object <- FLXreplaceParameters(object, parms)
  groupfirst <- if (length(object@group) > 1) groupFirst(object@group) else rep(TRUE, FLXgetObs(object@model[[1]]))
  postunscaled <- getPriors(object@concomitant, object@group, groupfirst) * Likfun_comp(object)(parms)
  if (is.null(object@weights)) return(sum(log(rowSums(postunscaled[groupfirst,,drop=FALSE]))))
  else return(sum(log(rowSums(postunscaled[groupfirst,,drop=FALSE]))*object@weights[groupfirst]))
})

setMethod("getPriors", signature(object="FLXP"),
function(object, group, groupfirst) {
  priors <- apply(object@coef, 2, function(x) object@x %*% x)
  ungroupPriors(priors/rowSums(priors), group, groupfirst)
})

setMethod("getPriors", signature(object="FLXPmultinom"),
function(object, group, groupfirst) {
  priors <- apply(object@coef, 2, function(x) exp(object@x %*% x))
  ungroupPriors(priors/rowSums(priors), group, groupfirst)
})

setMethod("FLXreplaceParameters", signature(object="FLXdist"),
function(object, parms) {
  comp_names <- names(object@components)
  components <- list()
  for (m in seq_along(object@model)) {
    indices <- grep(paste("^model.", m, sep = ""), names(parms))
    components[[m]] <- FLXreplaceParameters(object@model[[m]], lapply(object@components, "[[", m), parms[indices])
  }
  object@components <- lapply(seq_along(object@components), function(k) lapply(components, "[[", k))
  names(object@components) <- comp_names
  if (object@k > 1) {
    indices <- grep("^concomitant_", names(parms))
    object@concomitant <- FLXreplaceParameters(object@concomitant, parms[indices])
  }
  object
})

setMethod("FLXreplaceParameters", signature(object="FLXM"),
function(object, component, parms) {
  Parameters <- list()
  for (i in seq_along(component@parameters)) {
    Parameters[[i]] <- parms[1:length(component@parameters[[i]])]
    attributes(Parameters[[i]]) <- attributes(component@parameters[[i]])
    parms <- parms[-c(1:length(component@parameters[[i]]))]
  }
  names(Parameters) <- names(component@parameters)
  Parameters$df <- component@df
  variables <- c("x", "y")
  for (var in variables) 
    assign(var, slot(object, var))
  with(Parameters, eval(object@defineComponent))
})

setMethod("FLXreplaceParameters", signature(object="FLXMRglm"),
function(object, components, parms) {
  Design <- FLXgetDesign(object, components)
  lapply(seq_along(components), function(k) {
    Parameters <- list()
    parms_k <- parms[as.logical(Design[k,])]
    for (i in seq_along(components[[k]]@parameters)) {
      Parameters[[i]] <- parms_k[1:length(components[[k]]@parameters[[i]])]
      attributes(Parameters[[i]]) <- attributes(components[[k]]@parameters[[i]])
      parms_k <- parms_k[-c(1:length(components[[k]]@parameters[[i]]))]
    }
    names(Parameters) <- names(components[[k]]@parameters)
    if (object@family == "gaussian") {
      Parameters[["sigma"]] <- exp(Parameters[["sigma"]])
    }
    Parameters$df <- components[[k]]@df
    variables <- c("x", "y", "offset", "family")
    for (var in variables) {
      assign(var, slot(object, var))
    }
    with(Parameters, eval(object@defineComponent))
  })
})

setMethod("FLXreplaceParameters", signature(object="FLXP"),
function(object, parms) {
  parms <- exp(c(0, parms))
  parms <- parms/sum(parms)
  attributes(parms) <- attributes(object@coef)
  object@coef <- parms
  object
})

setMethod("FLXreplaceParameters", signature(object="FLXPmultinom"),
function(object, parms) {
  parms <- cbind(0, matrix(parms, nrow = nrow(object@coef)))
  attributes(parms) <- attributes(object@coef)
  object@coef <- parms
  object
})

setMethod("FLXgradlogLikfun", signature(object="flexmix"),
function(object, ...) {
  existFunction <- all(sapply(object@model, existGradient))
  if (object@k > 1) existFunction <- c(existFunction,
                                       existGradient(object@concomitant))
  if (any(!existFunction)) return(NULL)
  function(parms) {
    object <- FLXreplaceParameters(object, parms)
    groupfirst <- if (length(object@group) > 1) groupFirst(object@group) else rep(TRUE, FLXgetObs(object@model[[1]]))
    Lik_comp <- Likfun_comp(object)(parms)
    Priors <- getPriors(object@concomitant, object@group, groupfirst)
    Priors_Lik_comp <- Priors * Lik_comp
    Lik_sep <-  rowSums(Priors_Lik_comp)
    weights <- Priors_Lik_comp/Lik_sep
    ModelScores <- lapply(seq_along(object@model), function(m)
                          FLXgradlogLikfun(object@model[[m]],
                                           lapply(object@components, "[[", m), weights))
    ModelScores <- lapply(ModelScores, lapply, groupPosteriors, object@group)
    ConcomitantScores <- if (object@k > 1) {
      conc <- FLXgradlogLikfun(object@concomitant, Priors)
      lapply(conc, function(p) matrix(rowSums(sapply(seq_along(p), function(k)
                                                     p[[k]] * weights[,k])), ncol = ncol(p[[1]])))
    } else NULL
    colSums(cbind(do.call("cbind", lapply(ModelScores, function(x) do.call("cbind", x))),
                  do.call("cbind", ConcomitantScores))[groupfirst,,drop=FALSE])
  }
})

setMethod("existGradient", signature(object = "FLXMRglm"),
function(object) {
  if (object@family == "Gamma") return(FALSE)
  TRUE
})

setMethod("existGradient", signature(object = "FLXMRglmfix"),
function(object) FALSE)

setMethod("existGradient", signature(object = "FLXP"),
function(object) TRUE)

setMethod("FLXgradlogLikfun", signature(object="FLXMRglm"),
function(object, components, weights, ...) {
  lapply(seq_along(components), function(k) {
    res <- if (object@family == "binomial") as.vector(object@y[,1] - rowSums(object@y)*components[[k]]@predict(object@x))
    else as.vector(object@y - components[[k]]@predict(object@x))
    Scores <- weights[,k] * res * object@x
    if (object@family == "gaussian") {
      Scores <- cbind(Scores/components[[k]]@parameters$sigma^2, weights[,k]
                      * (-1 + res^2/components[[k]]@parameters$sigma^2))
    }
    Scores
  })
})                       

setMethod("FLXgradlogLikfun", signature(object="FLXP"),
function(object, fitted, ...) {
  Pi <- apply(fitted, 1, function(x) as.vector((diag(nrow = length(x))[,-1]) - rep(x[-1], each = length(x))))
  Pi <- lapply(1:(ncol(fitted)-1), function(i) t(Pi[(i-1)*ncol(fitted) + 1:ncol(fitted),]))
  lapply(Pi, function(p) lapply(1:ncol(p), function(i) apply(object@x, 2, "*", p[,i])))
})

setMethod("refit", signature(object = "flexmix"),
function(object, newdata, method = c("optim", "mstep"), ...) {
  method <- match.arg(method)
  if (method == "optim") {
    VarCov <- VarianceCovariance(object, ...)
    z <- new("FLXRoptim",
             call=sys.call(-1), k = object@k,
             coef = VarCov$coef, vcov = VarCov$vcov)
    z@components <- lapply(seq_along(object@model), function(m) {
      begin_name <- paste("^model", m, sep = ".")
      indices <- grep(begin_name, names(z@coef))
      refit_optim(object@model[[m]], components = lapply(object@components, "[[", m), coef = z@coef[indices], se = sqrt(diag(z@vcov)[indices]))
    })
    z@concomitant <- if (object@k > 1) {
      indices <- grep("^concomitant_", names(z@coef))
      refit_optim(object@concomitant, coef = z@coef[indices], se = sqrt(diag(z@vcov)[indices]))
    } else NULL
  } else {
    z <- new("FLXRmstep",
             call=sys.call(-1), k = object@k)
    z@components <- lapply(object@model, function(x) {
      x <- refit_mstep(x, weights=object@posterior$scaled)
      names(x) <- paste("Comp", 1:object@k, sep=".")
      x
    })
    z@concomitant <- if (object@k > 1) refit_mstep(object@concomitant, posterior = object@posterior$scaled,
                                                       group = object@group, w = object@weights) else NULL
  }
  z
})

setMethod("refit_optim", signature(object = "FLXM"),
function(object, components, coef, se) {
  design <- FLXgetDesign(object, components)
  x <- lapply(1:nrow(design), function(k) {
    rval <- cbind(Estimate = coef[as.logical(design[k,])],
                  "Std. Error" = se[as.logical(design[k,])])
    pars <- components[[k]]@parameters$coef
    rval <- rval[1:length(pars),,drop=FALSE]
    rownames(rval) <- names(pars)
    zval <- rval[,1]/rval[,2]
    new("Coefmat", cbind(rval, "z value" = zval, "Pr(>|z|)" = 2 * pnorm(abs(zval), lower.tail = FALSE)))
  })
  names(x) <- paste("Comp", seq_along(x), sep = ".")
  x
})

setMethod("refit_optim", signature(object = "FLXP"),
function(object, coef, se) {
  x <- lapply(2:ncol(object@coef), function(k) {
    indices <- grep(paste("Comp", k, sep = "."), names(coef))
    rval <- cbind(Estimate = coef[indices],
                  "Std. Error" = se[indices])
    rval <- rval[1:nrow(object@coef),,drop=FALSE]
    rownames(rval) <- rownames(object@coef)
    zval <- rval[,1]/rval[,2]
    new("Coefmat", cbind(rval, "z value" = zval, "Pr(>|z|)" = 2 * pnorm(abs(zval), lower.tail = FALSE)))
  })
  names(x) <- paste("Comp", 1 + seq_along(x), sep = ".")
  x
})

setMethod("FLXgetDesign", signature(object = "FLXM"),
function(object, components, ...) {
  parms <- lapply(components, function(x) unlist(slot(x, "parameters")))
  nr_parms <- sapply(parms, length)
  cumSum <- cumsum(c(0, nr_parms))
  Design <- t(sapply(1:(length(cumSum)-1), function(i) rep(c(0, 1, 0), c(cumSum[i], nr_parms[i], max(cumSum) - cumSum[i] - nr_parms[i]))))
  colnames(Design) <- paste(rep(paste("Comp", 1:nrow(Design), sep = "."), nr_parms),
                            unlist(lapply(parms, names)), sep = "_")
  Design
})

setMethod("FLXgetDesign", signature(object = "FLXMRglmfix"),
function(object, components, ...) {
  Design <- object@design
  if (object@family == "gaussian") {
    cumSum <- cumsum(c(0, object@variance))
    variance <- matrix(sapply(1:(length(cumSum)-1), function(i)
                              rep(c(0, 1, 0), c(cumSum[i], object@variance[i], length(components) - cumSum[i] - object@variance[i]))),
                       nrow = length(components))
    colnames(variance) <- paste("Comp", apply(variance, 2, function(x) which(x == 1)[1]), "sigma", sep= ".")
    Design <- cbind(Design, variance)
  }
  Design
})

###*********************************************************

setMethod("refit_mstep", signature(object="FLXM"),
function(object, newdata, weights, ...)
{
  lapply(seq_len(ncol(weights)), function(k) 
         object@fit(object@x,
                    object@y,
                    weights[,k])@parameters)
})

setMethod("refit_mstep", signature(object="FLXMRglm"),
function(object, newdata, weights, ...)
{
  lapply(seq_len(ncol(weights)), function(k) {
    fit <- object@refit(object@x,
                        object@y,
                        weights[,k])
    fit <- c(fit,
             list(formula = object@fullformula,
                  terms = object@terms,
                  contrasts = object@contrasts,
                  xlevels = object@xlevels))
    class(fit) <- c("glm", "lm")
    fit
  })
})

###**********************************************************

setMethod("fitted", signature(object="flexmix"),
function(object, drop=TRUE, aggregate = FALSE, ...)
{
    x<- list()
    for(m in seq_along(object@model)) {
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
      for (k in seq_len(object@k)) {
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

setMethod("predict", signature(object="FLXM"), function(object, newdata, components, ...)
{
  object <- FLXgetModelmatrix(object, newdata, formula = object@fullformula, lhs = FALSE) 
  z <- list()
  for(k in seq_along(components)) 
    z[[k]] <- components[[k]]@predict(object@x, ...)
  z
})
                           
###**********************************************************

setMethod("Lapply", signature(object="FLXRmstep"), function(object, FUN, model = 1, component = TRUE, ...) {
  X <- object@components[[model]]
  lapply(X[component], FUN, ...)
})

###*********************************************************

setMethod("refit_mstep", signature(object="flexmix", newdata="listOrdata.frame"),
function(object, newdata, ...)
{
    z <- new("FLXR",
            call=sys.call(-1), k = object@k)
    z@components <- lapply(object@model, function(x) {
      x <- refit_mstep(x, newdata = newdata,
                           weights=posterior(object, newdata = newdata))
      names(x) <- paste("Comp", 1:object@k, sep=".")
      x
    })
    z@concomitant <- if (object@k > 1) refit_mstep(object@concomitant, newdata, object@posterior$scaled, object@group, w = object@weights)
    else NULL
    z
})

setMethod("refit_mstep", signature(object="FLXMRglm", newdata="listOrdata.frame"),
function(object, newdata, weights, ...)
{
  w <- weights
  lapply(seq_len(ncol(w)), function(k) {
    newdata$weights <- weights <- w[,k]
    weighted.glm(formula = object@fullformula, data = newdata,
                 family = object@family, weights = weights, ...)
  })
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
