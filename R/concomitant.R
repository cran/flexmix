FLXPmultinom <- function(formula=~1) {
  z <- new("FLXPmultinom", name="FLXPmultinom", formula=formula)
  multinom.fit <- function(x, y, w, ...) {
    require(nnet)
    r <- ncol(x)
    p <- ncol(y)
    if (p < 2) stop("Multinom requires at least two components.")
    mask <- c(rep(0, r + 1), rep(c(0, rep(1, r)), p - 1))
    nnet.default(x, y, w, mask = mask, size = 0, 
                 skip = TRUE, softmax = TRUE, censored = FALSE, 
                 rang = 0, trace=FALSE,...)
  }
  z@fit <- function(x, y, w, ...) multinom.fit(x,y, w, ...)$fitted.values
  z@refit <- function(x, y, w, ...) {
    if (missing(w) || is.null(w)) w <- rep(1, nrow(y))
    fit <- multinom.fit(x, y, w, ...)
    fit$coefnames <- colnames(x)
    fit$weights <- w
    fit$vcoefnames <- fit$coefnames[1:ncol(x)]
    fit$lab <- 1:ncol(y)
    class(fit) <- c("multinom", "nnet")
    fit$Hessian <- nnet:::multinomHess(fit, x)
    Xr <- qr(x)$rank
    edf <- (ncol(y) - 1) * Xr
    fit$df.residual <- sum(w) - edf
    class(fit) <- NULL
    fit
  }
  z
}

FLXPconstant <- function() {
  new("FLXP", name="FLXPconstant", formula = ~1,
      fit = function(x, y, w, ...){
        if (missing(w) || is.null(w)) return(matrix(colMeans(y), ncol=ncol(y)))
        else return(matrix(colMeans(w*y)/mean(w), ncol=ncol(y)))
      },
      refit = function(x, y, w, ...) {
        if (missing(w) || is.null(w)) return(colMeans(y))
        else return(colMeans(w*y)/mean(w))
      })
}

###**********************************************************

setMethod("FLXgetModelmatrix", signature(model="FLXP"), function(model, data, groups, lhs, ...) {
  mt <- terms(model@formula, data=data)
  mf <- model.frame(delete.response(mt), data=data, na.action = NULL)
  X <- model.matrix(mt, data=mf)
  if (!checkGroup(X, groups$group)) stop("model variables have to be constant for grouping variable")
  model@x <- if (nrow(X)) X[groups$groupfirst,,drop=FALSE] else matrix(1, nrow=sum(groups$groupfirst))
  model
})

checkGroup <- function(x, group) {
  check <- TRUE
  for(g in levels(group)){
    gok <- group==g
    if(any(gok)){
      check <- all(c(check, apply(x[gok,,drop=FALSE], 2, function(z) length(unique(z)) == 1)))
    }
  }
  check
}

###**********************************************************

setMethod("refit", signature(object="FLXP", newdata="missing"),
function(object, newdata, posterior, group, ...) {
  groupfirst <- if (length(group)) groupFirst(group) else rep(TRUE, nrow(posterior))
  new("FLXRP", fitted = object@refit(object@x, posterior[groupfirst,,drop=FALSE], ...))
})

setMethod("refit", signature(object="FLXPmultinom", newdata="missing"),
function(object, newdata, posterior, group, ...) {
  z <- callNextMethod(object, newdata, posterior, group, ...)
  new("FLXRPmultinom", z)
})

setMethod("summary", signature(object="FLXRPmultinom"), function(object) {
  z <- object
  class(object@fitted) <- c("multinom", "nnet")
  fitted.summary <- summary(object@fitted)
  k <- nrow(coef(fitted.summary))+1
  coefs <- lapply(2:k, function(n) {
    coef.p <- fitted.summary$coefficients[n-1,,drop=FALSE]
    s.err <- fitted.summary$standard.errors[n-1,,drop=FALSE]
    tvalue <- coef.p/s.err
    pvalue <- 2 * pnorm(-abs(tvalue))
    coef.table <- t(rbind(coef.p, s.err, tvalue, pvalue))
    dimnames(coef.table) <- list(colnames(coef.p), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    new("Coefmat", coef.table)
  })
  names(coefs) <- paste("Comp", 2:k, sep=".")
  z@summary <- coefs
  z
})

###**********************************************************


fillConcomitant <- function(concomitant, ...) concomitant

setMethod("fillConcomitant", signature(concomitant="FLXPmultinom"), function(concomitant, posterior, weights) {
  concomitant@coef <- rbind("1" = 0, getS3method("coef", "multinom")(concomitant@refit(concomitant@x, posterior, weights)))
  concomitant
})


###**********************************************************

