setClass("FLXMRziglm", contains = "FLXMRglm")

FLXMRziglm <- function(formula = . ~ .,
                     family = c("binomial", "poisson"), ...)
{
  family <- match.arg(family)
  new("FLXMRziglm", FLXMRglm(formula, family, ...),
      name = paste("FLXMRziglm", family, sep=":"))
}

setMethod("FLXgetModelmatrix", signature(model="FLXMRziglm"),
          function(model, data, formula, lhs=TRUE, ...)
{
  model <- callNextMethod(model, data, formula, lhs)
  if (attr(terms(model@fullformula), "intercept")==0)
    stop("please include an intercept")
  new("FLXMRziglm", model)
})

setMethod("FLXremoveComponent", signature(model = "FLXMRziglm"),
          function(model, nok, ...)
{
  if (1 %in% nok) model <- as(model, "FLXMRglm")
  model
})

setMethod("FLXmstep", signature(model = "FLXMRziglm"),
          function(model, weights, ...)
{
  coef <- c(-Inf, rep(0, ncol(model@x)-1))
  names(coef) <- colnames(model@x)
  comp.1 <- with(list(coef = coef, df = 0, offset = NULL,
                      family = model@family),
                 eval(model@defineComponent))
  c(list(comp.1), FLXmstep(as(model, "FLXMRglm"),
                           weights[, -1, drop=FALSE]))
})

setClass("FLXRMRziglm", contains = "FLXRM")

setMethod("refit", signature(object="FLXMRziglm"),
          function(object, weights, ...)
{
  coef <- c(-Inf, rep(0, ncol(object@x)-1))
  names(coef) <- colnames(object@x)
  comp.1 <- new("FLXRMRziglm",
                fitted = list(coef))
  c(list(comp.1), callNextMethod(object, weights[, -1, drop=FALSE]))
})

