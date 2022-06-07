%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: FLXMRziglm.Rd 5229 2022-06-02 14:47:26Z gruen $
%
\name{FLXMRziglm}
\alias{FLXMRziglm}
\alias{FLXMRziglm-class}
\alias{refit,FLXMRziglm-method}
\alias{FLXreplaceParameters,FLXMRziglm-method}
\alias{FLXgradlogLikfun,FLXMRziglm-method}
\title{FlexMix Interface to Zero Inflated Generalized Linear Models}
\description{
  This is a driver which allows fitting of zero inflated poisson and
  binomial models.
}
\usage{
FLXMRziglm(formula = . ~ ., family = c("binomial", "poisson"), ...)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{flexmix} using
    \code{\link{update.formula}}. Default is to use the original
    \code{flexmix} model
    formula.}
  \item{family}{A character string naming a \code{\link{glm}}
    family function.}
  \item{\dots}{passed to \code{FLXMRglm}}
}
\value{
  Returns an object of class \code{FLXMRziglm} inheriting from \code{FLXMRglm}.
}
\author{Friedrich Leisch and Bettina Gruen}
\note{
  In fact this only approximates zero inflated models by fixing the
  coefficient of the intercept at -Inf and the other coefficients at
  zero for the first component.
}
\examples{
  data("dmft", package = "flexmix")
  Model <- FLXMRziglm(family = "poisson")
  Fitted <- flexmix(End ~ log(Begin + 0.5) + Gender + Ethnic + Treatment, 
                    model = Model, k = 2 , data = dmft, 
                    control = list(minprior = 0.01))
  summary(refit(Fitted))
}
\keyword{models}
