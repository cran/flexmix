\name{FLXMRlmer}
\alias{FLXMRlmer}
\alias{FLXMRlmer-class}
\alias{FLXMRlmm-class}
\alias{FLXMRlmmfix-class}
\alias{FLXdeterminePostunscaled,FLXMRlmer-method}
\alias{FLXdeterminePostunscaled,FLXMRlmm-method}
\alias{FLXmstep,FLXMRlmer-method}
\alias{FLXmstep,FLXMRlmm-method}
\alias{FLXgetModelmatrix,FLXMRlmer-method}
\alias{FLXgetModelmatrix,FLXMRlmm-method}
\alias{FLXMRlmm}
\alias{FLXgetObs,FLXMRlmm-method}
\alias{FLXmstep,FLXMRlmmfix-method}
\alias{predict,FLXMRlmm-method}
\title{FlexMix Interface to Linear Mixed Models}
\description{
  This is a driver which allows fitting of mixtures of linear models
  with random effects.
}
\usage{
FLXMRlmer(formula = . ~ ., random, weighted = FALSE,
          control = list(), eps = .Machine$double.eps)
FLXMRlmm(formula = . ~ ., random, lm.fit = c("lm.wfit",
          "smooth.spline"), varFix = c(Random = FALSE, Residual =
          FALSE), \dots)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{flexmix} using
    \code{\link{update.formula}}. Default is to use the original
    \code{flexmix} model
    formula.}
  \item{random}{A formula for specifying the random effects.}
  \item{weighted}{A logical indicating if the model can estimate
    weighted ML.}
  \item{control}{A list of control parameters. See
    \code{\link[lme4]{lmer}} for details.}
  \item{eps}{Observations with a component-specific posterior smaller
    than \code{eps} are omitted in the M-step for this component.}
  \item{lm.fit}{A character string indicating if the coefficients should
    be fitted using either a linear model or the function
    \code{smooth.spline}}
  \item{varFix}{Named logical vector of length 2 indicating if the
    variance of the random effects and the residuals are fixed over the
    components.}
  \item{\dots}{Additional arguments to be passed to \code{smooth.spline}.}
}
\details{
  \code{FLXMRlmm} allows only one random effect. \code{FLXMRlmer} allows
  an arbitrary number of random effects if \code{weighted=FALSE}; a
  certain structure of the model matrix of the random effects has to be
  given for weighted ML estimation, i.e. where \code{weighted=TRUE}.
}
\value{
  Returns an object of class \code{FLXMRlmer} and \code{FLXMRlmm}.
}
\section{warning}{
  For \code{FLXMRlmer} the weighted ML estimation is only correct if the
  covariate matrix of the random effects is the same for each
  observation. By default non-weighted ML estimation is made. If this is
  changed the condition on the covariate matrix of the random effects is
  checked.
}
\author{Bettina Gruen}
\examples{
id <- rep(1:100, each = 10)
x <- rep(1:10, 100)
sample <- data.frame(y = rep(rnorm(unique(id)/2, 0, c(5, 2)), each = 10) +
                         rnorm(length(id), rep(c(3, 8), each = 10)) +
                         rep(c(0, 3), each = 10) * x,
                     x = x,
                     id = factor(id))
fitted <- flexmix(.~.|id, k = 2, model = FLXMRlmer(y ~ x, random = ~ 1),
                  data = sample)
parameters(fitted)
xyplot(y ~ x | clusters(fitted), groups = id, data = sample, type = "l")
}
\keyword{models}