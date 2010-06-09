\name{FLXMRlmmc}
\alias{FLXMRlmmc}
\alias{FLXMRlmmc-class}
\alias{FLXMRlmmcfix-class}
\alias{FLXMRlmc-class}
\alias{FLXMRlmcfix-class}
\alias{FLXmstep,FLXMRlmmc-method}
\alias{FLXmstep,FLXMRlmmcfix-method}
\alias{FLXmstep,FLXMRlmc-method}
\alias{FLXmstep,FLXMRlmcfix-method}
\alias{FLXgetModelmatrix,FLXMRlmmc-method}
\alias{FLXgetModelmatrix,FLXMRlmc-method}
\alias{FLXgetObs,FLXMRlmc-method}
\alias{FLXdeterminePostunscaled,FLXMRlmmc-method}
\alias{FLXdeterminePostunscaled,FLXMRlmc-method}
\alias{predict,FLXMRlmc-method}
\title{FlexMix Interface to Linear Mixed Models with Left-Censored
  Observations}
\description{
  This is a driver which allows fitting of mixtures of linear models
  with random effects and left-censored observations.
}
\usage{
FLXMRlmmc(formula = . ~ ., random, censored, varFix, eps = 10^-6, ...)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{flexmix} using
    \code{\link{update.formula}}. Default is to use the original
    \code{flexmix} model
    formula.}
  \item{random}{A formula for specifying the random effects. If missing
    no random effects are fitted.}
  \item{varFix}{If random effects are specified a named logical vector
    of length 2 indicating if the variance of the random effects and the
    residuals are fixed over the components. Otherwise a logical
    indicating if the variance of the residuals are fixed over the
    components.}
  \item{censored}{A formula for specifying the censoring variable.}
  \item{eps}{Observations with an a-posteriori probability smaller or
    equal to eps are omitted in the M-step.}
  \item{\dots}{Additional arguments to be passed to \code{smooth.spline}.}
}
\value{
  Returns an object of class \code{FLXMRlmmc}, \code{FLXMRlmmcfix},
  \code{FLXMRlmc} or \code{FLXMRlmcfix}.
}
\author{Bettina Gruen}
\keyword{models}
