\name{BIC}
\docType{genericFunction}
\alias{BIC}
\alias{BIC,logLik-method}
\title{Bayesian Information Criterion}
\description{
  This generic function calculates the Bayesian information criterion,
  also known as Schwarz's Bayesian criterion (SBC), for one or several
  fitted model objects for which a log-likelihood value can be obtained,
  according to the formula \eqn{-2 \mbox{log-likelihood} + n_{par}
    \log(n_{obs})}{-2*log-likelihood + npar*log(nobs)}, where
  \eqn{n_{par}}{npar}  represents the
  number of parameters and \eqn{n_{obs}}{nobs} the number of
  observations in the fitted model. 
}
\usage{
BIC(object, \dots)
}
\arguments{
  \item{object}{An object of a suitable class for the BIC to be
    calculated - usually a \code{\link[base:logLik]{logLik}} object
    created by a call to the \code{\link[base:logLik]{logLik}} generic.
  }
  \item{\dots}{Some methods for this generic function may take
    additional, optional arguments.  At present none do.}
}
\value{
  if just one object is provided, returns a numeric value with the
  corresponding BIC; if more than one object are provided, returns a
  \code{data.frame} with rows corresponding to the objects and columns
  representing the number of parameters in the model (\code{df}) and the
  BIC.
}
\references{
  Schwarz, G. (1978) "Estimating the Dimension of a Model", Annals of
  Statistics, 6, 461-464.
}
\author{Function copied from package lme4 by Douglas Bates and Saikat DebRoy}
%\note{}
\seealso{\code{\link[base:logLik]{logLik}}, \code{\link{AIC}}}
%\examples{}
\keyword{models}
