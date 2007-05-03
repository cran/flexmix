%
%  Copyright (C) 2004-2005 Friedrich Leisch
%  $Id: refit.Rd 3221 2007-01-22 14:27:40Z gruen $
%
\name{refit}
\alias{refit,flexmix-method}
\alias{refit,FLXM-method}
\alias{refit,FLXMRglm-method}
\alias{plot,FLXR,missing-method}
\alias{show,FLXR-method}
\alias{show,FLXR-method}
\alias{show,FLXRM-method}
\alias{show,FLXRMRglm-method}
\alias{show,Coefmat-method}
\alias{summary,FLXR-method}
\alias{summary,FLXRM-method}
\alias{summary,FLXRMRglm-method}
\alias{FLXRM-class}
\alias{FLXRP-class}
\title{Refit a Fitted Model}
\description{
  Refits an estimated flexmix model to obtain additional information
  like coefficient significance p-values for GLM regression.}
}
\usage{
\S4method{refit}{flexmix}(object, model=1, ...)
}
\arguments{
  \item{object}{an object of class \code{"flexmix"}}
  \item{model}{The model (for a multivariate response) that shall be
    refitted.}
  \item{\dots}{currently not used}
}
\details{
  The \code{refit} method for \code{FLXMRglm} models in
  combination with the \code{summary} method can be
  used to obtain the usual tests for significance of coefficients. Note
  that the tests are valid only if \code{flexmix} returned the maximum
  likelihood estimator of the parameters.
}
\keyword{methods}
\author{Friedrich Leisch}
\references{
  Friedrich Leisch. FlexMix: A general framework for finite mixture
  models and latent class regression in R. Journal of Statistical
  Software, 11(8), 2004. http://www.jstatsoft.org/v11/i08/
}
\examples{
data("NPreg")
ex1 <- flexmix(yn~x+I(x^2), data=NPreg, k=2)
ex1r <- refit(ex1)

## in one component all coefficients should be highly significant,
## in the other component only the linear term
summary(ex1r)
}

