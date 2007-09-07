%
%  Copyright (C) 2004-2005 Friedrich Leisch
%  $Id: refit.Rd 3595 2007-06-29 16:08:37Z gruen $
%
\name{refit}
\alias{refit,flexmix,ANY-method}
\alias{refit,flexmix,missing-method}
\alias{refit,flexmix,listOrdata.frame-method}
\alias{refit,FLXM,missing-method}
\alias{refit,FLXMRglm,missing-method}
\alias{refit,FLXMRglm,listOrdata.frame-method}
\alias{refit,FLXMRglmfix,missing-method}
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
\S4method{refit}{flexmix,ANY}(object, newdata, model=1, which = c("model",
"concomitant"), summary=TRUE, ...)
}
\arguments{
  \item{object}{An object of class \code{"flexmix"}}
  \item{newdata}{Optional new data.}
  \item{model}{The model (for a multivariate response) that shall be
    refitted.}
  \item{which}{Specifies if a component specific model or the
    concomitant variable model is refitted.}
  \item{summary}{A logical if the summary output should also be
    calculated.}
  \item{\dots}{Currently not used}
}
\details{
  The \code{refit} method for \code{FLXMRglm} models in
  combination with the \code{summary} method can be
  used to obtain the usual tests for significance of coefficients. Note
  that the tests are valid only if \code{flexmix} returned the maximum
  likelihood estimator of the parameters. For this method the returned
  object contains a \code{glm} object for each component where the
  elements \code{model} which is the model frame and \code{data} which
  contains the original dataset are missing.
}
\keyword{methods}
\author{Friedrich Leisch}
\references{
  Friedrich Leisch. FlexMix: A general framework for finite mixture
  models and latent class regression in R. Journal of Statistical
  Software, 11(8), 2004. http://www.jstatsoft.org/v11/i08/
}
\section{Warning}{
  The standard deviations are determined separately for each of the
  components using the a-posteriori probabilities as weights without
  accounting for the fact that the components have been simultaneously
  estimated. The derived standard deviations are hence approximative and
  should only be used in an exploratory way, as they are underestimating
  the uncertainty given that the missing information of the component
  memberships are replaced by the expected values.

  The \code{newdata} argument can only be specified for refitting
  \code{FLXMRglm} components. A variant of \code{glm} for weighted ML
  estimation is used for fitting the components and full \code{glm}
  objects are returned. Please note that in this case the data and the
  model frame are stored for each component which can significantly
  increase the object size.
}
\examples{
data("NPreg")
ex1 <- flexmix(yn~x+I(x^2), data=NPreg, k=2)
ex1r <- refit(ex1)

## in one component all coefficients should be highly significant,
## in the other component only the linear term
summary(ex1r)
}

