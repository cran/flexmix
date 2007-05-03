%
%  Copyright (C) 2004-2005 Friedrich Leisch
%  $Id: FLXglm.Rd 3193 2007-01-04 10:35:17Z gruen $
%
\name{FLXMRglm}
\alias{FLXMRglm}
\alias{FLXglm}
\title{FlexMix Interface to Generalized Linear Models}
\description{
  This is the main driver for FlexMix interfacing the \code{\link{glm}}
  family of models.
}
\usage{
FLXMRglm(formula = . ~ .,
       family = c("gaussian", "binomial", "poisson", "Gamma"),
       offset = NULL)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{\link{flexmix}} using
    \code{\link{update.formula}}. Default is to use the original
    \code{\link{flexmix}} model
    formula.}
  \item{family}{A character string naming a \code{\link{glm}}
    family function.}
  \item{offset}{this can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.}
}
\details{
  See \code{\link{flexmix}} for examples.
}
\value{
  Returns an object of class \code{FLXMRglm}.
}
\author{Friedrich Leisch}
\references{
  Friedrich Leisch. FlexMix: A general framework for finite mixture
  models and latent class regression in R. Journal of Statistical
  Software, 11(8), 2004. http://www.jstatsoft.org/v11/i08/
}
\seealso{\code{\link{flexmix}}, \code{\link{glm}}}
\keyword{regression}
\keyword{models}
