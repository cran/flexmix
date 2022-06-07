%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: FLXglm.Rd 5229 2022-06-02 14:47:26Z gruen $
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
  \item{offset}{This can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.}
}
\details{
  See \code{\link{flexmix}} for examples.
}
\value{
  Returns an object of class \code{FLXMRglm} inheriting from \code{FLXMR}.
}
\author{Friedrich Leisch and Bettina Gruen}
\references{
  Friedrich Leisch. FlexMix: A general framework for finite mixture
  models and latent class regression in R. \emph{Journal of Statistical
  Software}, \bold{11}(8), 2004. doi:10.18637/jss.v011.i08
}
\seealso{\code{\link{flexmix}}, \code{\link{glm}}}
\keyword{regression}
\keyword{models}
