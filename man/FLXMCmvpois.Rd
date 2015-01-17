%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: FLXMCmvpois.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{FLXMCmvpois}
\alias{FLXMCmvpois}
\title{FlexMix Poisson Clustering Driver}
\description{
  This is a model driver for \code{\link{flexmix}} implementing
  model-based clustering of Poisson distributed data. 
}
\usage{
FLXMCmvpois(formula = . ~ .)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{\link{flexmix}} using
    \code{\link{update.formula}}. Only the left-hand side (response) of
    the formula is used. Default is to use the original
    \code{\link{flexmix}} model
    formula.}
}
\details{
  This can be used to cluster Poisson distributed data where given the
  component membership the variables are mutually independent.
}
\value{
  \code{FLXMCmvpois} returns an object of class \code{FLXMC}.
}
\author{Friedrich Leisch and Bettina Gruen}
\seealso{\code{\link{flexmix}}}
\keyword{cluster}
