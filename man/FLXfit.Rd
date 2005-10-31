%
%  Copyright (C) 2004-2005 Friedrich Leisch
%  $Id: FLXfit.Rd 1673 2005-06-16 16:15:42Z gruen $
%
\name{FLXfit}
\alias{FLXfit}
\alias{FLXfit,list-method}
\title{Fitter Function for FlexMix Models}
\description{
  This is the basic computing engine called by \code{\link{flexmix}}, it
  should usually not be used directly.
}
\usage{
FLXfit(model, control, postunscaled = NULL, group)
}
\arguments{
  \item{postunscaled}{Initial a posteriori probabilities of the
    observations at the start of the EM algorithm.}
  \item{model}{List of \code{FLXmodel} objects.}
  \item{control}{Object of class \code{FLXcontrol}.}
  \item{group}{Factor with optional grouping of observations.}
}
\value{
  Returns an object of class \code{flexmix}.
}
\author{Friedrich Leisch}
\seealso{\code{\link{flexmix}}, \code{\link{flexmix-class}}}
\keyword{regression}
\keyword{cluster}
