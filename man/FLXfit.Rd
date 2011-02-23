%
%  Copyright (C) 2004-2011 Friedrich Leisch and Bettina Gruen
%  $Id: FLXfit.Rd 4666 2011-02-23 15:52:35Z gruen $
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
FLXfit(model, concomitant, control, postunscaled = NULL, groups,
       weights)
}
\arguments{
  \item{model}{List of \code{FLXM} objects.}
  \item{concomitant}{Object of class \code{FLXP}.}
  \item{control}{Object of class \code{FLXcontrol}.}
    \item{weights}{A numeric vector of weights to be used in the fitting
      process.}
  \item{postunscaled}{Initial a-posteriori probabilities of the
    observations at the start of the EM algorithm.}
  \item{groups}{List with components \code{group} which is a factor with
    optional grouping of observations and \code{groupfirst} which is a
    logical vector for the first observation of each group.}
}
\value{
  Returns an object of class \code{flexmix}.
}
\author{Friedrich Leisch and Bettina Gruen}
\seealso{\code{\link{flexmix}}, \code{\link{flexmix-class}}}
\keyword{regression}
\keyword{cluster}
