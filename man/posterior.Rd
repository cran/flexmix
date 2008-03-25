%
%  Copyright (C) 2004-2008 Friedrich Leisch and Bettina Gruen
%  $Id: posterior.Rd 3934 2008-03-26 08:36:26Z gruen $
%
\name{posterior}
\alias{cluster,flexmix,missing-method}
\alias{cluster,FLXdist,ANY-method}
\alias{posterior,flexmix,missing-method}
\alias{posterior,FLXdist,listOrdata.frame-method}
\title{Determine cluster membership and posterior probabilities}
\description{Determine posterior probabilities or cluster memberships for a fitted
  \code{flexmix}  or unfitted \code{FLXdist} model.} 
\usage{
\S4method{posterior}{flexmix,missing}(object, newdata, unscaled = FALSE, ...)
\S4method{posterior}{FLXdist,listOrdata.frame}(object, newdata, unscaled = FALSE, ...)
\S4method{cluster}{flexmix,missing}(object, newdata, ...)
\S4method{cluster}{FLXdist,ANY}(object, newdata, ...)
}
\arguments{
  \item{object}{An object of class "flexmix" or "FLXdist".}
  \item{newdata}{Data frame or list containing new data. If missing the
    posteriors of the original observations are returned.}
  \item{unscaled}{Logical, if \code{TRUE} the component-specific
    likelihoods are returned.}
  \item{\dots}{Currently not used.}
}
\author{Friedrich Leisch and Bettina Gruen}
\keyword{methods}


