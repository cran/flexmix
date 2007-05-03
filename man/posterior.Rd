%
%  Copyright (C) 2007 Bettina Gruen and Friedrich Leisch
%  $Id: flexmix-class.Rd 3122 2006-11-05 15:10:25Z leisch $
%
\name{posterior}
\alias{cluster,flexmix,missing-method}
\alias{cluster,flexmix,ANY-method}
\alias{posterior,flexmix,missing-method}
\alias{posterior,FLXdist,listOrdata.frame-method}
\title{Determine cluster membership and posterior probabilities}
\description{Determine posterior probabilities for a fitted
  \code{flexmix}  or unfitted \code{FLXdist} model.} 
\usage{
\S4method{posterior}{flexmix,missing}(object, newdata, ...)
\S4method{posterior}{FLXdist,listOrdata.frame}(object, newdata, ...)
\S4method{cluster}{flexmix,missing}(object)
\S4method{cluster}{flexmix,ANY}(object, newdata, ...)
}
\arguments{
  \item{object}{An object of class "flexmix" or "FLXdist".}
  \item{newdata}{Dataframe or list containing new data. If missing the
    posteriors of the original observations are returned.}
  \item{\dots}{Currently not used.}
}
\author{Bettina Gruen and Friedrich Leisch}
\keyword{methods}


