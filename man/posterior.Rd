\name{posterior}
\alias{posterior}
\alias{posterior,flexmix,missing-method}
\alias{posterior,FLXdist,data.frame-method}
\title{Determine posterior probabilities}
\description{Determine posterior probabilities for a fitted
  \code{flexmix}  or unfitted \code{FLXdist} model.} 
\usage{
\S4method{posterior}{flexmix,missing}(object, newdata, ...)
\S4method{posterior}{FLXdist,data.frame}(object, newdata, ...)
}
\arguments{
  \item{object}{An object of class "flexmix" or "FLXdist".}
  \item{newdata}{Dataframe containing new data. If missing the
    posteriors of the original observations are returned.}
  \item{\dots}{Currently not used.}
}
\author{Bettina Gruen}
\keyword{methods}


