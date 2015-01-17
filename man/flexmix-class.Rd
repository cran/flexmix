%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: flexmix-class.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{flexmix-class}
\docType{class}
\alias{flexmix-class}
\title{Class "flexmix"}
\description{A fitted \code{\link{flexmix}} model.}
\section{Slots}{
  \describe{
    \item{\code{model}:}{List of \code{FLXM} objects.}
    \item{\code{prior}:}{Numeric vector with prior probabilities of clusters.}
    \item{\code{posterior}:}{Named list with elements \code{scaled}
      and \code{unscaled}, both matrices with one row per observation
      and one column per cluster.}
    \item{\code{iter}:}{Number of EM iterations.}
    \item{\code{k}:}{Number of clusters after EM.}
    \item{\code{k0}:}{Number of clusters at start of EM.}
    \item{\code{cluster}:}{Cluster assignments of observations.}
    \item{\code{size}:}{Cluster sizes.}
    \item{\code{logLik}:}{Log-likelihood at EM convergence.}
    \item{\code{df}:}{Total number of
      parameters of the model.}
    \item{\code{components}:}{List describing
      the fitted components using \code{FLXcomponent} objects.}
    \item{\code{formula}:}{Object of class \code{"formula"}.}
    \item{\code{control}:}{Object of class \code{"FLXcontrol"}.}
    \item{\code{call}:}{The function call
      used to create the object.}
    \item{\code{group}:}{Object of class \code{"factor"}.}
    \item{\code{converged}:}{Logical, \code{TRUE} if EM algorithm converged.}
    \item{\code{concomitant}:}{Object of class \code{"FLXP"}..}
    \item{\code{weights}:}{Optional weights of the observations.}
  }
}
\section{Extends}{
  Class \code{FLXdist}, directly.
}

\section{Accessor Functions}{
  The following functions should be used for accessing the corresponding
  slots:
  \describe{
    \item{\code{cluster}:}{Cluster assignments of observations.}
    \item{\code{posterior}:}{A matrix of posterior probabilities for each
      observation.}
  }
}
\author{Friedrich Leisch and Bettina Gruen}
\keyword{classes}
