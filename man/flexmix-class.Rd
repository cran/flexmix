\name{flexmix-class}
\docType{class}
\alias{flexmix-class}
\alias{cluster}
\alias{cluster,flexmix-method}
\alias{parameters}
\alias{parameters,flexmix-method}
\alias{posterior}
\alias{posterior,flexmix-method}
\title{Class "flexmix"}
\description{A fitted \code{\link{flexmix}} model.}
\usage{
\S4method{cluster}{flexmix}(object)
\S4method{parameters}{flexmix}(object, component=1, model=1)
\S4method{posterior}{flexmix}(object)
}
\arguments{
  \item{object}{An object of class "flexmix".}
  \item{component}{Number of component.}
  \item{model}{Number of model.}
}
\section{Slots}{
  \describe{
    \item{\code{model}:}{List of \code{FLXmodel} objects.}
    \item{\code{prior}:}{Numeric vector with prior probabilities of clusters.}
    \item{\code{posterior}:}{Named list with elements \code{scaled}
      and \code{unscaled}, both matrices with one row per observation
      and one column per cluster.}
    \item{\code{iter}:}{Number of EM
      iterations.}
    \item{\code{k}:}{Number of clusters.}
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
    \item{\code{converged}}{Logical, \code{TRUE} if EM algorithm converged.}
  }
}
\section{Accessor Functions}{
  The following functions should be used for accessing the corresponding
  slots:
  \describe{
    \item{\code{cluster}:}{Cluster assignments of observations.}
    \item{\code{parameters}:}{The parameters for each model and component,
      return value depends on the model.}
    \item{\code{posterior}:}{A matrix of posterior probabilities for each
      observation.}
  }
}
\author{Friedrich Leisch}
\keyword{classes}
