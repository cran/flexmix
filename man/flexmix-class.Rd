\name{flexmix-class}
\docType{class}
\alias{flexmix-class}
\title{Class "flexmix"}
\description{A fitted \code{\link{flexmix}} model.}
\section{Objects from the Class}{
Objects are typically created by calling function \code{\link{flexmix}}.
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
    \item{\code{cluster}:}{Cluster assignments of observations.}
    \item{\code{logLik}:}{Log-likelihood at EM convergence.}
    \item{\code{df}:}{Total number of
      parameters of the model.}
    \item{\code{components}:}{List describing
      the fitted components using \code{FLXcomponent} objects.}
    \item{\code{formula}:}{Object of class \code{"formula"}.}
    \item{\code{control}:}{Object of class \code{"FLXcontrol"}.}
    \item{\code{call}:}{The function call
      used to create the object.}
  }
}

\section{Methods}{
  \describe{
    \item{AIC}{\code{signature(object = "flexmix")}: calculate Akaike's
      information criterion.}
    \item{BIC}{\code{signature(object = "flexmix")}: calculate the
      Bayesian information criterion.}
    \item{logLik}{\code{signature(object = "flexmix")}: calculate the
      log-likelihood for the model.}
    \item{parameters}{\code{signature(object = "flexmix")}: extract
      parameters of model components.}
    \item{predict}{\code{signature(object = "flexmix")}: predict
      response for new data.}
  }
}
\author{Friedrich Leisch}
\keyword{classes}
