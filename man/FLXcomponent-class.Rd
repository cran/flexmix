\name{FLXcomponent-class}
\docType{class}
\alias{FLXcomponent-class}
\title{Class "FLXcomponent"}
\description{A fitted component of a \code{\link{flexmix}} model.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("FLXcomponent", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{df}:}{Number of paramters used by the component.}
    \item{\code{logLik}:}{Function computing the log-likelihood of
      observations.}
    \item{\code{parameters}:}{List with model parameters.}
    \item{\code{predict}:}{Function predicting response for new data.}
  }
}
\author{Friedrich Leisch}
\keyword{classes}
