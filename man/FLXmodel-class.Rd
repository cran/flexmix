\name{FLXmodel-class}
\docType{class}
\alias{FLXmodel-class}
\alias{show,FLXmodel-method}
\title{Class "FLXmodel"}
\description{FlexMix model specification.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("FLXmodel", ...)},
  typically inside driver functions like \code{\link{FLXglm}} or
  \code{\link{FLXmclust}}. 
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Function returning an \code{FLXmodel} object.}
    \item{\code{weighted}:}{Logical indicating whether \code{fit} can do
    weighted likelihood maximization.}
    \item{\code{name}:}{Character string used in print methods.}
    \item{\code{formula}:}{Formula describing the model.}
    \item{\code{x}:}{Model matrix.}
    \item{\code{y}:}{Model response.}
  }
}

\author{Friedrich Leisch}
\keyword{classes}
