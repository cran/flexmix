\name{FLXnested-class}
\docType{class}
\alias{FLXnested-class}
\alias{coerce,list,FLXnested-method}
\alias{coerce,NULL,FLXnested-method}
\alias{coerce,numeric,FLXnested-method}
\alias{initialize,FLXnested-method}

\title{Class "FLXnested"}
\description{Specification of nesting structure for regression coefficients.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("FLXnested",
  formula, k, ...)}. In addition, named lists can be coerced to
\code{FLXnested} objects, names are completed if unique.
}
\section{Slots}{
  \describe{
    \item{\code{formula}:}{Object of class \code{"list"} containing the
      formula for determining the model matrix for each nested effect.}
    \item{\code{k}:}{Object of class \code{"numeric"} specifying the
      number of components in each group.}
  }
}
\author{Bettina Gruen}
\keyword{classes}
