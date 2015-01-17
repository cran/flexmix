%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: FLXP-class.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{FLXP-class}
\docType{class}
\alias{FLXP-class}
\title{Class "FLXP"}
\alias{initialize,FLXP-method}
\alias{FLXPconstant-class}
\alias{FLXPmultinom-class}
\description{
  Concomitant model class.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("FLXP", ...)},
  typically inside driver functions like \code{\link{FLXPconstant}} or
  \code{\link{FLXPmultinom}}. 
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Character string used in print methods.}
    \item{\code{formula}:}{Formula describing the model.}
    \item{\code{x}:}{Model matrix.}
    \item{\code{fit}:}{Function returning the fitted prior probabilities.}
    \item{\code{refit}:}{Function returning the fitted concomitant model.}
    \item{\code{coef}:}{Matrix containing the fitted parameters.}
    \item{\code{df}:}{Function for determining the number of degrees of
      freedom used.}
  }
}
\author{Friedrich Leisch and Bettina Gruen}
\keyword{classes}

