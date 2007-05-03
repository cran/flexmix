%
%  Copyright (C) 2004-2005 Friedrich Leisch
%  $Id: FLXmodel-class.Rd 3181 2006-12-18 10:43:28Z gruen $
%
\name{FLXM-class}
\docType{class}
\alias{FLXM-class}
\alias{FLXMC-class}
\alias{FLXMR-class}
\alias{show,FLXM-method}
\title{Class "FLXM"}
\description{FlexMix model specification.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("FLXM", ...)},
  typically inside driver functions like \code{\link{FLXMRglm}} or
  \code{\link{FLXMCmvnorm}}. 
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Function returning an \code{FLXM} object.}
    \item{\code{weighted}:}{Logical indicating whether \code{fit} can do
    weighted likelihood maximization.}
    \item{\code{name}:}{Character string used in print methods.}
    \item{\code{formula}:}{Formula describing the model.}
    \item{\code{x}:}{Model matrix.}
    \item{\code{y}:}{Model response.}
    \item{\code{preproc.x}:}{Function for preprocessing matrix \code{x}
      before the EM algorithm starts, by default the identity function.}
    \item{\code{preproc.y}:}{Function for preprocessing matrix \code{y}
      before the EM algorithm starts, by default the identity function.}
  }
}

\author{Friedrich Leisch}
\keyword{classes}
