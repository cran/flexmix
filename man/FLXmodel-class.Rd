%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: FLXmodel-class.Rd 5229 2022-06-02 14:47:26Z gruen $
%
\name{FLXM-class}
\docType{class}
\alias{FLXM-class}
\alias{FLXMC-class}
\alias{FLXMR-class}
\alias{FLXMCsparse-class}
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
    \item{\code{fit}:}{Function returning an \code{FLXcomponent} object.}
    \item{\code{defineComponent}:}{Function or expression to determine the
      \code{FLXcomponent} object given the parameters.}
    \item{\code{weighted}:}{Logical indicating whether \code{fit} can do
    weighted likelihood maximization.}
    \item{\code{name}:}{Character string used in print methods.}
    \item{\code{formula}:}{Formula describing the model.}
    \item{\code{fullformula}:}{Resulting formula from updating the model
      formula with the formula specified in the call to \code{flexmix}.}
    \item{\code{x}:}{Model matrix.}
    \item{\code{y}:}{Model response.}
    \item{\code{terms}, \code{xlevels}, \code{contrasts}:}{Additional
      information for model matrix.}
    \item{\code{preproc.x}:}{Function for preprocessing matrix \code{x}
      before the EM algorithm starts, by default the identity function.}
    \item{\code{preproc.y}:}{Function for preprocessing matrix \code{y}
      before the EM algorithm starts, by default the identity function.}
  }
}
\details{
  The most generic class is the virtual class \code{FLXM}. The classes
  \code{FLXMC} for model-based clustering and \code{FLXMR} for
  clusterwise regression extend the virtual class. Both have further
  more specific model classes which inherit from them.

  Model class \code{FLXMCsparse} allows for model-based clustering with
  a sparse matrix as data input.
}
\author{Friedrich Leisch and Bettina Gruen}
\keyword{classes}
