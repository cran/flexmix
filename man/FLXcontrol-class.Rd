%
%  Copyright (C) 2004-2005 Friedrich Leisch
%  $Id: FLXcontrol-class.Rd 3692 2007-09-07 08:56:43Z gruen $
%
\name{FLXcontrol-class}
\docType{class}
\alias{FLXcontrol-class}
\alias{coerce,list,FLXcontrol-method}
\alias{coerce,NULL,FLXcontrol-method}
\title{Class "FLXcontrol"}
\description{Hyperparameters for the EM algorithm.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("FLXcontrol",
    ...)}. In addition, named lists can be coerced to \code{FLXcontrol}
  objects, names are completed if unique (see examples).
}
\section{Slots}{
  \describe{
    \item{\code{iter.max}:}{Maximum number of iterations.}
    \item{\code{minprior}:}{Minimum prior probability of clusters,
      components falling below this threshold are removed during the
      iteration.}
    \item{\code{tolerance}:}{The EM algorithm is stopped when the
      (relative) change of log-likelood is smaller than \code{tolerance}.}
    \item{\code{verbose}:}{If a positive integer, then the
      log-likelihood is reported every \code{verbose} iterations. If 0,
      no output is generated during model fitting.}
    \item{\code{classify}:}{Character string, one of \code{"auto"},
      \code{"weighted"}, \code{"hard"}, \code{"random"}.}
    \item{\code{nrep}:}{Reports the number of random initializations
      used to determine the mixture with \code{\link{stepFlexmix}()}.}
  }
  Run \code{new("FLXcontrol")} to see the default settings of all slots.
}
\author{Friedrich Leisch}
\keyword{classes}
\examples{
## have a look at the defaults
new("FLXcontrol")

## corce a list
mycont = list(iter=200, tol=0.001, class="r")
as(mycont, "FLXcontrol")
}
