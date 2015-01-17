%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: Lapply-methods.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{Lapply-methods}
\docType{methods}
\title{Methods for Function Lapply}
\alias{Lapply,FLXRmstep-method}
\description{Apply a function to each component of a finite mixture}
\usage{
\S4method{Lapply}{FLXRmstep}(object, FUN, model = 1, component = TRUE, ...)
}
\arguments{
  \item{object}{S4 class object.}
  \item{FUN}{The function to be applied.}
  \item{model}{The model (for a multivariate response) that shall be
    used.}
  \item{component}{Index vector for selecting the components.}
  \item{\dots}{Optional arguments to \code{FUN}.}
}
\section{Methods}{
  \describe{
    \item{object = FLXRmstep:}{Apply a function to each component of a
      refitted \code{flexmix} object using method = \code{"mstep"}.}
  }  
}
\details{
  \code{FUN} is found by a call to \code{match.fun} and typically is specified
  as a function or a symbol (e.g. a backquoted name) or a character
  string specifying a function to be searched for from the
  environment of the call to \code{Lapply}.
}
\value{
  A list of the length equal to the number of components specified is
  returned, each element of which is the result of applying \code{FUN} to the
  specified component of the refitted mixture model.
}
\keyword{methods}
\author{Friedrich Leisch and Bettina Gruen}
\examples{
data("NPreg", package = "flexmix")
ex2 <- flexmix(yn ~ x, data = NPreg, k = 2, model = list(FLXMRglm(yn ~ 
    . + I(x^2)), FLXMRglm(yp ~ ., family = "poisson")))
ex2r <- refit(ex2, method = "mstep")
Lapply(ex2r, "vcov", 2)
}
