%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: fitted.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{fitted-methods}
\docType{methods}
\alias{fitted,flexmix-method}
\alias{fitted,FLXM-method}
\alias{fitted,FLXR-method}
\alias{fitted,FLXRMRglm-method}
\title{Extract Model Fitted Values}
\description{
  Extract fitted values for each component from a flexmix object.
}
\usage{
\S4method{fitted}{flexmix}(object, drop = TRUE, aggregate = FALSE, ...)
}
\arguments{
  \item{object}{an object of class \code{"flexmix"} or
    \code{"FLXR"}}
  \item{drop}{logical, if \code{TRUE} then the function tries to
    simplify the return object by combining lists of length 1 into
    matrices.}
  \item{aggregate}{logical, if \code{TRUE} then the fitted values for
    each model aggregated over the components are returned.}
  \item{\dots}{currently not used}
}
\keyword{methods}
\author{Friedrich Leisch and Bettina Gruen}
\examples{
data("NPreg", package = "flexmix")
ex1 <- flexmix(yn ~ x + I(x^2), data = NPreg, k = 2)
matplot(NPreg$x, fitted(ex1), pch = 16, type = "p")
points(NPreg$x, NPreg$yn)
}

