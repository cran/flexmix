\name{fitted}
\docType{methods}
\alias{fitted-methods}
\alias{fitted,flexmix-method}
\alias{fitted,FLXrefit-method}
\alias{fitted,FLXrefitglm-method}
\title{Extract Model Fitted Values}
\description{
  Extract fitted values for each component from a flexmix object.
}
\usage{
\S4method{fitted}{flexmix}(object, drop=TRUE, ...)
}
\arguments{
  \item{object}{an object of class \code{"flexmix"} or
    \code{"FLXrefit"}}
  \item{drop}{logical, if \code{TRUE} then the function tries to
    simplify the return object by combining lists of length 1 into
    matrices.}
  \item{\dots}{currently not used}
}
\keyword{methods}
\author{Friedrich Leisch}
\examples{
data(NPreg)
ex1 <- flexmix(yn~x+I(x^2), data=NPreg, k=2)
matplot(NPreg$x, fitted(ex1), pch=16, type="p")
points(NPreg$x, NPreg$yn)
}

