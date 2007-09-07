\name{Lapply-methods}
\docType{methods}
\title{Methods for Function Lapply}
\alias{Lapply,FLXR-method}
\alias{Lapply,FLXRM-method}
\description{Apply a function to each component of a finite mixture}
\usage{
\S4method{Lapply}{FLXR}(object, FUN, component = TRUE, ...)
}
\arguments{
  \item{object}{S4 class object.}
  \item{FUN}{The function to be applied.}
  \item{component}{Index vector for selecting the components.}
  \item{\dots}{Optional arguments to 'FUN'.}
}
\section{Methods}{
  \describe{
    \item{object = FLXR:}{Apply a function to each component of a
      refitted \code{flexmix} object}
  }  
}
\details{
  'FUN' is found by a call to 'match.fun' and typically is specified
  as a function or a symbol (e.g. a backquoted name) or a character
  string specifying a function to be searched for from the
  environment of the call to 'Lapply'.
}
\value{
  A list of the length equal to the number of components specified is
  returned, each element of which is the result of applying 'FUN' to the
  specified component of the refitted mixture model.
}
\keyword{methods}
\author{Bettina Gruen}
\examples{
data("NPreg")
ex2 <- flexmix(yn ~ x, data = NPreg, k = 2, model = list(FLXMRglm(yn ~ 
    . + I(x^2)), FLXMRglm(yp ~ ., family = "poisson")))
ex2r <- refit(ex2)
Lapply(ex2r, "vcov", 2)
}
