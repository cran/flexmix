\name{refit}
\alias{refit}
\alias{refit-methods}
\alias{refit,flexmix-method}
\alias{refit,FLXglmmodel-method}
\alias{show,FLXrefit-method}
\alias{summary,FLXrefit-method}
\alias{summary,FLXrefitglm-method}
\title{Refit a Fitted Model}
\description{
  Refits an estimated flexmix model to obtain additional information
  like coefficient significance p-values for GLM regression.}
}
\usage{
\S4method{refit}{flexmix}(object, model=1, ...)
}
\arguments{
  \item{object}{an object of class \code{"flexmix"}}
  \item{model}{The model (for a multivariate response) that shall be
    refitted.}
  \item{\dots}{currently not used}
}
\details{
  Currently there is only a \code{refit} method for \code{FLXglm}
  models, which in combination with the \code{summary} method can be
  used to obtain the usual tests for significance of coefficients. Note
  that the tests are valid only if \code{flexmix} returned the maximum
  likelihood estimator of the parameters.
}
\keyword{methods}
\author{Friedrich Leisch}
\examples{
data(NPreg)
ex1 <- flexmix(yn~x+I(x^2), data=NPreg, k=2)
ex1r <- refit(ex1)

## in one component all coefficients should be highly significant,
## in the other component only the linear term
summary(ex1r)
}

