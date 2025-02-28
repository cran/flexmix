\name{flxglht}
\docType{methods}
\alias{flxglht}
\alias{flxglht,flexmix,character-method}
\alias{flxglht,FLXRoptim,character-method}
\title{
  Simultaneous Inference in Finite Mixtures of Regression Models
}
\description{
  Extracts coefficient estimates and their covariance estimate to
  perform tests for zero coefficients or constant effects re-using
  functionality from package \pkg{multcomp}.
}
\usage{
flxglht(model, linfct, ...)
\S4method{flxglht}{flexmix,character}(model, linfct, ...)
\S4method{flxglht}{FLXRoptim,character}(model, linfct, ...)
}
\arguments{

  \item{model}{
    Either a fitted mixture model of class \code{"flexmix"}
    or a re-fitted mixture model of class \code{"FLXRoptim"}.
  }
  \item{linfct}{
    A character string which can be either \code{"zero"}
    or \code{"tukey"}.
  }
  \item{\dots}{
    Currently not used.
  }
}
\details{
  Only tested for finite mixture models fitted with driver
  \code{\link{FLXMRglm}}.
}
\value{
  An object of class \code{"glht"}.
}
\references{
  Friedrich Leisch and Torsten Hothorn. Simultaneous Inference in Finite
  Mixtures of Regression Models. \emph{Austrian Journal of Statistics},
  forthcoming.
}
\author{
  Friedrich Leisch
}

\seealso{
  \code{\link[multcomp]{glht}}
}
\examples{
data("NPreg", package = "flexmix")

ex1 <- flexmix(yn ~ x + I(x^2), data = NPreg, k = 2,
               control = list(verb = 5, iter = 100))

zero_effect <- flxglht(ex1, "zero")
zero_effect
summary(zero_effect)

comp_effect <- flxglht(ex1, "tukey")
comp_effect
summary(comp_effect)
}
\keyword{cluster}
\keyword{regression}
\keyword{methods}
