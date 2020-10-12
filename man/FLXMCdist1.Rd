%
%  Copyright (C) 2016 Bettina Gruen
%  $Id: FLXMCdist1.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{FLXMCdist1}
\alias{FLXMCdist1}
\title{FlexMix Clustering of Univariate Distributions}
\description{
  These are drivers for \code{\link{flexmix}} implementing model-based
  clustering of univariate data using different distributions for the
  component-specific models.
}
\usage{
FLXMCdist1(formula = . ~ ., dist, ...)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{\link{flexmix}} using
    \code{\link{update.formula}}. Only the left-hand side (response) of
    the formula is used. Default is to use the original
    \code{\link{flexmix}} model
    formula.}
  \item{dist}{Character string indicating the component-specific
  univariate distribution.}
  \item{...}{Arguments for the specific model drivers.}
}
\details{
    Currently drivers for the following distributions are available:
    \enumerate{
        \item Lognormal (\code{"lnorm"})
        \item inverse Gaussian (\code{"invGauss"} using \code{\link[SuppDists:invGauss]{dinvGauss}})
        \item gamma (\code{"gamma"})
        \item exponential (\code{"exp"})
        \item Weibull (\code{"weibull"})
        \item Burr (\code{"burr"} using \code{\link[actuar:Burr]{dburr}})
        \item Inverse Burr (\code{"invburr"} using \code{\link[actuar:Burr]{dinvburr}})
    }
}
\value{
  \code{FLXMCdist1} returns an object of class \code{FLXMC}.
}
\author{Friedrich Leisch and Bettina Gruen}
\references{
  Tatjana Miljkovic and Bettina Gruen. Modeling loss data using
  mixtures of distributions. \emph{Insurance: Mathematics and
    Economics}, \bold{70}, 387-396, 2016.
  doi:10.1016/j.insmatheco.2016.06.019
}
\seealso{\code{\link{flexmix}}}
\keyword{cluster}
\examples{
if (require("actuar")) {
    set.seed(123)
    y <- c(rexp(100, 10), rexp(100, 1))
    ex <- flexmix(y ~ 1, cluster = rep(1:2, each = 100), model = FLXMCdist1(dist = "exp"))    
    parameters(ex)
}
}
