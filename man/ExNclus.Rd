%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: ExNclus.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{ExNclus}
\alias{ExNclus}
\alias{Nclus}
\title{Artificial Example with 4 Gaussians}
\description{
  A simple artificial example for normal clustering with 4 latent
  classes, all of
  them having a Gaussian distribution. See the function definition for
  true means and covariances.
}
\usage{
ExNclus(n)
data("Nclus")
}
\arguments{
  \item{n}{Number of observations in the two small latent classes.}
}
\details{
  The \code{Nclus} data set can be re-created by \code{ExNclus(100)}
  using \code{set.seed(2602)}, it has been saved as a data set for
  simplicity of examples only.
}
\examples{
data("Nclus", package = "flexmix")
require("MASS")
eqscplot(Nclus, col = rep(1:4, c(100, 100, 150, 200)))
}
\keyword{datasets}
