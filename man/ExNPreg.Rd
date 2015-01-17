%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: ExNPreg.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{ExNPreg}
\alias{ExNPreg}
\alias{NPreg}
\title{Artificial Example for Normal, Poisson and Binomial Regression}
\description{
  A simple artificial regression example with 2 latent classes, one
  independent variable
  (uniform on \eqn{[0,10]}), and three dependent variables with
  Gaussian, Poisson and Binomial distribution, respectively.
}
\usage{
ExNPreg(n)
data("NPreg")
}
\arguments{
  \item{n}{Number of observations per latent class.}
}
\details{
  The \code{NPreg} data frame can be re-created by \code{ExNPreg(100)}
  using \code{set.seed(2602)}, it has been saved as a data set for
  simplicity of examples only.
}
\examples{
data("NPreg", package = "flexmix")
plot(yn ~ x, data = NPreg, col = class)
plot(yp ~ x, data = NPreg, col = class)
plot(yb ~ x, data = NPreg, col = class)
}
\keyword{datasets}
