%
%  Copyright (C) 2004-2005 Friedrich Leisch
%  $Id: ExNPreg.Rd 3221 2007-01-22 14:27:40Z gruen $
%
\name{ExNPreg}
\alias{ExNPreg}
\alias{NPreg}
%- Also NEED an '\alias' for EACH other topic documented here.
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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of observations per latent class.}
}
\details{
  The \code{NPreg} data frame can be re-created by \code{ExNPreg(100)}
  using \code{set.seed(2602)}, it has been saved as a data set for
  simplicity of examples only.
}
\examples{
data("NPreg")
plot(yn~x, data=NPreg, col=class)
plot(yp~x, data=NPreg, col=class)
plot(yb~x, data=NPreg, col=class)
}
\keyword{datasets}
