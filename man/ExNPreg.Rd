%
%  Copyright (C) 2004-2011 Friedrich Leisch and Bettina Gruen
%  $Id: ExNPreg.Rd 4666 2011-02-23 15:52:35Z gruen $
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
data("NPreg")
plot(yn~x, data=NPreg, col=class)
plot(yp~x, data=NPreg, col=class)
plot(yb~x, data=NPreg, col=class)
}
\keyword{datasets}
