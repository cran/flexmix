\name{ExNclus}
\alias{ExNclus}
\alias{Nclus}
\title{Artificial Example for Normal Clustering}
\description{
  A simple artificial regression example with 4 latent classes, all of
  them having a Gaussian distribution. See the function definition for
  true means and covariances.
}
\usage{
ExNclus(n)
data(Nclus)
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
data(Nclus)
require("MASS")
eqscplot(Nclus, col=rep(1:4, c(100, 100, 150, 200)))
}
\keyword{datasets}
