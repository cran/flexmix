\name{KLdiv}
\alias{KLdiv}
\alias{KLdiv-methods}
\alias{KLdiv,matrix-method}
\alias{KLdiv,flexmix-method}
\title{Kullbach Leibler Divergence}
\description{
  Estimate the Kullbach Leibler divergence of several distributions.}
\usage{
KLdiv(object, ...)
\S4method{KLdiv}{matrix}(object, eps=1e-4, ...)
}
\arguments{
  \item{object}{see Methods section below}
  \item{eps}{probabilities below this treshold are discarded for
    numerical stability}
  \item{...}{Passed to the matrix method.}
}
\section{Methods}{
  \describe{
    \item{object = "matrix":}{Takes as input a matrix of
      density values with one row per observation and one column per
      distribution.}
    \item{object = "flexmix":}{Returns the Kullbach Leibler divergence
      of the mixture components.}
}}
\details{
  Estimates \deqn{\int f(x) (\log f(x) - \log g(x)) dx}
  for distributions with densities \eqn{f()} and \eqn{g()}.
}
\value{
  A matrix of of KL divergences where the rows correspond to using the
  respective distribution as \eqn{f()} in the formula above.
}
\keyword{methods}
\author{Friedrich Leisch}
\references{
  Friedrich Leisch. Exploring the structure of mixture model
  components. In Jaromir Antoch, editor, Compstat 2004 - Proceedings in
  Computational Statistics, pages 1405-1412. Physika Verlag, Heidelberg,
  Germany, 2004. ISBN 3-7908-1554-3.
}
\examples{
x = (1:100)/100
## Gaussian and Student t are much closer to each other than
## to the uniform:
KLdiv(cbind(u=dunif(x), n=dnorm(x), t=dt(x, df=10)))
}
