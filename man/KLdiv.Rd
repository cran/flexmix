%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: KLdiv.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{KLdiv}
\alias{KLdiv,matrix-method}
\alias{KLdiv,flexmix-method}
\alias{KLdiv,FLXMRglm-method}
\alias{KLdiv,FLXMC-method}
\title{Kullback-Leibler Divergence}
\description{
  Estimate the Kullback-Leibler divergence of several distributions.}
\usage{
\S4method{KLdiv}{matrix}(object, eps = 10^-4, overlap = TRUE,...)
\S4method{KLdiv}{flexmix}(object, method = c("continuous", "discrete"), ...)
}
\arguments{
  \item{object}{See Methods section below.}
  \item{method}{The method to be used; "continuous" determines
    the Kullback-Leibler divergence between the unweighted theoretical
    component distributions and the unweighted posterior probabilities
    at the observed points are used by "discrete".}
  \item{eps}{Probabilities below this threshold are replaced by this
    threshold for numerical stability.}
  \item{overlap}{Logical, do not determine the KL divergence for
    those pairs where for each point at least one of the densities has a
    value smaller than \code{eps}.}
  \item{...}{Passed to the matrix method.}
}
\section{Methods}{
  \describe{
    \item{object = "matrix":}{Takes as input a matrix of
      density values with one row per observation and one column per
      distribution.}
    \item{object = "flexmix":}{Returns the Kullback-Leibler divergence
      of the mixture components.}
}}
\details{
  Estimates \deqn{\int f(x) (\log f(x) - \log g(x)) dx}
  for distributions with densities \eqn{f()} and \eqn{g()}.
}
\value{
  A matrix of KL divergences where the rows correspond to using the
  respective distribution as \eqn{f()} in the formula above.
}
\note{
  The density functions are modified to have equal support.
  A weight of at least \code{eps} is given to each
  observation point for the modified densities.
}
\keyword{methods}
\author{Friedrich Leisch and Bettina Gruen}
\references{
  S. Kullback and R. A. Leibler. On information and sufficiency.\emph{The
  Annals of Mathematical Statistics}, \bold{22}(1), 79--86, 1951.
  
  Friedrich Leisch. Exploring the structure of mixture model
  components. In Jaromir Antoch, editor, Compstat 2004--Proceedings in
  Computational Statistics, 1405--1412. Physika Verlag, Heidelberg,
  Germany, 2004. ISBN 3-7908-1554-3.
}
\examples{
## Gaussian and Student t are much closer to each other than
## to the uniform:

x <- seq(-3, 3, length = 200)
y <- cbind(u = dunif(x), n = dnorm(x), t = dt(x, df = 10))

matplot(x, y, type = "l")
KLdiv(y)

if (require("mlbench")) {
set.seed(2606)
x <-  mlbench.smiley()$x
model1 <- flexmix(x ~ 1, k = 9, model = FLXmclust(diag = FALSE),
                  control  =  list(minprior = 0))
plotEll(model1, x)
KLdiv(model1)
}
}