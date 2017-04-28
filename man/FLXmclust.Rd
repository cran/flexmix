%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: FLXmclust.Rd 5115 2017-04-07 08:18:13Z gruen $
%
\name{FLXMCmvnorm}
\alias{FLXMCmvnorm}
\alias{FLXMCnorm1}
\alias{FLXmclust}
\title{FlexMix Clustering Demo Driver}
\description{
  These are demo drivers for \code{\link{flexmix}} implementing
  model-based clustering of Gaussian data. 
}
\usage{
FLXMCmvnorm(formula = . ~ ., diagonal = TRUE)
FLXMCnorm1(formula = . ~ .)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{\link{flexmix}} using
    \code{\link{update.formula}}. Only the left-hand side (response) of
    the formula is used. Default is to use the original
    \code{\link{flexmix}} model
    formula.}
  \item{diagonal}{If \code{TRUE}, then the covariance matrix of the
    components is restricted to diagonal matrices.}
}
\details{
  This is mostly meant as a demo for FlexMix driver programming, you
  should also look at package \pkg{mclust} for real
  applications. \code{FLXMCmvnorm} clusters multivariate data,
  \code{FLXMCnorm1} univariate data. In the latter case smart
  initialization is important, see the example below. 
}
\value{
  \code{FLXMCmvnorm} returns an object of class \code{FLXMC}.
}
\author{Friedrich Leisch and Bettina Gruen}
\references{
  Friedrich Leisch. FlexMix: A general framework for finite mixture
  models and latent class regression in R. \emph{Journal of Statistical
  Software}, \bold{11}(8), 2004. doi:10.18637/jss.v011.i08
}
\seealso{\code{\link{flexmix}}}
\keyword{cluster}
\examples{
data("Nclus", package = "flexmix")

require("MASS")
eqscplot(Nclus)

## This model is wrong (one component has a non-diagonal cov matrix)
ex1 <- flexmix(Nclus ~ 1, k = 4, model = FLXMCmvnorm())
print(ex1)
plotEll(ex1, Nclus)

## True model, wrong number of components
ex2 <- flexmix(Nclus ~ 1, k = 6, model = FLXMCmvnorm(diagonal = FALSE))  
print(ex2)

plotEll(ex2, Nclus)

## Get parameters of first component
parameters(ex2, component = 1)

## Have a look at the posterior probabilies of 10 random observations
ok <- sample(1:nrow(Nclus), 10)
p  <- posterior(ex2)[ok, ]
p

## The following two should be the same
max.col(p)
clusters(ex2)[ok]
\testonly{
  stopifnot(all.equal(max.col(p), clusters(ex2)[ok]))
}
## Now try the univariate case
plot(density(Nclus[, 1]))

ex3 <- flexmix(Nclus[, 1] ~ 1, cluster = cut(Nclus[, 1], 3),
               model = FLXMCnorm1())
ex3
parameters(ex3)
}
