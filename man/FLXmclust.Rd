\name{FLXmclust}
\alias{FLXmclust}
\alias{plotEll}
\title{FlexMix Clustering Demo Driver}
\description{
  This is a demo driver for \code{\link{flexmix}} implementing
  model-based clustering of Gaussian data. 
}
\usage{
FLXmclust(formula = . ~ ., diagonal = TRUE)
plotEll(object, data, ...)
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
  \item{object}{An object of class \code{flexmix} using \code{FLXmclust}
    model.}
  \item{data}{The data that were clustered.}
  \item{\dots}{Passed to \code{\link[MASS]{eqscplot}}.}
}
\details{
  This is meant as a demo for FlexMix driver programming, use package
  \pkg{mclust} for real applications. 
}
\value{
  \code{FLXmclust} returns an object of class \code{FLXmodel}.
}
\author{Friedrich Leisch}
\seealso{\code{\link{flexmix}}}
\keyword{cluster}
\examples{
data(Nclus)

require("MASS")
eqscplot(Nclus)

## This model is wrong (one component has a non-diagonal cov matrix)
ex1 <- flexmix(Nclus~1, k=4, model=FLXmclust())
print(ex1)
plotEll(ex1, Nclus)

## True model, wrong number of components
ex2 <- flexmix(Nclus~1, k=6, model=FLXmclust(diag=FALSE))  
print(ex2)

plotEll(ex2, Nclus)

parameters(ex2, component=1)
}
