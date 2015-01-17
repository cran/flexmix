%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: ExNclus.Rd 3912 2008-03-13 15:10:24Z gruen $
%
\name{ExLinear}
\alias{ExLinear}
\title{Artificial Data from a Generalized Linear Regression Mixture}
\description{
Generate random data mixed from k generalized linear regressions (GLMs).
}
\usage{
ExLinear(beta, n, xdist = "runif", xdist.args = NULL,
         family = c("gaussian","poisson"), sd = 1, ...)
}
\arguments{
  \item{beta}{A matrix of regression coefficients. Each row corresponds
    to a variable, each column to a mixture component. The first row is
    used as intercept.}
  \item{n}{Integer, the number of observations per component.}
  \item{xdist}{Character, name of a random number function for the
    explanatory variables.}
  \item{xdist.args}{List, arguments for the random number functions.}
  \item{family}{A character string naming a GLM
    family.Only \code{"gaussian"} and \code{"poisson"} are implemented
    at the moment.}
  \item{sd}{Numeric, the error standard deviation for each component for
    Gaussian responses.}
  \item{\dots}{Used as default for \code{xdist.args} if that is not
    specified.}
}
\details{
  First, arguments \code{n} (and \code{sd} for Gaussian response) are
  recycled to the number of mixture components \code{ncol(beta)}, and
  arguments \code{xdist} and \code{xdist.args} are recycled to the
  number of explanatory variables \code{nrow(beta)-1}. Then a design
  matrix is created for each mixture component by drawing random numbers
  from \code{xdist}. For each component, the design matrix is multiplied
  by the regression coefficients to form the linear predictor. For
  Gaussian responses the identity link is used, for Poisson responses
  the log link.

  The true cluster memberships are returned as attribute \code{"clusters"}.
}
\examples{
## simple example in 2d
beta <- matrix(c(1, 2, 3, -1), ncol = 2)
sigma <- c(.5, 1)
df1 <- ExLinear(beta, 100, sd = sigma, min = -1, max = 2)
plot(y~x1, df1, col = attr(df1, "clusters"))

## add a second explanatory variable with exponential distribution
beta2 <- rbind(beta, c(-2, 2))
df2 <-  ExLinear(beta2, 100, sd = c(.5, 1),
                 xdist = c("runif", "rexp"),
                 xdist.args = list(list(min = -1, max = 2),
                                   list(rate = 3)))

summary(df2)

opar = par("mfrow")
par(mfrow = 1:2)
hist(df2$x1)
hist(df2$x2)
par(opar)

f1 <- flexmix(y ~ ., data = df2, k = 2)

## sort components on Intercept
f1 <- relabel(f1, "model", "Intercept")

## the parameters should be close to the true beta and sigma
round(parameters(f1), 3)
rbind(beta2, sigma)

### A simple Poisson GLM
df3 <- ExLinear(beta/2, 100, min = -1, max = 2, family = "poisson")
plot(y ~ x1, df3, col = attr(df3, "clusters"))

f3 <- flexmix(y ~ ., data = df3, k = 2,
              model = FLXMRglm(family = "poisson"))
round(parameters(relabel(f3, "model", "Intercept")), 3)
beta/2
}
\keyword{datasets}
