%
%  Copyright (C) 2008 Friedrich Leisch and Bettina Gruen
%  $Id: FLXMRrobglm.Rd 5229 2022-06-02 14:47:26Z gruen $
%
\name{FLXMRrobglm}
\alias{FLXMRrobglm}
\alias{FLXMRrobglm-class}
\title{FlexMix Driver for Robust Estimation of Generalized Linear Models}
\description{
  This driver adds a noise component to the mixture model which can be
  used to model background noise in the data. See the Compstat paper
  Leisch (2008) cited below for details.
}
\usage{
FLXMRrobglm(formula = . ~ ., family = c("gaussian", "poisson"),
            bgw = FALSE, ...)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{flexmix} using
    \code{\link{update.formula}}. Default is to use the original
    \code{flexmix} model
    formula.}
  \item{family}{A character string naming a \code{\link{glm}}
    family function.}
  \item{bgw}{Logical, controls whether the parameters of the background
    component are fixed to multiples of location and scale of the
    complete data (the default), or estimated by EM with normal weights for the
    background (\code{bgw = TRUE}).}
  \item{\dots}{passed to \code{FLXMRglm}}
}
\value{
  Returns an object of class \code{FLXMRrobglm} inheriting from \code{FLXMRglm}.
}
\author{Friedrich Leisch and Bettina Gruen}
\note{
  The implementation of this model class is currently under development,
  and some methods like \code{refit} are still missing.
}
\references{
  Friedrich Leisch. Modelling background noise in finite mixtures of
  generalized linear regression models. In Paula Brito, editor, Compstat
  2008--Proceedings in Computational Statistics, 385--396. Physica
  Verlag, Heidelberg, Germany, 2008.\cr
  Preprint available at http://epub.ub.uni-muenchen.de/6332/.
}
\examples{
## Example from Compstat paper, see paper for detailed explanation:
data("NPreg", package = "flexmix")
DATA <- NPreg[, 1:2]
set.seed(3)
DATA2 <- rbind(DATA, cbind(x = -runif(3), yn = 50 + runif(3)))

## Estimation without (f2) and with (f3) background component
f2 <- flexmix(yn ~ x + I(x^2), data = DATA2, k = 2)
f3 <- flexmix(yn ~ x + I(x^2), data = DATA2, k = 3,
              model = FLXMRrobglm(), 
              control = list(minprior = 0))

## Predict on new data for plots
x <- seq(-5,15, by = .1)
y2 <- predict(f2, newdata = data.frame(x = x))
y3 <- predict(f3, newdata = data.frame(x = x))

## f2 was estimated without background component:
plot(yn ~ x, data = DATA2, pch = clusters(f2), col = clusters(f2))
lines(x, y2$Comp.1, col = 1)
lines(x, y2$Comp.2, col = 2)

## f3 is with background component:
plot(yn ~ x, data = DATA2, pch = 4 - clusters(f3),
     col = 4 - clusters(f3))
lines(x, y3$Comp.2, col = 2)
lines(x, y3$Comp.3, col = 1)
}
\keyword{models}
