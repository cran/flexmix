\name{NregFix}
\alias{NregFix}
\title{Artificial Example for Normal Regression}
\description{
  A simple artificial regression example with 3 latent classes, two
  independent variables, one concomitant variable and a dependent
  variable which follows a Gaussian distribution.  
}
\usage{
data("NregFix")
}
\format{
  A data frame with 200 observations on the following 5 variables.
  \describe{
    \item{\code{x1}}{independent variable, numeric variable}
    \item{\code{x2}}{independent variable, a factor with two levels \code{0} \code{1}}
    \item{\code{w}}{concomitant variable, a factor with two levels \code{0} \code{1}}
    \item{\code{y}}{dependent variable}
    \item{\code{class}}{latent class memberships}
  }
}
\examples{
data("NregFix")
library("lattice")
xyplot(y ~ x1 | x2*w, data = NregFix, groups = class)
Model <- FLXMRglmfix(~1, fixed = ~x2, 
                   nested = list(k = c(2,1), formula = c(~x1, ~0)))
fittedModel <- stepFlexmix(y ~ 1, model = Model, data = NregFix, k = 3,
                          concomitant = FLXPmultinom(~ w), nrep=5)
fittedModel
summary(refit(fittedModel))
}
\keyword{datasets}
