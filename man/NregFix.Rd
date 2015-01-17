%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: NregFix.Rd 5008 2015-01-13 20:30:25Z gruen $
%  
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
    \item{\code{x1}}{Independent variable: numeric variable.}
    \item{\code{x2}}{Independent variable: a factor with two levels:
      \code{0} and \code{1}.}
    \item{\code{w}}{Concomitant variable: a factor with two levels:
      \code{0} and \code{1}.}
    \item{\code{y}}{Dependent variable.}
    \item{\code{class}}{Latent class memberships.}
  }
}
\examples{
data("NregFix", package = "flexmix")
library("lattice")
xyplot(y ~ x1 | x2 * w, data = NregFix, groups = class)
Model <- FLXMRglmfix(~ 1, fixed = ~ x2, 
                     nested = list(k = c(2, 1),
                     formula = c(~x1, ~0)))
fittedModel <- initFlexmix(y ~ 1, model = Model, data = NregFix, k = 3,
                           concomitant = FLXPmultinom(~ w), nrep = 5)
fittedModel
}
\keyword{datasets}
