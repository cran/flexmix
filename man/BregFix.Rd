\name{BregFix}
\alias{BregFix}
\docType{data}
\title{Artificial Example for Binomial Regression}
\description{
	A simple artificial regression example data set with 3 latent
	classes, one independent variable \code{x} and a concomitant
	variable \code{w}.
}
\usage{data("BregFix")}
\format{
  A data frame with 200 observations on the following 5 variables.
  \describe{
    \item{\code{yes}}{number of successes}
    \item{\code{no}}{number of failures}
    \item{\code{x}}{independent variable}
    \item{\code{w}}{concomitant variable, a factor with levels \code{0} \code{1}}
    \item{\code{class}}{latent class memberships}
  }
}
\examples{
data("BregFix")
Model <- FLXMRglmfix(family="binomial",
                   nested = list(formula = c(~x, ~0), k = c(2,1)))
Conc <- FLXPmultinom(~w)
FittedBin <- stepFlexmix(cbind(yes, no) ~ 1, data = BregFix,
                         k = 3, model = Model, concomitant = Conc)
summary(FittedBin)
}
\keyword{datasets}
