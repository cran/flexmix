\name{FLXglm}
\alias{FLXglm}
\title{FlexMix Interface to Generalized Linear Models}
\description{
  This is the main driver for FlexMix interfacing the \code{\link{glm}}
  family of models.
}
\usage{
FLXglm(formula = . ~ ., family = c("gaussian", "binomial", "poisson", "Gamma"))
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{\link{flexmix}} using
    \code{\link{update.formula}}. Default is to use the original
    \code{\link{flexmix}} model
    formula.}
  \item{family}{A character string naming a \code{\link{glm}}
    family function.}
}
\details{
  See \code{\link{flexmix}} for examples.
}
\value{
  Returns an object of class \code{FLXmodel}.
}
\author{Friedrich Leisch}
\references{
  Friedrich Leisch. FlexMix: A general framework for finite mixture
  models and latent class regression in R. Journal of Statistical
  Software, 11(8), 2004. http://www.jstatsoft.org/v11/i08/
}
\seealso{\code{\link{flexmix}}, \code{\link{glm}}}
\keyword{regression}
\keyword{models}
