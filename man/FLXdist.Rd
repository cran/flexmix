\name{FLXdist}
\alias{FLXdist}
\alias{simulate,FLXdist-method}
\alias{show,FLXdist-method}
\title{Finite Mixtures of Distributions}
\description{
  Constructs objects of class \code{FLXdist} which represent unfitted finite mixture
  models.
}
\usage{
 FLXdist(formula, k = NULL, model = FLXMRglm(), components,
         concomitant = FLXPconstant())
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit. The
    general form is \code{~x|g} where \code{x} is the set of predictors
    and \code{g} an optional grouping factor for repeated measurements.}
  \item{k}{Integer specifying the number of cluster or a numeric vector of
    length equal to the length of components, specifying the prior
    probabilities of clusters.}
  \item{model}{Object of class \code{FLXM} or a list of \code{FLXM}
    objects. Default is the object returned by calling
    \code{FLXMRglm()}.}
  \item{components}{A list of length equal to the number of components
    containing a list of length equal to the number of models which
    again contains a list of named elements for defining the parameters
    of the component-specific model.}
  \item{concomitant}{Object of class \code{FLXconcomitant} specifying
    the model for concomitant variables.}
}
\value{
  Returns an object of class \code{FLXdist}.
}
\author{Bettina Gruen}
\seealso{\code{FLXdist-class}}
\keyword{utilities}
