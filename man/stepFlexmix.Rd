%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: stepFlexmix.Rd 5115 2017-04-07 08:18:13Z gruen $
%
\name{stepFlexmix}
\alias{stepFlexmix}
\alias{initFlexmix}
\alias{initMethod}
\alias{stepFlexmix-class}
\alias{initMethod-class}
\alias{plot,stepFlexmix,missing-method}
\alias{show,stepFlexmix-method}
\alias{getModel,stepFlexmix-method}
\alias{unique,stepFlexmix-method}
\title{Run FlexMix Repeatedly}
\description{
  Runs flexmix repeatedly for different numbers of components and returns
  the maximum likelihood solution for each.
}
\usage{
initFlexmix(..., k, init = list(), control = list(), nrep = 3L,
            verbose = TRUE, drop = TRUE, unique = FALSE)
initMethod(name = c("tol.em", "cem.em", "sem.em"),
           step1 = list(tolerance = 10^-2),
           step2 = list(), control = list(), nrep = 3L)

stepFlexmix(..., k = NULL, nrep = 3, verbose = TRUE, drop = TRUE,
            unique = FALSE)

\S4method{plot}{stepFlexmix,missing}(x, y, what = c("AIC", "BIC", "ICL"),
  xlab = NULL, ylab = NULL, legend = "topright", ...)

\S4method{getModel}{stepFlexmix}(object, which = "BIC")

\S4method{unique}{stepFlexmix}(x, incomparables = FALSE, ...)
}
\arguments{
  \item{\dots}{Passed to \code{\link{flexmix}} (or \code{\link{matplot}}
    in the \code{plot} method).}
  \item{k}{A vector of integers passed in turn to the \code{k} argument
    of \code{\link{flexmix}}.}
  \item{init}{An object of class \code{"initMethod"} or a named list
    where \code{initMethod} is called with it as arguments in addition
    to the \code{control} argument.}
  \item{name}{A character string indication which initialization
    strategy should be employed: short runs of EM followed by a long
    (\code{"tol.em"}), short runs of CEM followed by a long EM run
    (\code{"cem.em"}), short runs of SEM followed by a long EM run
    (\code{"sem.em"}).}
  \item{step1}{A named list which combined with the \code{control}
    argument is coercable to a \code{"FLXcontrol"} object. This control
    setting is used for the short runs.}
  \item{step2}{A named list which combined with the \code{control}
    argument is coercable to a \code{"FLXcontrol"} object. This control
    setting is used for the long run.}
  \item{control}{A named list which combined with the \code{step1} or
    the \code{step2} argument is coercable to a \code{"FLXcontrol"}
    object.}  
  \item{nrep}{For each value of \code{k} run \code{\link{flexmix}}
    \code{nrep} times and keep only the solution with maximum
    likelihood. If \code{nrep} is set for the long run, it is ignored,
    because the EM algorithm is deterministic using the best solution
    discovered in the short runs for initialization.}
  \item{verbose}{If \code{TRUE}, show progress information during
    computations.}
  \item{drop}{If \code{TRUE} and \code{k} is of length 1, then a single
    flexmix object is returned instead of a \code{"stepFlexmix"}
    object.}
  \item{unique}{If \code{TRUE}, then \code{unique()} is called on the
    result, see below.}
  \item{x, object}{An object of class \code{"stepFlexmix"}.}
  \item{y}{Not used.}
  \item{what}{Character vector naming information criteria to
    plot. Functions of the same name must exist, which take a
    \code{stepFlexmix} object as input and return a numeric vector like
    \code{AIC,stepFlexmix-method} (see examples below).}
  \item{xlab,ylab}{Graphical parameters.}
  \item{legend}{If not \code{FALSE} and \code{what} contains more
    than 1 element, a legend is placed at the specified location, see
    \code{\link{legend}} for details.}
  \item{which}{Number of model to get. If character, interpreted as
    number of components or name of an information criterion.}
  \item{incomparables}{A vector of values that cannot be
          compared. Currently, \code{FALSE} is the only possible value,
          meaning that all values can be compared.}
}
\value{
  An object of class \code{"stepFlexmix"} containing the best models
  with respect to the log likelihood for the different number of
  components in a slot if \code{length(k)>1}, else directly an object of
  class \code{"flexmix"}.

  If \code{unique = FALSE}, then the resulting object contains one
  model per element of \code{k} (which is the number of clusters the EM
  algorithm started with). If \code{unique = TRUE}, then the result
  is resorted according to the number of clusters contained in the
  fitted models (which may be less than the number with which the EM
  algorithm started), and only the maximum likelihood solution for each
  number of fitted clusters is kept. This operation can also be done
  manually by calling \code{unique()} on objects of class
  \code{"stepFlexmix"}. 
}
\author{Friedrich Leisch and Bettina Gruen}
\references{
  Friedrich Leisch. FlexMix: A general framework for finite mixture
  models and latent class regression in R. \emph{Journal of Statistical
  Software}, \bold{11}(8), 2004. doi:10.18637/jss.v011.i08

  Christophe Biernacki, Gilles Celeux and Gerard Govaert. Choosing
  starting values for the EM algorithm for getting the highest
  likelihood in multivariate Gaussian mixture models. \emph{Computational
  Statistics & Data Analysis}, \bold{41}(3--4), 561--575, 2003.

  Theresa Scharl, Bettina Gruen and Friedrch Leisch.  Mixtures of
  regression models for time-course gene expression data: Evaluation of
  initialization and random effects. \emph{Bioinformatics},
  \bold{26}(3), 370--377, 2010.
}
\examples{
data("Nclus", package = "flexmix")

## try 2 times for k = 4
set.seed(511)
ex1 <- initFlexmix(Nclus~1, k = 4, model = FLXMCmvnorm(diagonal = FALSE),
                   nrep = 2)
ex1

## now 2 times each for k = 2:5, specify control parameter
ex2 <- initFlexmix(Nclus~1, k = 2:5, model = FLXMCmvnorm(diagonal = FALSE),
                   control = list(minprior = 0), nrep = 2)
ex2
plot(ex2)

## get BIC values
BIC(ex2)

## get smallest model
getModel(ex2, which = 1)

## get model with 3 components
getModel(ex2, which = "3")

## get model with smallest ICL (here same as for AIC and BIC: true k = 4)
getModel(ex2, which = "ICL")

## now 1 time each for k = 2:5, with larger minimum prior
ex3 <- initFlexmix(Nclus~1, k = 2:5,
                   model = FLXMCmvnorm(diagonal = FALSE),
                   control = list(minprior = 0.1), nrep = 1)
ex3

## keep only maximum likelihood solution for each unique number of
## fitted clusters:
unique(ex3)
}
\keyword{cluster}
\keyword{regression}
