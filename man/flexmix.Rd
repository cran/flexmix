%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: flexmix.Rd 5115 2017-04-07 08:18:13Z gruen $
%
\name{flexmix}
\alias{flexmix}
\alias{flexmix,formula,ANY,ANY,ANY,missing-method}
\alias{flexmix,formula,ANY,ANY,ANY,list-method}
\alias{flexmix,formula,ANY,ANY,ANY,FLXM-method}
\alias{prior,flexmix-method}
\alias{show,flexmix-method}
\alias{summary,flexmix-method}
\alias{show,summary.flexmix-method}
\title{Flexible Mixture Modeling}
\description{
  FlexMix implements a general framework for finite
  mixtures of regression models. Parameter estimation is performed using
  the EM algorithm: the E-step is implemented by \code{flexmix}, while
  the user can specify the M-step.
}
\usage{
flexmix(formula, data = list(), k = NULL, cluster = NULL, 
        model = NULL, concomitant = NULL, control = NULL,
        weights = NULL)
\S4method{summary}{flexmix}(object, eps = 1e-4, ...)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit. The
    general form is \code{y~x|g} where \code{y} is the response,
    \code{x} the set of predictors and \code{g} an optional grouping
    factor for repeated measurements.}
  \item{data}{An optional data frame containing the variables in the model.}
  \item{k}{Number of clusters (not needed if \code{cluster} is specified).}
  \item{cluster}{Either a matrix with \code{k} columns of initial
    cluster membership probabilities for each observation; or a factor
    or integer vector with the initial cluster
    assignments of observations at the start of the EM
    algorithm. Default is random assignment into \code{k}
    clusters.}
  \item{weights}{An optional vector of replication weights to be used in
    the fitting process.  Should be \code{NULL}, an integer vector or a
    formula.}
  \item{model}{Object of class \code{FLXM} or list of \code{FLXM}
    objects. Default is the object returned by calling
    \code{\link{FLXMRglm}()}.}
  \item{concomitant}{Object of class \code{FLXP}. Default is the
    object returned by calling \code{\link{FLXPconstant}}.}
  \item{control}{Object of class \code{FLXcontrol} or a named list.}
  \item{object}{Object of class \code{flexmix}.}
  \item{eps}{Probabilities below this threshold are treated as zero in the
    summary method.}
  \item{\dots}{Currently not used.}
}
\details{
  FlexMix models are described by objects of class \code{FLXM},
  which in turn are created by driver functions like
  \code{\link{FLXMRglm}} or \code{\link{FLXMCmvnorm}}. Multivariate
  responses with independent components can be specified using a
  list of \code{FLXM} objects.

  The \code{summary} method lists for each component the prior
  probability, the number of observations assigned to the corresponding
  cluster, the number of observations with a posterior probability
  larger than \code{eps} and the ratio of the latter two numbers (which
  indicates how separated the cluster is from the others).
}
\value{
  Returns an object of class \code{flexmix}.
}
\author{Friedrich Leisch and Bettina Gruen}
\seealso{\code{\link[flexmix]{plot-methods}}}
\references{
  Friedrich Leisch. FlexMix: A general framework for finite mixture
  models and latent class regression in R. \emph{Journal of Statistical
    Software}, \bold{11}(8), 2004. doi:10.18637/jss.v011.i08

  Bettina Gruen and Friedrich Leisch. Fitting finite mixtures of
  generalized linear regressions in R. Computational Statistics & Data
  Analysis, 51(11), 5247-5252, 2007. doi:10.1016/j.csda.2006.08.014

  Bettina Gruen and Friedrich Leisch. FlexMix Version 2: Finite
  mixtures with concomitant variables and varying and constant
  parameters Journal of Statistical Software, 28(4), 1-35, 2008.
  doi:10.18637/jss.v028.i04
  
}
\keyword{regression}
\keyword{cluster}
\examples{

data("NPreg", package = "flexmix")

## mixture of two linear regression models. Note that control parameters
## can be specified as named list and abbreviated if unique.
ex1 <- flexmix(yn ~ x + I(x^2), data = NPreg, k = 2,
               control = list(verb = 5, iter = 100))

ex1
summary(ex1)
plot(ex1)

## now we fit a model with one Gaussian response and one Poisson
## response. Note that the formulas inside the call to FLXMRglm are
## relative to the overall model formula.
ex2 <- flexmix(yn ~ x, data = NPreg, k = 2,
               model = list(FLXMRglm(yn ~ . + I(x^2)), 
                            FLXMRglm(yp ~ ., family = "poisson")))
plot(ex2)

ex2
table(ex2@cluster, NPreg$class)

## for Gaussian responses we get coefficients and standard deviation
parameters(ex2, component = 1, model = 1)

## for Poisson response we get only coefficients
parameters(ex2, component = 1, model = 2)

## fitting a model only to the Poisson response is of course
## done like this
ex3 <- flexmix(yp ~ x, data = NPreg, k = 2,
               model = FLXMRglm(family = "poisson"))

## if observations are grouped, i.e., we have several observations per
## individual, fitting is usually much faster:
ex4 <- flexmix(yp~x|id1, data = NPreg, k = 2,
               model = FLXMRglm(family = "poisson"))

## And now a binomial example. Mixtures of binomials are not generically
## identified, here the grouping variable is necessary:
set.seed(1234)
ex5 <- initFlexmix(cbind(yb,1 - yb) ~ x, data = NPreg, k = 2,
                   model = FLXMRglm(family = "binomial"), nrep = 5)
table(NPreg$class, clusters(ex5))

ex6 <- initFlexmix(cbind(yb, 1 - yb) ~ x | id2, data = NPreg, k = 2,
                   model = FLXMRglm(family = "binomial"), nrep = 5)
table(NPreg$class, clusters(ex6))

}

