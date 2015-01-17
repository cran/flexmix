\name{boot}
\alias{boot}
\alias{boot,flexmix-method}
\alias{LR_test}
\alias{LR_test,flexmix-method}
\alias{boot,flexmix-method}
\alias{show,FLXboot-method}
\alias{FLXboot-class}
\alias{plot,FLXboot,missing-method}
\alias{parameters,FLXboot-method}
\alias{clusters,FLXboot,listOrdata.frame-method}
\alias{predict,FLXboot-method}
\alias{posterior,FLXboot,listOrdata.frame-method}
\title{Bootstrap a flexmix Object}
\description{
  Given a \code{flexmix} object perform parametric or empirical bootstrap.
}
\usage{
boot(object, ...)
\S4method{boot}{flexmix}(object, R, sim = c("ordinary", "empirical", "parametric"),
    initialize_solution = FALSE, keep_weights = FALSE,
    keep_groups = TRUE, verbose = 0, control,
    k, model = FALSE, ...)
LR_test(object, ...)
\S4method{LR_test}{flexmix}(object, R, alternative = c("greater", "less"), control, ...)
}
\arguments{
  \item{object}{A fitted finite mixture model of class \code{flexmix}.} 
  \item{R}{The number of bootstrap replicates.}
  \item{sim}{A character string indicating the type of simulation
    required.  Possible values are \code{"ordinary"} (the default),
    \code{"parametric"}, or \code{"empirical"}.}
  \item{initialize_solution}{A logical. If \code{TRUE} the EM algorithm is
    initialized in the given solution.}
  \item{keep_weights}{A logical. If \code{TRUE} the weights are kept.}
  \item{keep_groups}{A logical. If \code{TRUE} the groups are kept.}
  \item{verbose}{If a positive integer, then progress information
    is reported every \code{verbose} iterations. If 0,
    no output is generated during the bootstrap replications.}
  \item{control}{Object of class \code{FLXcontrol} or a named list.
    If missing the control of the fitted \code{object} is taken.}
  \item{k}{Vector of integers specifying for which number of components
    finite mixtures are fitted to the bootstrap samples.  If missing the
    number of components of the fitted \code{object} are taken.}
  \item{alternative}{A character string specifying the alternative
    hypothesis, must be either \code{"greater"} (default) or
    \code{"less"} indicating if the alternative hypothesis is that the
    mixture has one more component or one less.}
  \item{model}{A logical. If \code{TRUE} the model and the weights slot for
    each sample are stored and returned.}
  \item{\dots}{Further arguments to be passed to or from methods.}
}
\value{
  \code{boot} returns an object of class \code{FLXboot} which
  contains the fitted parameters, the fitted priors, the log
  likelihoods, the number of components of the fitted mixtures and the
  information if the EM algorithm has converged.

  \code{LR_test} returns an object of class \code{htest} containing the
  number of valid bootstrap replicates, the p-value, the - twice log
  likelihood ratio test statistics for the original data and the
  bootstrap replicates.
}
\author{Bettina Gruen}
\examples{
data("NPreg", package = "flexmix")
fitted <- initFlexmix(yn ~ x + I(x^2) | id2, data = NPreg, k = 2)
\dontrun{
lrtest <- LR_test(fitted, alternative = "greater", R = 20,
                  verbose = 1)
}
}
\keyword{methods}


