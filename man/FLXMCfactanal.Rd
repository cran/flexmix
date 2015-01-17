\name{FLXMCfactanal}
\alias{FLXMCfactanal}
\alias{rFLXM,FLXMCfactanal,FLXcomponent-method}
\title{Driver for Mixtures of Factor Analyzers}
\description{
  This driver for \code{\link{flexmix}} implements estimation of mixtures of
  factor analyzers using ML estimation of factor analysis implemented in
  \code{factanal} in each M-step. 
}
\usage{
FLXMCfactanal(formula = . ~ ., factors = 1, ...)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{\link{flexmix}} using
    \code{\link{update.formula}}. Only the left-hand side (response) of
    the formula is used. Default is to use the original
    \code{\link{flexmix}} model
    formula.}
  \item{factors}{Integer specifying the number of factors in each component.}
  \item{\dots}{Passed to \code{factanal}}
}
\value{
  \code{FLXMCfactanal} returns an object of class \code{FLXM}.
}
\references{
  G. McLachlan and D. Peel. \emph{Finite Mixture Models}, 2000.
  John Wiley and Sons Inc.
}
\author{Bettina Gruen}
\section{Warning}{

  This does not implement the AECM framework presented in McLachlan and
  Peel (2000, p.245), but uses the available functionality in R for ML
  estimation of factor analyzers. The implementation therefore is only
  experimental and has not been well tested.

  Please note that in general a good initialization is crucial for the
  EM algorithm to converge to a suitable solution for this model class.
}
\seealso{\code{\link{flexmix}}}
\examples{
## Reproduce (partly) Table 8.1. p.255 (McLachlan and Peel, 2000)
if (require("gclus")) {
  data("wine", package = "gclus")
  wine_data <- as.matrix(wine[,-1])
  set.seed(123)
  wine_fl_diag <- initFlexmix(wine_data ~ 1, k = 3, nrep = 10,
                              model = FLXMCmvnorm(diagonal = TRUE))
  wine_fl_fact <- lapply(1:4, function(q) flexmix(wine_data ~ 1, model =
                          FLXMCfactanal(factors = q, nstart = 3),
                          cluster = posterior(wine_fl_diag)))
  sapply(wine_fl_fact, logLik)
  ## FULL
  set.seed(123)
  wine_full <- initFlexmix(wine_data ~ 1, k = 3, nrep = 10,
                           model = FLXMCmvnorm(diagonal = FALSE))
  logLik(wine_full)
  ## TRUE
  wine_true <- flexmix(wine_data ~ 1, cluster = wine$Class,
                       model = FLXMCmvnorm(diagonal = FALSE))
  logLik(wine_true)
}
}
\keyword{models}
