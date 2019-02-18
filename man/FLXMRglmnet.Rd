\name{FLXMRglmnet}
\alias{FLXMRglmnet}
\alias{FLXMRglmnet-class}

\title{FlexMix Interface for Adaptive Lasso / Elastic Net with GLMs}
\description{

    This is a driver which allows fitting of mixtures of GLMs where
    the coefficients are penalized using the (adaptive) lasso or the
    elastic net by reusing functionality from package \pkg{glmnet}.
    
}
\usage{
    FLXMRglmnet(formula = . ~ ., family = c("gaussian", "binomial", "poisson"),
              adaptive = TRUE, select = TRUE, offset = NULL, ...)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{\link{flexmix}} using
    \code{\link{update.formula}}. Default is to use the original
    \code{\link{flexmix}} model
    formula.}
  \item{family}{A character string naming a \code{\link{glm}}
    family function.}
  \item{adaptive}{A logical indicating if the adaptive lasso should be
    used. By default equal to \code{TRUE}.}
  \item{select}{A logical vector indicating which variables in the
    model matrix should be included in the penalized part.  By default
    equal to \code{TRUE} implying that all variables are penalized.}
  \item{offset}{This can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.}
  \item{\dots}{Additional arguments to be passed to
    \code{\link[glmnet]{cv.glmnet}} fitter.}
}

\details{

  Some care is needed to ensure convergence of the algorithm, which is
  computationally more challenging than a standard EM.  In the proposed
  method, not only are cluster allocations identified and component
  parameters estimated as commonly done in mixture models, but there is
  also variable selection via penalized regression using $k$-fold
  cross-validation to choose the penalty parameter.  For the algorithm
  to converge, it is necessary that the same cross-validation
  partitioning be used across the EM iterations, i.e., the subsamples
  for cross-validation must be defined at the beginning This is
  accomplished using the \code{foldid} option as an additional parameter
  to be passed to \code{\link[glmnet]{cv.glmnet}} (see
  \pkg{glmnet} package documentation).
}

\value{
  Returns an object of class \code{FLXMRglm}.
}
\author{
  Frederic Mortier and Nicolas Picard.
}
\seealso{
  \code{\link{FLXMRglm}}
}
\references{
  
    Frederic Mortier, Dakis-Yaoba Ouedraogo, Florian Claeys, Mahlet
    G. Tadesse, Guillaume Cornu, Fidele Baya, Fabrice Benedet, Vincent
    Freycon, Sylvie Gourlet-Fleury and Nicolas Picard.
    Mixture of inhomogeneous matrix models for species-rich ecosystems.
    \emph{Environmetrics}, \bold{26}(1), 39-51, 2015.
    doi:10.1002/env.2320

}

\examples{
    set.seed(12)
    p <- 10
    beta <- matrix(0, nrow = p + 1, ncol = 2)
    beta[1,] <- c(-1, 1)
    beta[cbind(c(5, 10), c(1, 2))] <- 1

    nobs <- 100
    X <- matrix(rnorm(nobs * p), nobs, p)
    mu <- cbind(1, X) \%*\% beta
    z <- sample(1:ncol(beta), nobs, replace = TRUE)
    y <- mu[cbind(1:nobs, z)] + rnorm(nobs)
    data <- data.frame(y, X)
    ## The maximum number of iterations is reduced to
    ## avoid a long running time.
    fo <- sample(rep(seq(10), length = nrow(data)))
    ex1 <- flexmix(y ~ ., data = data, k = 2, cluster = z,
                   model = FLXMRglmnet(foldid = fo),
                   control = list(iter.max = 2))
    parameters(ex1)
}
\keyword{regression}
\keyword{cluster}
