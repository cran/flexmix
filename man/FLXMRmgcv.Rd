\name{FLXMRmgcv}
\alias{FLXMRmgcv}
\alias{FLXMRmgcv-class}

\title{FlexMix Interface to GAMs}
\description{
  This is a driver which allows fitting of mixtures of GAMs.
}
\usage{
FLXMRmgcv(formula = . ~ ., family = c("gaussian", "binomial", "poisson"),
          offset = NULL, control = NULL, optimizer = c("outer", "newton"),
          in.out = NULL, eps = .Machine$double.eps, ...)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{\link{flexmix}} using
    \code{\link{update.formula}}. Default is to use the original
    \code{\link{flexmix}} model
    formula.}
  \item{family}{A character string naming a \code{\link{glm}}
    family function.}
  \item{offset}{This can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.}
  \item{control}{A list of fit control parameters returned by
    \code{gam.control}.}
  \item{optimizer}{An array specifying the numerical optimization method
    to use to optimize the smoothing parameter estimation criterion; for
    more details see \code{\link[mgcv]{gam}}.}
  \item{in.out}{Optional list for initializing outer iteration; for more
    details see \code{\link[mgcv]{gam}}.}
  \item{eps}{Observations with an a-posteriori probability smaller or
    equal to \code{eps} are omitted in the M-step.}
  \item{\dots}{Additional arguments to be pased to the GAM fitter.}
}
\value{
  Returns an object of class \code{FLXMRmgcv} inheriting from \code{FLXMRglm}.
}
\author{
  Bettina Gruen
}
\seealso{
  \code{\link{FLXMRglm}}
}
\examples{
set.seed(2012)
x <- seq(0, 1, length.out = 100)
z <- sample(0:1, length(x), replace = TRUE)
y <- rnorm(length(x), ifelse(z, 5 * sin(x * 2 * pi), 10 * x - 5))
fitted_model <- flexmix(y ~ s(x), model = FLXMRmgcv(),
                        cluster = z + 1,
                        control = list(tolerance = 10^-3))
plot(y ~ x, col = clusters(fitted_model))
matplot(x, fitted(fitted_model), type = "l", add = TRUE)
}
\keyword{regression}
\keyword{cluster}
