\name{compatibility-aliases}
\alias{FLXglm}
\alias{FLXglmFix}
\alias{FLXmclust}
\alias{FLXbclust}
\alias{FLXconstant}
\alias{FLXmultinom}
\title{Backward Compatibility Aliases for Model Drivers}
\description{
  These functions are renamed versions of model drivers included for 
  backward compatibility.
}
\usage{
FLXglm(formula = . ~ .,
       family = c("gaussian", "binomial", "poisson", "Gamma"),
       offset = NULL)
FLXglmFix(formula = . ~ ., fixed = ~0, varFix = FALSE, nested = NULL, 
          family = c("gaussian", "binomial", "poisson", "Gamma"),
          offset = NULL)
FLXmclust(formula = . ~ ., diagonal = TRUE)
FLXbclust(formula = . ~ ., truncated = FALSE)
FLXconstant()
FLXmultinom(formula = ~1)
}
\arguments{
  \item{formula}{For component model drivers, a formula which is
    interpreted relative to the formula specified in the call to
    \code{\link{flexmix}} using \code{\link{update.formula}}. Default is
    to use the original \code{\link{flexmix}} model formula.  For
    concomitant model drivers, a formula for determining the model
    matrix of the concomitant variables.}
  \item{family}{A character string naming a \code{\link{glm}}
    family function.}
  \item{offset}{This can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.}
  \item{fixed}{A formula which specifies the additional
    regressors for the fixed (constant) coefficients.}
  \item{varFix}{A logical indicating if the variance estimate for
    Gaussian components should be constrained to be equal for all
    components. It can be also a vector specifying the number of
    components with equal variance.}
  \item{nested}{An object of class \code{FLXnested} or a list specifying
    the nested structure.}
  \item{diagonal}{If \code{TRUE}, then the covariance matrix of the
    components is restricted to diagonal matrices.}
  \item{truncated}{A logical, if \code{TRUE} the observations for the
    pattern with only zeros are missing and the truncated likelihood is
    optimized using an EM-algorithm.}
  
}
\details{
  The functions map to the current implementations as follows:
  \itemize{
    \item \code{FLXglm} aliases \code{\link{FLXMRglm}}.
    \item \code{FLXglmFix} aliases \code{\link{FLXMRglmfix}}.
    \item \code{FLXmclust} aliases \code{\link{FLXMCmvnorm}}.
    \item \code{FLXbclust} aliases \code{\link{FLXMCmvbinary}}.
    \item \code{FLXmultinom} aliases \code{\link{FLXPmultinom}}.
    \item \code{FLXconstant} aliases \code{\link{FLXPconstant}}.
  }
}
\keyword{internal}
