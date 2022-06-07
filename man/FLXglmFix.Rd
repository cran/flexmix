%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: FLXglmFix.Rd 5229 2022-06-02 14:47:26Z gruen $
%
\name{FLXMRglmfix}
\alias{FLXMRglmfix}
\alias{FLXglmFix}
\title{FlexMix Interface to GLMs with Fixed Coefficients}
\description{
  This implements a driver for FlexMix which interfaces the \code{glm}
  family of models and where it is possible to specify fixed (constant) or
  nested varying coefficients or to ensure that in the Gaussian case the
  variance estimate is equal for all components.
}
\usage{
FLXMRglmfix(formula = . ~ ., fixed = ~0, varFix = FALSE, nested = NULL, 
       family = c("gaussian", "binomial", "poisson", "Gamma"),
       offset = NULL)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{flexmix} using
    \code{update.formula}. Default is to use the original
    \code{flexmix} model
    formula.}
  \item{fixed}{A formula which specifies the additional
    regressors for the fixed (constant) coefficients.}
  \item{varFix}{A logical indicating if the variance estimate for
    Gaussian components should be constrained to be equal for all
    components. It can be also a vector specifying the number of
    components with equal variance.}
  \item{nested}{An object of class \code{FLXnested} or a list specifying
    the nested structure.}
  \item{family}{A character string naming a \code{glm}
    family function.}
  \item{offset}{This can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.}
}
\value{
  Returns an object of class \code{FLXMRglmfix} 
  inheriting from \code{FLXMRglm} and \code{FLXMRfix}.
}
\author{Friedrich Leisch and Bettina Gruen}
\seealso{\code{FLXMRglm}}
\examples{
data("NPreg", package = "flexmix")
ex <- flexmix(yn ~ x | id2, data = NPreg, k = 2,
              cluster = NPreg$class,
              model = FLXMRglm(yn ~ . + I(x^2)))
ex.fix <- flexmix(yn ~ x | id2, data = NPreg,
                  cluster = posterior(ex),
                  model = FLXMRglmfix(nested = list(k = c(1, 1),
                                      formula = c(~0, ~I(x^2)))))
summary(refit(ex))
\dontrun{
summary(refit(ex.fix))
}
}
\keyword{regression}
\keyword{models}
