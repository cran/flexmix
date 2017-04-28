%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: refit.Rd 5115 2017-04-07 08:18:13Z gruen $
%
\name{refit-methods}
\alias{refit,flexmix-method}
\alias{FLXRmstep-class}
\alias{FLXRoptim-class}
\alias{show,FLXR-method}
\alias{show,Coefmat-method}
\alias{summary,FLXRoptim-method}
\alias{summary,FLXRmstep-method}
\alias{plot,FLXRoptim,missing-method}
\title{Refit a Fitted Model}
\description{
  Refits an estimated flexmix model to obtain additional information
  like coefficient significance p-values for GLM regression.
}
\usage{
\S4method{refit}{flexmix}(object, newdata, method = c("optim",
"mstep"), ...)
\S4method{summary}{FLXRoptim}(object, model = 1, which = c("model",
"concomitant"), ...)
\S4method{summary}{FLXRmstep}(object, model = 1, which = c("model",
"concomitant"), ...)

\S4method{plot}{FLXRoptim,missing}(x, y, model = 1, which = c("model", "concomitant"),
         bycluster = TRUE, alpha = 0.05, components, labels = NULL,
         significance = FALSE, xlab = NULL, ylab = NULL, ci = TRUE,
         scales = list(), as.table = TRUE, horizontal = TRUE, ...)
}
\arguments{
  \item{object}{An object of class \code{"flexmix"}}
  \item{newdata}{Optional new data.}
  \item{method}{Specifies if the variance covariance matrix is
    determined using \code{\link{optim}} or if the posteriors are
    assumed as given and an M-step is performed.}
  \item{model}{The model (for a multivariate response) that shall be
    used.}
  \item{which}{Specifies if a component specific model or the
    concomitant variable model is used.}
  \item{x}{An object of class \code{"FLXRoptim"}}
  \item{y}{Missing object.}
  \item{bycluster}{A logical if the parameters should be group by
    cluster or by variable.}
  \item{alpha}{Numeric indicating the significance level.}
  \item{components}{Numeric vector specifying which components are
    plotted. The default is to plot all components.}
  \item{labels}{Character vector specifying the variable names used.}
  \item{significance}{A logical indicating if non-significant
    coefficients are shaded in a lighter grey.}
  \item{xlab}{String for the x-axis label.}
  \item{ylab}{String for the y-axis label.}
  \item{ci}{A logical indicating if significant and insignificant
    parameter estimates are shaded differently.}
  \item{scales}{See argument of the same name for
    function \code{\link[lattice]{xyplot}}.}
  \item{as.table}{See arguments of the same name for
    function \code{\link[lattice]{xyplot}}.}
  \item{horizontal}{See arguments of the same name for
    function \code{\link[lattice]{xyplot}}.}
  \item{\dots}{Currently not used}
}
\value{
  An object inheriting form class \code{FLXR} is returned. For the
  method using \code{optim} the object has class \code{FLXRoptim} and
  for the M-step method it has class \code{FLXRmstep}. Both classes give
  similar results for their \code{summary} methods. 
  Objects of class \code{FLXRoptim} have their own \code{plot} method.
  \code{Lapply} can be used to further analyse the refitted component
  specific models of objects of class \code{FLXRmstep}.
}
\details{
  The \code{refit} method for \code{FLXMRglm} models in
  combination with the \code{summary} method can be
  used to obtain the usual tests for significance of coefficients. Note
  that the tests are valid only if \code{flexmix} returned the maximum
  likelihood estimator of the parameters.  If \code{refit} is used with
  \code{method = "mstep"} for these component specific models the
  returned object contains a \code{glm} object for each component where
  the elements \code{model} which is the model frame and \code{data}
  which contains the original dataset are missing.
}
\keyword{methods}
\author{Friedrich Leisch and Bettina Gruen}
\references{
  Friedrich Leisch. FlexMix: A general framework for finite mixture
  models and latent class regression in R. \emph{Journal of Statistical
  Software}, \bold{11}(8), 2004. doi:10.18637/jss.v011.i08
}
\section{Warning}{
  For \code{method = "mstep"} the standard deviations are determined
  separately for each of the components using the a-posteriori
  probabilities as weights without accounting for the fact that the
  components have been simultaneously estimated. The derived standard
  deviations are hence approximative and should only be used in an
  exploratory way, as they are underestimating the uncertainty given
  that the missing information of the component memberships are replaced
  by the expected values.

  The \code{newdata} argument can only be specified when using
  \code{method = "mstep"} for refitting \code{FLXMRglm} components. A
  variant of \code{glm} for weighted ML estimation is used for fitting
  the components and full \code{glm} objects are returned. Please note
  that in this case the data and the model frame are stored for each
  component which can significantly increase the object size.  }

\examples{

data("NPreg", package = "flexmix")
ex1 <- flexmix(yn ~ x + I(x^2), data = NPreg, k = 2)
ex1r <- refit(ex1)

## in one component all coefficients should be highly significant,
## in the other component only the linear term
summary(ex1r)
}

