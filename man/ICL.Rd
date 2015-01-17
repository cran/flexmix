%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: ICL.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{ICL}
\alias{ICL,flexmix-method}
\alias{ICL,stepFlexmix-method}
\title{Integrated Completed Likelihood Criterion}
\description{
  Compute the Integrated Completed Likelihood criterion for model selection.
}
\usage{
\S4method{ICL}{flexmix}(object, \dots)
\S4method{ICL}{stepFlexmix}(object, \dots)
}
\arguments{
  \item{object}{see Methods section below}
  \item{\dots}{Some methods for this generic function may take additional,
          optional arguments. At present none do.}
}
\section{Methods}{
  \describe{
    \item{object = "flexmix":}{Compute the ICL of a \code{flexmix} object.}
    \item{object = "stepFlexmix":}{Compute the ICL of all
      models contained in the \code{stepFlexmix} object.}
}}
\value{
  Returns a numeric vector with the corresponding ICL value(s). 
}
\keyword{methods}
\author{Friedrich Leisch and Bettina Gruen}
\references{
  C. Biernacki, G. Celeux, and G. Govaert. Assessing a mixture model
  for clustering with the integrated completed likelihood. \emph{IEEE
  Transactions on Pattern Analysis and Machine Intelligence},
  \emph{22}(7), 719--725, 2000.
}
\examples{
data("NPreg", package = "flexmix")
ex1 <- flexmix(yn ~ x + I(x^2), data = NPreg, k = 2)
ICL(ex1)
}
