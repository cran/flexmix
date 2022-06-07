%
%  Copyright (C) 2004-2009 Friedrich Leisch and Bettina Gruen
%  $Id: EIC.Rd 3912 2008-03-13 15:10:24Z gruen $
%
\name{EIC}
\alias{EIC}
\alias{EIC,flexmix-method}
\alias{EIC,stepFlexmix-method}
\title{Entropic Measure Information Criterion}
\description{
  Compute the entropic measure information criterion for model selection.
}
\usage{
\S4method{EIC}{flexmix}(object, \dots)
\S4method{EIC}{stepFlexmix}(object, \dots)
}
\arguments{
  \item{object}{See Methods section below}
  \item{\dots}{Some methods for this generic function may take additional,
          optional arguments. At present none do.}
}
\section{Methods}{
  \describe{
    \item{object = "flexmix":}{Compute the EIC of a \code{flexmix} object.}
    \item{object = "stepFlexmix":}{Compute the EIC of all
      models contained in the \code{stepFlexmix} object.}
}}
\value{
  Returns a numeric vector with the corresponding EIC value(s). 
}
\keyword{methods}
\author{Bettina Gruen}
\references{
  V. Ramaswamy, W. S. DeSarbo, D. J. Reibstein, and W. T. Robinson. An
  empirical pooling approach for estimating marketing mix elasticities
  with PIMS data. \emph{Marketing Science}, \bold{12}(1), 103--124, 1993.
}
\examples{
data("NPreg", package = "flexmix")
ex1 <- flexmix(yn ~ x + I(x^2), data = NPreg, k = 2)
EIC(ex1)
}
