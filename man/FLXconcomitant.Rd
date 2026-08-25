%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: FLXconcomitant.Rd 5346 2026-08-25 06:32:19Z gruen $
%
\name{FLXP}
\docType{class}
\alias{FLXPconstant}
\alias{FLXPmultinom}
\alias{show,FLXP-method}
\title{Creates the Concomitant Variable Model}
\description{
  Creator functions for the concomitant variable
  model. \code{FLXPconstant} specifies constant priors
  and \code{FLXPmultinom} multinomial logit models for the priors.  
}
\usage{
FLXPconstant()
FLXPmultinom(formula = ~1)
}
\arguments{
  \item{formula}{A formula for determining the model matrix of the
    concomitant variables.}
}
\details{
  \code{FLXPmultinom} uses \code{nnet.default} from \pkg{nnet} to fit the
  multinomial logit model.
}
\value{
  Object of class \code{FLXP}. \code{FLXPmultinom} returns an
  object of class \code{FLXPmultinom} which extends class
  \code{FLXP} directly and is used for method dispatching.
}
\author{Friedrich Leisch and Bettina Gruen}
\keyword{models}
