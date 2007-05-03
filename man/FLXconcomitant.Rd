\name{FLXP}
\docType{class}
\alias{FLXPconstant}
\alias{FLXPmultinom}
\alias{FLXconstant}
\alias{FLXmultinom}
\alias{FLXP-class}
\alias{FLXPconstant-class}
\alias{FLXPmultinom-class}
\alias{show,FLXP-method}
\alias{refit,FLXP-method}
\alias{refit,FLXPmultinom-method}
\alias{summary,FLXRP-method}
\alias{summary,FLXRPmultinom-method}
\title{Creates the concomitant variable model.}
\description{
  Creator functions for the concomitant variable
  model. \code{FLXPconstant} specifies constant priors
  and \code{FLXPmultinom} multinomial logit models for the priors.  
}
\usage{
FLXPconstant()
FLXPmultinom(formula=~1)
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
\author{Bettina Gruen}
\keyword{models}
