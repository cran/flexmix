%
%  Copyright (C) 2004-2005 Friedrich Leisch
%  $Id: BIC-methods.Rd 3122 2006-11-05 15:10:25Z leisch $
%
\name{BIC-methods}
\docType{methods}
\title{Methods for Function BIC}
\alias{BIC,flexmix-method}
\alias{BIC,stepFlexmix-method}
\description{Compute the Bayesian Information Criterion.}
\section{Methods}{
  \describe{
    \item{object = flexmix:}{Compute the BIC of a \code{flexmix} object}
    \item{object = stepFlexmix:}{Compute the BIC of all
      models contained in the \code{stepFlexmix} object.}
  }
}
\keyword{methods}
