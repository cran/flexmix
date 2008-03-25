%
%  Copyright (C) 2004-2008 Friedrich Leisch and Bettina Gruen
%  $Id: FLXdist-class.Rd 3934 2008-03-26 08:36:26Z gruen $
%
\name{FLXdist-class}
\docType{class}
\alias{FLXdist-class}
\alias{predict,FLXdist-method}
\alias{predict,FLXM-method}
\alias{predict,FLXMRglm-method}
\alias{parameters}
\alias{parameters,FLXdist-method}
\alias{prior}
\alias{prior,FLXdist-method}
\title{Class "FLXdist"}
\description{
  Objects of class \code{FLXdist} represent unfitted finite mixture
  models. 
}
\usage{
\S4method{parameters}{FLXdist}(object, component=NULL, model=NULL, simplify=TRUE,
          drop=TRUE)
\S4method{predict}{FLXdist}(object, newdata=list(), aggregate=FALSE, ...)
}
\arguments{
  \item{object}{An object of class "FLXdist".}
  \item{component}{Number of component(s), if \code{NULL} all components
    are returned.}
  \item{model}{Number of model(s), if \code{NULL} all models are
    returned.}
  \item{simplify}{Logical, if \code{TRUE} the returned values are
    simplified to a vector or matrix if possible.}
  \item{drop}{Logical, if \code{TRUE} the function tries to simplify the
    return object by omitting lists of length one.}
  \item{newdata}{Dataframe containing new data.}
  \item{aggregate}{Logical, if \code{TRUE} then the predicted values for
    each model aggregated over the components are returned.}
  \item{\dots}{Passed to the method of the model class.}
}

\section{Slots}{
  \describe{
    \item{model}{List of \code{FLXM} objects.}
    \item{prior}{Numeric vector with prior probabilities of clusters.}
    \item{components}{List describing the components using 
      \code{FLXcomponent} objects.}
    \item{\code{concomitant}:}{Object of class \code{"FLXP"}.}
    \item{formula}{Object of class \code{"formula"}.}
    \item{call}{The function call used to create the object.}
    \item{k}{Number of clusters.}
  }  
}

\section{Accessor Functions}{
  The following functions should be used for accessing the corresponding
  slots:
  \describe{
    \item{\code{parameters}:}{The parameters for each model and component,
      return value depends on the model.}
    \item{\code{prior}:}{Numeric vector of prior class probabilities/component weights}

  }
}
\author{Friedrich Leisch and Bettina Gruen}
\seealso{\code{FLXdist}}
\keyword{classes}

