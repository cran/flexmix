%
%  Copyright (C) 2005 Bettina Gruen
%  $Id: FLXdist-class.Rd 1679 2005-06-28 14:29:41Z gruen $
%
\name{FLXdist-class}
\docType{class}
\alias{FLXdist-class}
\alias{predict,FLXdist-method}
\alias{predict,FLXglmmodel-method}
\alias{parameters}
\alias{parameters,FLXdist-method}
\title{Class "FLXdist"}
\description{
  Objects of class \code{FLXdist} represent unfitted finite mixture
  models. 
}
\usage{
\S4method{parameters}{FLXdist}(object, component=1, model=1)
\S4method{predict}{FLXdist}(object, newdata=list())
}
\arguments{
  \item{object}{An object of class "FLXdist".}
  \item{component}{Number of component.}
  \item{model}{Number of model.}
  \item{newdata}{Dataframe containing new data.}
}

\section{Slots}{
  \describe{
    \item{model}{List of \code{FLXmodel} objects.}
    \item{prior}{Numeric vector with prior probabilities of clusters.}
    \item{components}{List describing the components using 
      \code{FLXcomponent} objects.}
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
  }
}


\author{Bettina Gruen}
\seealso{\code{FLXdist}}
\keyword{classes}

