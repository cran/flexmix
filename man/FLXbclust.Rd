%
%  Copyright (C) 2004-2005 Friedrich Leisch
%  $Id: FLXbclust.Rd 1943 2005-12-21 10:36:16Z leisch $
%
\name{FLXbclust}
\alias{FLXbclust}
\title{FlexMix Binary Clustering Driver}
\description{
  This is a demo driver for \code{\link{flexmix}} implementing
  model-based clustering of binary data. 
}
\usage{
FLXbclust(formula = . ~ .)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{\link{flexmix}} using
    \code{\link{update.formula}}. Only the left-hand side (response) of
    the formula is used. Default is to use the original
    \code{\link{flexmix}} model
    formula.}
}
\details{
  This is mostly meant as a demo for FlexMix driver programming, but can
  be used to cluster binary data. The only parameter is the column-wise
  mean of the data, which equals the probability of observing a 1.
}
\value{
  \code{FLXbclust} returns an object of class \code{FLXmodel}.
}
\author{Friedrich Leisch}
\seealso{\code{\link{flexmix}}}
\keyword{cluster}
